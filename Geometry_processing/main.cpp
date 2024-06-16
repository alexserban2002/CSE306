// For this project I collaborated with my colleague Cezara Petrui

#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>
#include <iostream>
#include <random>
#include <list>
#include <chrono>
#include <thread>


#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define M_PI 3.14159265358979323846

#include "lbfgs.c"

class Vector {
public:
    explicit Vector(double x = 0, double y = 0, double z = 0) {
        data[0] = x;
        data[1] = y;
        data[2] = z;
    }
    double norm2() const {
        return data[0] * data[0] + data[1] * data[1] + data[2] * data[2];
    }
    double norm() const {
        return sqrt(norm2());
    }
    void normalize() {
        double n = norm();
        data[0] /= n;
        data[1] /= n;
        data[2] /= n;
    }
    double operator[](int i) const { return data[i]; };
    double& operator[](int i) { return data[i]; };
    double data[3];
};

Vector operator+(const Vector& a, const Vector& b) {
    return Vector(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
}
Vector operator-(const Vector& a, const Vector& b) {
    return Vector(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}
Vector operator*(const double a, const Vector& b) {
    return Vector(a*b[0], a*b[1], a*b[2]);
}
Vector operator*(const Vector& a, const double b) {
    return Vector(a[0]*b, a[1]*b, a[2]*b);
}
Vector operator*(const Vector& a, const Vector& b){
    return Vector(a[0]*b[0], a[1]*b[1], a[2]*b[2]);
}
Vector operator/(const Vector& a, const double b) {
    return Vector(a[0] / b, a[1] / b, a[2] / b);
}
double dot(const Vector& a, const Vector& b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
Vector cross(const Vector& a, const Vector& b) {
    return Vector(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]);
}

// if the Polygon class name conflicts with a class in wingdi.h on Windows, use a namespace or change the name
class Polygon {  
public:

    Polygon(){}

    Polygon(std::vector<Vector> vertices){
        this->vertices = vertices;
    }

    void add_vertex(Vector vertex){
        vertices.push_back(vertex);
    }

    double compute_area(){
        double area = 0;
        int n = vertices.size();
        for (int i = 0; i < n; i++){
            int j = (i + 1) % n;
            area += vertices[i][0] * vertices[j][1];
            area -= vertices[j][0] * vertices[i][1];
        } 
        return std::abs(area) / 2.0;
    }

    std::vector<Vector> vertices;
};  

Polygon clipByBisector(Polygon& polygon, Vector& P_i, Vector& P_j, double w_i, double w_j){
    // the code implements the Sutherland-Hodgman algorithm as presented in the lecture notes
    Polygon new_cell;
    int n = polygon.vertices.size();
    Vector M = (P_i + P_j) / 2;
    // M_prime is introduced for power diagram
    Vector M_prime = M + ((w_i - w_j) / (2 * (P_i - P_j).norm2())) * (P_j - P_i);
    for (int i = 0; i < n; i++){
        Vector A = polygon.vertices[i % n];
        Vector B = polygon.vertices[(i + 1) % n];
        double t = dot(M_prime - A, P_i - P_j) / dot(B - A, P_i - P_j);
        Vector P = A + t * (B - A);
        // if (dot(B - M_prime, P_j - P_i) < 0){
        //     if (dot(A - M_prime, P_j - P_i) >= 0){
        //         new_cell.add_vertex(P);
        //     }
        //     new_cell.add_vertex(B);
        // }
        // else if (dot(A - M_prime, P_j - P_i) < 0){
        //     new_cell.add_vertex(P);
        // }
        if (((B - P_i).norm2()) - w_i <= ((B - P_j).norm2()) - w_j){ //inside
            if (((A - P_i).norm2() - w_i) > ((A - P_j).norm2() - w_j)){ //outside
                new_cell.add_vertex(P);
            }
            new_cell.add_vertex(B);
        }
        else if (((A - P_i).norm2() - w_i) <= ((A - P_j).norm2() - w_j)){
            new_cell.add_vertex(P);
        }
    }
    return new_cell;
}

class Voronoi{
public:
    Voronoi(std::vector<Vector> points, std::vector<double> weights, std::vector<double> lambdas){
        this->points = points;
        this->weights = weights;
        this->lambdas = lambdas;
    }

    void add_point(Vector point){
        points.push_back(point);
    }

    void compute_voronoi(){
        // we just implemented Voronoï Parallel Linear Enumeration algorithm in 2D
        // as seen in the notes
        Polygon square;
        square.add_vertex(Vector(0, 0, 0));
        square.add_vertex(Vector(1, 0, 0));
        square.add_vertex(Vector(1, 1, 0));
        square.add_vertex(Vector(0, 1, 0)); // the introduction of the square is
                                            // similar to the on epresented during the TD
        int n = points.size();
        cells.resize(n);
        for (int i = 0; i < n; i++){
            Polygon cell = square;
            Vector curPoint = points[i];
            double curWeight = this->weights[i];
            for (int j = 0; j < n; j++){
                if (i == j) continue;
                Vector nextPoint = points[j];
                double nextWeight = this->weights[j];
                cell = clipByBisector(cell, curPoint, nextPoint, curWeight, nextWeight);
            }
            cells[i] = cell;
        }
    }

    std::vector<Vector> points;
    std::vector<Polygon> cells;
    std::vector<double> weights;
    std::vector<double> lambdas;
};

// for the evaluate fuunction I collaborated with my colleague Cezara Petrui
// since both of us were running into convergence issues initially
static lbfgsfloatval_t evaluate(void *instance, const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const int n, const lbfgsfloatval_t step){
    Voronoi *voronoi = static_cast<Voronoi*>(instance);
    lbfgsfloatval_t fx = 0.0;
    for (int i = 0; i < n; i++){
        voronoi->weights[i] = std::max(x[i], 0.); // we added 0 in the debugging process
                                                // to make sure that the weights don't go to 
                                                // - infinity, which was causing some initial
                                                // issues
    }
    voronoi->compute_voronoi();
    for (int i = 0; i < n; i++){
        // just used the formulae from the lecture notes for the integral
        Polygon cell = voronoi->cells[i];
        double t1 = 0;
        for (int j = 0; j < cell.vertices.size(); j++){
            int k = (j + 1) % (cell.vertices.size());
            double X_j = cell.vertices[j][0];
            double Y_j = cell.vertices[j][1];
            double X_k = cell.vertices[k][0];
            double Y_k = cell.vertices[k][1];
            t1 = t1 + ((X_j * Y_k - X_k * Y_j) * (std::pow(X_j, 2) + X_j * X_k + std::pow(X_k, 2) + std::pow(Y_j, 2) + Y_j * Y_k + std::pow(Y_k, 2) - 4 * (voronoi->points[i][0] * (X_j + X_k) + voronoi->points[i][1] * (Y_j + Y_k)) + 6 * voronoi->points[i].norm2()));
        }
        t1 /= 12;
        double t2 = -x[i] * cell.compute_area();
        double t3 = x[i] * voronoi->lambdas[i];
        fx += t1 + t2 + t3;
        g[i] = - voronoi->lambdas[i] + cell.compute_area(); // megative gradient
    }
    return -fx; // negative funtion for maximizing
}

static int progress(
    // the progress function was taken from the lbfgs library provided,
    // namely from https://github.com/chokkan/liblbfgs/tree/master/sample
    void *instance,
    const lbfgsfloatval_t *x,
    const lbfgsfloatval_t *g,
    const lbfgsfloatval_t fx,
    const lbfgsfloatval_t xnorm,
    const lbfgsfloatval_t gnorm,
    const lbfgsfloatval_t step,
    int n,
    int k,
    int ls
    )
{
    printf("Iteration %d:\n", k);
    printf("  fx = %f, x[0] = %f, x[1] = %f\n", fx, x[0], x[1]);
    printf("  xnorm = %f, gnorm = %f, step = %f\n", xnorm, gnorm, step);
    printf("\n");
    return 0;
}

 
// saves a static svg file. The polygon vertices are supposed to be in the range [0..1], and a canvas of size 1000x1000 is created
    void save_svg(const std::vector<Polygon> &polygons, std::string filename, std::string fillcol = "none") {
        FILE* f = fopen(filename.c_str(), "w+"); 
        fprintf(f, "<svg xmlns = \"http://www.w3.org/2000/svg\" width = \"1000\" height = \"1000\">\n");
        for (int i=0; i<polygons.size(); i++) {
            fprintf(f, "<g>\n");
            fprintf(f, "<polygon points = \""); 
            for (int j = 0; j < polygons[i].vertices.size(); j++) {
                fprintf(f, "%3.3f, %3.3f ", (polygons[i].vertices[j][0] * 1000), (1000 - polygons[i].vertices[j][1] * 1000));
            }
            fprintf(f, "\"\nfill = \"%s\" stroke = \"black\"/>\n", fillcol.c_str());
            fprintf(f, "</g>\n");
        }
        fprintf(f, "</svg>\n");
        fclose(f);
    }
 
 
// Adds one frame of an animated svg file. frameid is the frame number (between 0 and nbframes-1).
// polygons is a list of polygons, describing the current frame.
// The polygon vertices are supposed to be in the range [0..1], and a canvas of size 1000x1000 is created
    void save_svg_animated(const std::vector<Polygon> &polygons, std::string filename, int frameid, int nbframes) {
        FILE* f;
        if (frameid == 0) {
            f = fopen(filename.c_str(), "w+");
            fprintf(f, "<svg xmlns = \"http://www.w3.org/2000/svg\" width = \"1000\" height = \"1000\">\n");
            fprintf(f, "<g>\n");
        } else {
            f = fopen(filename.c_str(), "a+");
        }
        fprintf(f, "<g>\n");
        for (int i = 0; i < polygons.size(); i++) {
            fprintf(f, "<polygon points = \""); 
            for (int j = 0; j < polygons[i].vertices.size(); j++) {
                fprintf(f, "%3.3f, %3.3f ", (polygons[i].vertices[j][0] * 1000), (1000-polygons[i].vertices[j][1] * 1000));
            }
            fprintf(f, "\"\nfill = \"none\" stroke = \"black\"/>\n");
        }
        fprintf(f, "<animate\n");
        fprintf(f, "    id = \"frame%u\"\n", frameid);
        fprintf(f, "    attributeName = \"display\"\n");
        fprintf(f, "    values = \"");
        for (int j = 0; j < nbframes; j++) {
            if (frameid == j) {
                fprintf(f, "inline");
            } else {
                fprintf(f, "none");
            }
            fprintf(f, ";");
        }
        fprintf(f, "none\"\n    keyTimes = \"");
        for (int j = 0; j < nbframes; j++) {
            fprintf(f, "%2.3f", j / (double)(nbframes));
            fprintf(f, ";");
        }
        fprintf(f, "1\"\n   dur = \"5s\"\n");
        fprintf(f, "    begin = \"0s\"\n");
        fprintf(f, "    repeatCount = \"indefinite\"/>\n");
        fprintf(f, "</g>\n");
        if (frameid == nbframes - 1) {
            fprintf(f, "</g>\n");
            fprintf(f, "</svg>\n");
        }
        fclose(f);
    }

struct Facet {
    std::vector<Vector> vertices;
};

    void save_frame(const std::vector<Facet> &cells, std::string filename, int frameid = 0) {
        int W = 1000, H = 1000;
        std::vector<unsigned char> image(W*H * 3, 255);
        #pragma omp parallel for schedule(dynamic)
        for (int i = 0; i < cells.size(); i++) {
 
            double bminx = 1E9, bminy = 1E9, bmaxx = -1E9, bmaxy = -1E9;
            for (int j = 0; j < cells[i].vertices.size(); j++) {
                bminx = std::min(bminx, cells[i].vertices[j][0]);
                bminy = std::min(bminy, cells[i].vertices[j][1]);
                bmaxx = std::max(bmaxx, cells[i].vertices[j][0]);
                bmaxy = std::max(bmaxy, cells[i].vertices[j][1]);
            }
            bminx = std::min(W-1., std::max(0., W * bminx));
            bminy = std::min(H-1., std::max(0., H * bminy));
            bmaxx = std::max(W-1., std::max(0., W * bmaxx));
            bmaxy = std::max(H-1., std::max(0., H * bmaxy));
 
            for (int y = bminy; y < bmaxy; y++) {
                for (int x = bminx; x < bmaxx; x++) {
                    int prevSign = 0;
                    bool isInside = true;
                    double mindistEdge = 1E9;
                    for (int j = 0; j < cells[i].vertices.size(); j++) {
                        double x0 = cells[i].vertices[j][0] * W;
                        double y0 = cells[i].vertices[j][1] * H;
                        double x1 = cells[i].vertices[(j + 1) % cells[i].vertices.size()][0] * W;
                        double y1 = cells[i].vertices[(j + 1) % cells[i].vertices.size()][1] * H;
                        double det = (x - x0)*(y1-y0) - (y - y0)*(x1-x0);
                        //int sign = sgn(det);
                        int sign = (0 < det) - (det < 0);
                        if (prevSign == 0) prevSign = sign; else
                            if (sign == 0) sign = prevSign; else
                            if (sign != prevSign) {
                                isInside = false;
                                break;
                            }
                        prevSign = sign;
                        double edgeLen = sqrt((x1 - x0)*(x1 - x0) + (y1 - y0)*(y1 - y0));
                        double distEdge = std::abs(det)/ edgeLen;
                        double dotp = (x - x0)*(x1 - x0) + (y - y0)*(y1 - y0);
                        if (dotp<0 || dotp>edgeLen*edgeLen) distEdge = 1E9;
                        mindistEdge = std::min(mindistEdge, distEdge);
                    }
                    if (isInside) {
                        //if (i < N) {   // the N first particles may represent fluid, displayed in blue
                        //  image[((H - y - 1)*W + x) * 3] = 0;
                        //  image[((H - y - 1)*W + x) * 3 + 1] = 0;
                        //  image[((H - y - 1)*W + x) * 3 + 2] = 255;
                        //}
                        if (mindistEdge <= 2) {
                            image[((H - y - 1)*W + x) * 3] = 0;
                            image[((H - y - 1)*W + x) * 3 + 1] = 0;
                            image[((H - y - 1)*W + x) * 3 + 2] = 0;
                        }
 
                    }
                    
                }
            }
        }
        std::ostringstream os;
        os << filename << frameid << ".png";
        stbi_write_png(os.str().c_str(), W, H, 3, &image[0], 0);
    }

int main(){
    // the initialization of lbfgs was taken from the lbfgs library provided,
    // namely from https://github.com/chokkan/liblbfgs/tree/master/sample
    auto start_code = std::chrono::system_clock::now();
    int ret = 0;
    lbfgsfloatval_t fx;
    int num_points = 100;
    lbfgsfloatval_t *x = lbfgs_malloc(num_points);
    lbfgs_parameter_t param;
    Vector C = Vector(0.5, 0.5, 0);

    if (x == NULL) {
        printf("ERROR: Failed to allocate a memory block for variables.\n");
        return 1;
    }

    lbfgs_parameter_init(&param);

    std::vector<Vector> points;
    std::vector<double> weights;
    std::vector<double> lambdas;
    for (int i = 0; i < num_points; i++){
        double first_coordinate = ((double)rand()) / RAND_MAX;
        double second_coordinate = ((double)rand()) / RAND_MAX;
        Vector point(first_coordinate, second_coordinate, 0);
        points.push_back(point);
        double weight= ((double)rand()) / RAND_MAX;
        weights.push_back(weight);
        x[i] = weight;
        lambdas.push_back(1. / num_points);
    }
    Voronoi vor = Voronoi(points, weights, lambdas);

    ret = lbfgs(num_points, x, &fx, evaluate, progress, static_cast<void*>(&vor), &param);

    printf("L-BFGS optimization terminated with status code = %d\n", ret);
    printf("  fx = %f, x[0] = %f, x[1] = %f\n", fx, x[0], x[1]);


    std::cout<<vor.weights[0]<<" "<<vor.weights[1]<<std::endl;
    vor.compute_voronoi();
    auto end_code = std::chrono::system_clock::now();
    std::string filename = "lbfgs_diagram" + std::to_string(num_points) + ".svg";
    save_svg(vor.cells, filename);
    std::chrono::duration<double> elapsed_time = end_code - start_code;
    lbfgs_free(x);
    std::cout<<"Time elapsed for computing the power diagram for "<<num_points<<" points is "<<elapsed_time.count()<<std::endl;
    return 0;
}