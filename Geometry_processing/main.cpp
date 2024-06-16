// For this project I collaborated with my colleague Cezara Petrui

#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>
#include <iostream>
#include <random>
#include <list>
#include <chrono>
#include <thread>
#include <utility>


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

const Vector g = Vector(0., -0.8, 0.);

class Edge {
    public:
        explicit Edge(Vector u, Vector v){
            edge = std::pair<Vector, Vector>{u, v};
            N = Vector(v[1] - u[1], u[0] - v[0], 0);
        }

        Vector intersect(Vector A, Vector B){
            Vector u = edge.first;
            Vector v = edge.second;
            double t = (dot(u - A, N) / dot(B - A, N));
            if (t < 0 || t > 1){
                throw std::runtime_error("No intersection exists!");
            }
            Vector P = A + t * (B - A);
            return P;
        }

        bool is_inside(Vector P){
            if (dot(P - edge.first, N) <= 0) return true;
            return false;
        }

    std::pair<Vector, Vector> edge;
    Vector N;
};


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

    double compute_signed_area(){
        if (vertices.size() <= 2) return 0.;
        double area = 0;
        int n = vertices.size();
        for (int i = 0; i < n; i++){
            int j = (i + 1) % n;
            area += vertices[i][0] * vertices[j][1];
            area -= vertices[j][0] * vertices[i][1];
        } 
        return area / 2.0;
    }

    double compute_area(){
        return std::abs(compute_signed_area());
    }

    Vector compute_centroid(){
        double signed_area = compute_signed_area();
        double centroid_x, centroid_y;
        int n = vertices.size();
        for (int i = 0; i < n; i++){
            int j = (i + 1) % n;
            centroid_x += (vertices[i][0] + vertices[j][0]) * (vertices[i][0] * vertices[j][1] - vertices[j][0] * vertices[i][1]);
            centroid_y += (vertices[i][1] + vertices[j][1]) * (vertices[i][0] * vertices[j][1] - vertices[j][0] * vertices[i][1]);
        }
        centroid_x = (1. / (6. * signed_area)) * centroid_x;
        centroid_y = (1. / (6. * signed_area)) * centroid_y;
        return Vector(centroid_x, centroid_y, 0.);
    }

    void compute_edges(){
        int n = vertices.size(); 
        for (int i = 0; i < n - 2; i++){
            edges.push_back(Edge(vertices[i], vertices[i + 1]));
        }
        edges.push_back(Edge(vertices[n - 1], vertices[0]));
    }

    std::vector<Edge> edges;
    std::vector<Vector> vertices;
};  

Polygon clipByPolygon(Polygon& subjectPolygon, Polygon& clipPolygon){
    int n = subjectPolygon.vertices.size();
    int m = clipPolygon.vertices.size();
    clipPolygon.compute_edges();
    for (auto& clipEdge: clipPolygon.edges){
        Polygon* outPolygon = new Polygon();
        for (int j = 0; j < n; j++){
            Vector curVertex = subjectPolygon.vertices[j];
            int k = n - 1;
            if (j > 0) k = j - 1;
            Vector prevVertex = subjectPolygon.vertices[k];
            Vector intersection = clipEdge.intersect(prevVertex, curVertex);
            if (clipEdge.is_inside(curVertex)){
                if (! clipEdge.is_inside(prevVertex)){
                    outPolygon->add_vertex(intersection);
                }
                outPolygon->add_vertex(curVertex);
            }
            else if (clipEdge.is_inside(prevVertex)){
                outPolygon->add_vertex(intersection);
            }
        }
        subjectPolygon = *outPolygon;
    }
    return subjectPolygon;
}

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
        else if (((A - P_i).norm2() - w_i) <= ((A - P_j).norm2() - w_j)){ //inside
            new_cell.add_vertex(P);
        }
    }
    return new_cell;
}

Polygon create_disk(Vector& center, double radius, int num_points = 30){
    std::vector<Vector> vertices_disk(num_points);
    for (int i = 0; i < num_points; i++){
        Vector vertex(std::cos(2 * M_PI * i / num_points), std::sin(2 * M_PI * i / num_points), 0.);
        vertices_disk[i] = center + radius * vertex;
    }
    return Polygon(vertices_disk);
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

class VoronoiFluid{
public:
    VoronoiFluid(std::vector<Vector> particles, std::vector<double> weights, std::vector<double> lambdas, double desired_fluid_volume){
        this->particles = particles;
        int n = particles.size();
        weight_air = weights[n];
        weights.pop_back();
        weights_fluids = weights;
        lambda_air = lambdas[n];
        lambdas.pop_back();
        lambdas_fluids = lambdas;
        this->desired_fluid_volume = desired_fluid_volume;
    }

    void add_particle(Vector particle){
        particles.push_back(particle);
    }

    void compute_voronoi_fluids(){
        Polygon square;
        square.add_vertex(Vector(0, 0, 0));
        square.add_vertex(Vector(1, 0, 0));
        square.add_vertex(Vector(1, 1, 0));
        square.add_vertex(Vector(0, 1, 0));
        int n = particles.size();
        cells.resize(n);
        for (int i = 0; i < n; i++){
            Polygon cell = square;
            Vector curPoint = particles[i];
            double curWeight = weights_fluids[i];
            for (int j = 0; j < n; j++){
                if (i == j) continue;
                Vector nextPoint = particles[j];
                double nextWeight = weights_fluids[j];
                cell = clipByBisector(cell, curPoint, nextPoint, curWeight, nextWeight);
            }
            if (weights_fluids[i] <= weight_air){
                cells[i] = Polygon();
            }
            else{
                double R = std::sqrt(weights_fluids[i] - weight_air);
                Polygon clippind_disk = create_disk(curPoint, R);
                cell = clipByPolygon(cell, clippind_disk);
                cells[i] = cell;
            }   
        }
    }

    std::vector<Vector> particles;
    std::vector<Polygon> cells;
    std::vector<double> weights_fluids;
    std::vector<double> lambdas_fluids;
    double weight_air;
    double lambda_air;
    double desired_fluid_volume;
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
        t1 /= 12.;
        double t2 = -x[i] * cell.compute_area();
        double t3 = x[i] * voronoi->lambdas[i];
        fx += t1 + t2 + t3;
        g[i] = - voronoi->lambdas[i] + cell.compute_area(); // megative gradient
    }
    return -fx; // negative funtion for maximizing
}

static lbfgsfloatval_t evaluate_fluid(void *instance, const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const int n, const lbfgsfloatval_t step){
    VoronoiFluid *vor_fluid = static_cast<VoronoiFluid*>(instance);
    lbfgsfloatval_t fx = 0.0;
    for (int i = 0; i < n - 1; i++){
        vor_fluid->weights_fluids[i] = x[i];
    }
    vor_fluid->weight_air = x[n - 1];
    vor_fluid->compute_voronoi_fluids();
    double total_area = 0.0;
    for (int i = 0; i < n - 1; i++){
        Polygon cell = vor_fluid->cells[i];
        double t1 = 0;
        for (int j = 0; j < cell.vertices.size(); j++){
            int k = (j + 1) % (cell.vertices.size());
            double X_j = cell.vertices[j][0];
            double Y_j = cell.vertices[j][1];
            double X_k = cell.vertices[k][0];
            double Y_k = cell.vertices[k][1];
            t1 = t1 + ((X_j * Y_k - X_k * Y_j) * (std::pow(X_j, 2) + X_j * X_k + std::pow(X_k, 2) + std::pow(Y_j, 2) + Y_j * Y_k + std::pow(Y_k, 2) - 4 * (vor_fluid->particles[i][0] * (X_j + X_k) + vor_fluid->particles[i][1] * (Y_j + Y_k)) + 6 * vor_fluid->particles[i].norm2()));
        }
        t1 /= 12.;
        double area = cell.compute_area();
        total_area += area;
        double t2 = -x[i] * area;
        double t3 = vor_fluid->desired_fluid_volume / ((double) n - 1) * x[i];
        fx += t1 + t2 + t3;
        g[i] = vor_fluid->desired_fluid_volume / ((double) n - 1) - area;
    }
    double desired_air_volume = 1. - vor_fluid->desired_fluid_volume;
    double estimated_air_volume = 1 - total_area;
    fx += x[n - 1] * (desired_air_volume - estimated_air_volume);
    g[n - 1] = desired_air_volume / ((double) n - 1) - estimated_air_volume;
    return fx; 
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

    void save_frame(const std::vector<Polygon> &cells, std::string filename, int frameid = 0) {
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

std::pair<std::vector<Vector>, std::vector<Vector>> gallouet_merigot_scheme(VoronoiFluid vor_fluid, std::vector<Vector> v, std::vector<double> m, double eps = 0.004, double dt = 0.002){
    int n = vor_fluid.particles.size();
    vor_fluid.compute_voronoi_fluids();
    std::vector<Vector>new_positions(n);
    std::vector<Vector>new_velocities(n);
    for (int i = 0; i < n; i++){
        Polygon cell = vor_fluid.cells[i];
        Vector centroid_cell = cell.compute_centroid();
        Vector F_spring_particle = (1. / (eps * eps)) * (centroid_cell - vor_fluid.particles[i]);
        Vector F_particle = F_spring_particle + m[i] * g;
        Vector new_v_particle = v[i] + dt/m[i] * F_particle;
        Vector new_position_particle = vor_fluid.particles[i] + dt * v[i];
        new_positions[i] = new_position_particle;
        new_velocities[i] = new_v_particle;
    }
    return std::pair<std::vector<Vector>, std::vector<Vector>>{new_positions, new_velocities};
}

void compute_frame(std::vector<Vector> fluid_particles, std::vector<double> weights, std::vector<double> lambdas, double desired_fluid_volume, std::vector<Vector> v, std::vector<double> m, int frameid){
    int num_particles = fluid_particles.size();
    VoronoiFluid vor_fluid = VoronoiFluid(fluid_particles, weights, lambdas, desired_fluid_volume);
    int ret = 0;
    lbfgsfloatval_t fx;
    lbfgsfloatval_t *x = lbfgs_malloc(num_particles + 1);
    lbfgs_parameter_t param;

    if (x == NULL) {
        printf("ERROR: Failed to allocate a memory block for variables.\n");
    }

    lbfgs_parameter_init(&param);

    ret = lbfgs(num_particles + 1, x, &fx, evaluate_fluid, progress, static_cast<void*>(&vor_fluid), &param);

    printf("L-BFGS optimization terminated with status code = %d\n", ret);
    printf("  fx = %f, x[0] = %f, x[1] = %f\n", fx, x[0], x[1]);

    lbfgs_free(x);

    std::pair<std::vector<Vector>, std::vector<Vector>> outputs = gallouet_merigot_scheme(vor_fluid, v, m);
    v = outputs.second;
    fluid_particles = outputs.first;
    vor_fluid.particles = outputs.first;
    vor_fluid.compute_voronoi_fluids();
    std::string filename = "fluids_frame" + std::to_string(frameid) + ".svg";
    save_frame(vor_fluid.cells, filename, frameid);
}

int main(){
    // the initialization of lbfgs was taken from the lbfgs library provided,
    // namely from https://github.com/chokkan/liblbfgs/tree/master/sample
    int test_lbfgs = 0; //for testing the power diagram and the LBFGS optimization, make this 1
    auto start_code = std::chrono::system_clock::now();
    if (test_lbfgs){
        int num_points = 100;
        int ret = 0;
        lbfgsfloatval_t fx;
        lbfgsfloatval_t *x = lbfgs_malloc(num_points);
        lbfgs_parameter_t param;

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
        
        std::string filename = "lbfgs_diagram" + std::to_string(num_points) + ".svg";
        save_svg(vor.cells, filename);
        
        lbfgs_free(x);
        
    }
    else{ // handle the fluids
        int num_particles = 700;
        double desired_fluid_volume = 0.25;
        int num_frames = 10;
       // int num_fluid_particles = 700;
        double desired_air_volume = 1 - desired_fluid_volume;
        std::vector<Vector> fluid_particles(num_particles);
        std::vector<double> weights(num_particles + 1);
        std::vector<double> lambdas(num_particles + 1);
        for (int i = 0; i < num_particles; i++){
            double first_coordinate_fluid = ((double)rand()) / RAND_MAX;
            double second_coordinate_fluid = ((double)rand()) / RAND_MAX;
            Vector fluid_particle(first_coordinate_fluid, second_coordinate_fluid, 0);
            fluid_particles[i] = fluid_particle;
            double fluid_weight = ((double)rand()) / RAND_MAX;
            weights[i] = fluid_weight;
            lambdas[i] = desired_fluid_volume / num_particles;
        }
        //fluid_particles[num_particles] = Vector(((double)rand()) / RAND_MAX, ((double)rand()) / RAND_MAX, 0);
        weights[num_particles] = ((double)rand()) / RAND_MAX;
        lambdas[num_particles] = desired_air_volume;

        std::vector<Vector> v(num_particles);
        std::vector<double> m(num_particles);
        for (int i = 0; i < num_particles; i++){
            v[i] = Vector(0., 0., 0.);
            m[i] = 200.;
        }

        for (int frame = 0; frame < num_frames; frame++){
            compute_frame(fluid_particles, weights, lambdas, desired_fluid_volume, v, m, frame);
        }
    }
    auto end_code = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_time = end_code - start_code;
    //std::cout<<"Time elapsed for computing the power diagram for "<<num_points<<" points is "<<elapsed_time.count()<<std::endl;
    return 0;
}