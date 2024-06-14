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

Polygon clip_polygon(Polygon& subjectPolygon, Polygon& clipPolygon){
    int n = subjectPolygon.vertices.size();
    int m = clipPolygon.vertices.size();
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

Polygon clipByBisector(Polygon& polygon, Vector& P_i, Vector& P_j){
    Polygon* new_cell = new Polygon();
    int n = polygon.vertices.size();
    Vector M = (P_i + P_j) / 2;
    for (int i = 0; i < n; i++){
        Vector A = polygon.vertices[i % n];
        Vector B = polygon.vertices[(i + 1) % n];
        double t = dot(M - A, P_i - P_j) / dot(B - A, P_i - P_j);
        Vector P = A + t * (B - A);
        if (dot(B - M, P_j - P_i) < 0){
            if (dot(A - M, P_j - P_i) >= 0){
                new_cell->add_vertex(P);
            }
            new_cell->add_vertex(B);
        }
        else if (dot(A - M, P_j - P_i) < 0){
            new_cell->add_vertex(A);
        }
    }
    return *new_cell;
}

class Voronoi{
public:
    Voronoi(std::vector<Vector> points){
        this->points = points;
    }

    void add_point(Vector point){
        points.push_back(point);
    }

    void compute_voronoi(){
        Polygon square;
        square.add_vertex(Vector(0, 0, 0));
        square.add_vertex(Vector(1, 0, 0));
        square.add_vertex(Vector(1, 1, 0));
        square.add_vertex(Vector(0, 1, 0));
        int n = points.size();
        cells.resize(n);
        for (int i = 0; i < n; i++){
            Polygon cell = square;
            Vector curPoint = points[i];
            for (int j = 0; j < n; j++){
                if (i == j) continue;
                Vector nextPoint = points[j];
                Polygon cell = clipByBisector(cell, curPoint, nextPoint);
            }
            cells[i] = cell;
        }
    }

    std::vector<Vector> points;
    std::vector<Polygon> cells;
};

 
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

int main(){
    return 0;
}