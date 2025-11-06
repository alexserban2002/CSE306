// For the implementation of this algorithm I collaborated with my colleague Cezara Petrui
// For inspiration, I looked at the repository of Vrushank Agrawal: 
// https://github.com/vrushank-agrawal/CSE306/blob/main/Geometry%20Processing/main.cpp.
#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>
#include <iostream>
#include <random>
#include <list>
#include <chrono>
#include <thread>
#include <utility>
#define M_PI 3.14159265358979323846
#include <unordered_map>

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

std::vector<Vector> tutte_embedding(std::vector<Vector> mesh_vertices, std::vector<int> indices, std::unordered_map<Vector, std::vector<Vector>> map, int niter){
    int n = indices.size();
    std::vector<Vector> boundary_vertices;
    std::vector<Vector> new_vertices(n);
    for (int i = 0; i < n; i++){
        boundary_vertices.push_back(mesh_vertices[indices[i]]);
    }
    double s = 0;
    for (int i = 0; i < boundary_vertices.size(); i++){
        Vector b_i = boundary_vertices[i];
        Vector b_i_plus_1 = boundary_vertices[(i + 1) % boundary_vertices.size()];
        s += (b_i_plus_1 - b_i).norm();
    }
    double cs = 0;
    for (int i = 0; i < boundary_vertices.size(); i++){
        double theta_i = 2 * M_PI * cs / s;
        new_vertices[i] = Vector(std::cos(theta_i), std::sin(theta_i), 0);
        cs += (boundary_vertices[(i + 1) % boundary_vertices.size()] - boundary_vertices[i]).norm();
    }
    for (int i = 0; i < niter; i++){
        for (int j = 0; j < n; j++){
            if (std::find(indices.begin(), indices.end(), j) != indices.end()){
                continue;
            }
            Vector new_vertex(0, 0, 0);
            for (auto v: map[mesh_vertices[j]]){
                new_vertex = new_vertex + v;
            }
            new_vertices[j] = new_vertex / map[mesh_vertices[j]].size();
        }
    }
    return new_vertices;
}

int main(){
    return 0;
}