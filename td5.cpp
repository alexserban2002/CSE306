// for this code I collaborated with my colleague Cezara Petrui
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

void optimal_transport_match(std::vector<Vector> image, std::vector<Vector> model, int niter){
    int n = image.size();
    for (int iter = 0; iter < niter; iter++){
        std::vector<std::pair<double, int>> projection_image(n);
        std::vector<std::pair<double, int>> projection_model(n);
        double r1 = ((double)rand()) / RAND_MAX;
        double r2 = ((double)rand()) / RAND_MAX;
        Vector v(std::cos(2 * M_PI * r1) * std::sqrt(r2 * (1 - r2)), std::sin(2 * M_PI * r1) * std::sqrt(r2 * (1 - r2)),1 - 2 * r2);
        for (int i = 0; i < n; i++){
            projection_image[i] = std::pair<double, int>(dot(image[i], v), i);
            projection_model[i] = std::pair<double, int>(dot(model[i], v), i);
        }
        std::sort(projection_image.begin(), projection_image.end());
        std::sort(projection_model.begin(), projection_model.end());
        for (int i = 0; i < n; i++){
            image[projection_image[i].second] =  image[projection_image[i].second] + (projection_model[i].first - projection_image[i].first) * v;
        }
    }
};

//I asked ChatGPT to provide random values for the matrices I and M to see how the algorithm performs
int main(){ 
    std::vector<Vector> I = {
        Vector(0, 0, 0), Vector(255, 255, 255), Vector(128, 0, 0),
        Vector(0, 128, 0), Vector(0, 0, 128), Vector(128, 128, 128),
        Vector(192, 192, 192), Vector(64, 0, 64), Vector(0, 64, 64)
    };

    std::vector<Vector> M = {
        Vector(32, 32, 64), Vector(64, 64, 128), Vector(96, 96, 192),
        Vector(160, 0, 0), Vector(0, 160, 0), Vector(0, 0, 160),
        Vector(160, 160, 0), Vector(0, 160, 160), Vector(160, 0, 160)
    };

    int niter = 10;

    std::cout << "Original I:\n";
    for (const auto& vec : I) {
        std::cout << "[" << vec[0] << ", " << vec[1] << ", " << vec[2] << "]\n";
    }

    optimal_transport_match(I, M, niter);

    std::cout << "\nModified I after optimal transport:\n";
    for (const auto& vec : I) {
        std::cout << "[" << vec[0] << ", " << vec[1] << ", " << vec[2] << "]\n";
    }
    return 0;
}