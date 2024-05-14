#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>
 
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
 
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#include <cmath>

#define M_PI =3.14 ///to complete PI

///double sqr(double )
double I = 2E10;

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
Vector operator/(const Vector& a, const double b) {
    return Vector(a[0] / b, a[1] / b, a[2] / b);
}
double dot(const Vector& a, const Vector& b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
Vector cross(const Vector& a, const Vector& b) {
    return Vector(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]);
}
 

class Ray {
public:
    Ray(Vector &O, const Vector& u): O(O), u(u) {} )

    Vector O, u;
};


class Sphere {
public:
    Sphere(const Vector& C, double R), const :C(C)
    bool intersect(const& Ray r, Vector &P, Vector &N) const{
        double delta = sqrt(dot(r.u, r.O -C)) - ((r.O-C).norm2() - R*R);
        if (delta < 0) return false;

        double t1 = dot(u, C-r.O) - sqrt(delta);
        double t2 = dot(u, C-r.O) + sqrt(delta);

        if (t2<0) return false;
        double t;
        if (t1>0){
            t = t1;
        }
        else{
            t=t2
        }
    P =r.O +t*r.u
    N = P-C
    return true;
    }
    
};

class Scene{
    public:
        bool intersect(const Ray& r, Vector& P, Vector& N ) const{
            double bestt = 1E9;
            bool has_inter =false;
            for (int i=0; i<objects.size();i++){
                if (objects[i].intersect(r,Ptmp,Ntmp,t)){
                    if (t<bestt){
                        bestt =t;
                        P=Ptmp;
                        N=Ntmp;
                    }
                }
            }
        }
        std::vector<Sphere> objects;
    
}
 
int main() {
    int W = 512;
    int H = 512;

    
    double fov = 60*M_PI/180;
    double z= -W/(2*tan(fov/2.))
    Sphere S(Vector(0,0,0),10,Vector(1.,0.5,0.3));
    Vector C(0,0,55);
    Vector L

    std::vector<unsigned char> image(W * H * 3, 0);
    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {
            
            Vector u(j-W/2+0.5,H/2-i-0.5, z); 
            u.normalize();
            Ray r(C,u);
            Vector P;
            bool inter = S.intersect(r,P);
            Vector color(0,0,0);
            if (inter){
                Vector wlight = L-P;
                wlight.normalize();
                double l = I/(4*M_PI*(L-P).norm2())*dot(N,wlight);
                color = l*S.albedo/M_PI;
                color = Vector(255,255,255);
            }
            image[(i * W + j) * 3 + 0] = std::min(255,color[0]);
            image[(i * W + j) * 3 + 1] = std::min(255,color[1]);
            image[(i * W + j) * 3 + 2] = std::min(255,color[2]);
        }
    }
    stbi_write_png("image.png", W, H, 3, &image[0], 0);
 
    return 0;
}
