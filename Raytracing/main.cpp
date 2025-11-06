#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>
#include <iostream>
#include <random>
#include <list>
#include <chrono>
#include <thread>

using namespace std;

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

static std :: default_random_engine engine(10) ; 
static std::uniform_real_distribution<double> uniform(0, 1);

void delayOneSecond() {
    std::this_thread::sleep_for(std::chrono::seconds(5));
}

class Ray {
public:
    Ray(const Vector &O, const Vector& u): O(O), u(u) {};
    Vector O, u;
};


class Intersection {
public:
    Intersection(const Vector& albedo = Vector(0., 0., 0.),  bool reflection = false, bool refraction = false, double refractive_index = 1.0, bool hit = false, const Vector& point = Vector(), const Vector& normal = Vector(), double distance = 0.)
        : albedo(albedo), reflection(reflection), refraction(refraction), refractive_index(refractive_index), hit(hit), point(point), normal(normal), distance(distance)  {}

    bool hit;
    Vector point;
    Vector normal;
    Vector albedo;
    double distance;
    double refractive_index;
    bool reflection;
    bool refraction;
};

class Geometry{
    public:
        virtual Intersection intersect(const Ray& ray) = 0;
};

class BoundingBox{
    public:
        Vector B_min;
        Vector B_max;

        BoundingBox(Vector B_min = Vector(), Vector B_max = Vector()): B_min(B_min), B_max(B_max){};
};

struct Node{
        BoundingBox boundingbox;
        int starting_triangle;
        int ending_triangle;
        Node* left;
        Node* right;
};


class TriangleIndices {
public:
	TriangleIndices(int vtxi = -1, int vtxj = -1, int vtxk = -1, int ni = -1, int nj = -1, int nk = -1, int uvi = -1, int uvj = -1, int uvk = -1, int group = -1, bool added = false) : vtxi(vtxi), vtxj(vtxj), vtxk(vtxk), uvi(uvi), uvj(uvj), uvk(uvk), ni(ni), nj(nj), nk(nk), group(group) {
	};
	int vtxi, vtxj, vtxk; // indices within the vertex coordinates array
	int uvi, uvj, uvk;  // indices within the uv coordinates array
	int ni, nj, nk;  // indices within the normals array
	int group;       // face group
};


class TriangleMesh: public Geometry {
private:
    Vector translation;
    double scaling_factor;
    Vector albedo; 

public:
    std::vector<TriangleIndices> indices;
	std::vector<Vector> vertices;
	std::vector<Vector> normals;
	std::vector<Vector> uvs;
	std::vector<Vector> vertexcolors;
	Node* root;

    ~TriangleMesh() {}
	TriangleMesh(const Vector& translation, double scaling_factor, const Vector& albedo){// :translation(translation), scaling_factor(scaling_factor), albedo(albedo) {
        this->translation = translation;
        this->scaling_factor = scaling_factor;
        this->albedo = albedo;
        this->root = new Node();
    }; 
	
	void readOBJ(const char* obj) {

		char matfile[255];
		char grp[255];

		FILE* f;
		f = fopen(obj, "r");
		int curGroup = -1;
		while (!feof(f)) {
			char line[255];
			if (!fgets(line, 255, f)) break;

			std::string linetrim(line);
			linetrim.erase(linetrim.find_last_not_of(" \r\t") + 1);
			strcpy(line, linetrim.c_str());

			if (line[0] == 'u' && line[1] == 's') {
				sscanf(line, "usemtl %[^\n]\n", grp);
				curGroup++;
			}

			if (line[0] == 'v' && line[1] == ' ') {
				Vector vec;

				Vector col;
				if (sscanf(line, "v %lf %lf %lf %lf %lf %lf\n", &vec[0], &vec[1], &vec[2], &col[0], &col[1], &col[2]) == 6) {
					col[0] = std::min(1., std::max(0., col[0]));
					col[1] = std::min(1., std::max(0., col[1]));
					col[2] = std::min(1., std::max(0., col[2]));

					vertices.push_back(vec);
					vertexcolors.push_back(col);

				} else {
					sscanf(line, "v %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
					vertices.push_back(vec);
				}
			}
			if (line[0] == 'v' && line[1] == 'n') {
				Vector vec;
				sscanf(line, "vn %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
				normals.push_back(vec);
			}
			if (line[0] == 'v' && line[1] == 't') {
				Vector vec;
				sscanf(line, "vt %lf %lf\n", &vec[0], &vec[1]);
				uvs.push_back(vec);
			}
			if (line[0] == 'f') {
				TriangleIndices t;
				int i0, i1, i2, i3;
				int j0, j1, j2, j3;
				int k0, k1, k2, k3;
				int nn;
				t.group = curGroup;

				char* consumedline = line + 1;
				int offset;

				nn = sscanf(consumedline, "%u/%u/%u %u/%u/%u %u/%u/%u%n", &i0, &j0, &k0, &i1, &j1, &k1, &i2, &j2, &k2, &offset);
				if (nn == 9) {
					if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
					if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
					if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
					if (j0 < 0) t.uvi = uvs.size() + j0; else	t.uvi = j0 - 1;
					if (j1 < 0) t.uvj = uvs.size() + j1; else	t.uvj = j1 - 1;
					if (j2 < 0) t.uvk = uvs.size() + j2; else	t.uvk = j2 - 1;
					if (k0 < 0) t.ni = normals.size() + k0; else	t.ni = k0 - 1;
					if (k1 < 0) t.nj = normals.size() + k1; else	t.nj = k1 - 1;
					if (k2 < 0) t.nk = normals.size() + k2; else	t.nk = k2 - 1;
					indices.push_back(t);
				} else {
					nn = sscanf(consumedline, "%u/%u %u/%u %u/%u%n", &i0, &j0, &i1, &j1, &i2, &j2, &offset);
					if (nn == 6) {
						if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
						if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
						if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
						if (j0 < 0) t.uvi = uvs.size() + j0; else	t.uvi = j0 - 1;
						if (j1 < 0) t.uvj = uvs.size() + j1; else	t.uvj = j1 - 1;
						if (j2 < 0) t.uvk = uvs.size() + j2; else	t.uvk = j2 - 1;
						indices.push_back(t);
					} else {
						nn = sscanf(consumedline, "%u %u %u%n", &i0, &i1, &i2, &offset);
						if (nn == 3) {
							if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
							if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
							if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
							indices.push_back(t);
						} else {
							nn = sscanf(consumedline, "%u//%u %u//%u %u//%u%n", &i0, &k0, &i1, &k1, &i2, &k2, &offset);
							if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
							if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
							if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
							if (k0 < 0) t.ni = normals.size() + k0; else	t.ni = k0 - 1;
							if (k1 < 0) t.nj = normals.size() + k1; else	t.nj = k1 - 1;
							if (k2 < 0) t.nk = normals.size() + k2; else	t.nk = k2 - 1;
							indices.push_back(t);
						}
					}
				}

				consumedline = consumedline + offset;

				while (true) {
					if (consumedline[0] == '\n') break;
					if (consumedline[0] == '\0') break;
					nn = sscanf(consumedline, "%u/%u/%u%n", &i3, &j3, &k3, &offset);
					TriangleIndices t2;
					t2.group = curGroup;
					if (nn == 3) {
						if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
						if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
						if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
						if (j0 < 0) t2.uvi = uvs.size() + j0; else	t2.uvi = j0 - 1;
						if (j2 < 0) t2.uvj = uvs.size() + j2; else	t2.uvj = j2 - 1;
						if (j3 < 0) t2.uvk = uvs.size() + j3; else	t2.uvk = j3 - 1;
						if (k0 < 0) t2.ni = normals.size() + k0; else	t2.ni = k0 - 1;
						if (k2 < 0) t2.nj = normals.size() + k2; else	t2.nj = k2 - 1;
						if (k3 < 0) t2.nk = normals.size() + k3; else	t2.nk = k3 - 1;
						indices.push_back(t2);
						consumedline = consumedline + offset;
						i2 = i3;
						j2 = j3;
						k2 = k3;
					} else {
						nn = sscanf(consumedline, "%u/%u%n", &i3, &j3, &offset);
						if (nn == 2) {
							if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
							if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
							if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
							if (j0 < 0) t2.uvi = uvs.size() + j0; else	t2.uvi = j0 - 1;
							if (j2 < 0) t2.uvj = uvs.size() + j2; else	t2.uvj = j2 - 1;
							if (j3 < 0) t2.uvk = uvs.size() + j3; else	t2.uvk = j3 - 1;
							consumedline = consumedline + offset;
							i2 = i3;
							j2 = j3;
							indices.push_back(t2);
						} else {
							nn = sscanf(consumedline, "%u//%u%n", &i3, &k3, &offset);
							if (nn == 2) {
								if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
								if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
								if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
								if (k0 < 0) t2.ni = normals.size() + k0; else	t2.ni = k0 - 1;
								if (k2 < 0) t2.nj = normals.size() + k2; else	t2.nj = k2 - 1;
								if (k3 < 0) t2.nk = normals.size() + k3; else	t2.nk = k3 - 1;								
								consumedline = consumedline + offset;
								i2 = i3;
								k2 = k3;
								indices.push_back(t2);
							} else {
								nn = sscanf(consumedline, "%u%n", &i3, &offset);
								if (nn == 1) {
									if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
									if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
									if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
									consumedline = consumedline + offset;
									i2 = i3;
									indices.push_back(t2);
								} else {
									consumedline = consumedline + 1;
								}
							}
						}
					}
				}

			}

		}
		fclose(f);
	}

    BoundingBox boundingBox(int starting_triangle, int ending_triangle){
        double min_x = std::numeric_limits<double>::infinity();
        double max_x = - min_x;
        double min_y = std::numeric_limits<double>::infinity();
        double max_y = - min_y;
        double min_z = std::numeric_limits<double>::infinity();
        double max_z = - min_z;

        for (int i = starting_triangle; i < ending_triangle; i++){
            Vector v1 = scaling_factor * vertices[indices[i].vtxi] + translation;
            Vector v2 = scaling_factor * vertices[indices[i].vtxj] + translation;
            Vector v3 = scaling_factor * vertices[indices[i].vtxk] + translation;
            max_x = std::max(std::max(std::max(v1[0], v2[0]), v3[0]), max_x);
            max_y = std::max(std::max(std::max(v1[1], v2[1]), v3[1]), max_y);
            max_z = std::max(std::max(std::max(v1[2], v2[2]), v3[2]), max_z);
            min_x = std::min(std::min(std::min(v1[0], v2[0]), v3[0]), min_x);
            min_y = std::min(std::min(std::min(v1[1], v2[1]), v3[1]), min_y);
            min_z = std::min(std::min(std::min(v1[2], v2[2]), v3[2]), min_z);
        }
        Vector B_min = Vector(min_x, min_y, min_z);
        Vector B_max = Vector(max_x, max_y, max_z);
        return BoundingBox(B_min, B_max);
    }


    bool intersect_bounding_box(const Ray &r, BoundingBox boundingbox, double* inter_distance){
        double t0_x = (boundingbox.B_min[0] - r.O[0]) / (r.u[0]); 
        double t1_x = (boundingbox.B_max[0] - r.O[0]) / (r.u[0]); 

        double t0_y = (boundingbox.B_min[1] - r.O[1]) / (r.u[1]); 
        double t1_y = (boundingbox.B_max[1] - r.O[1]) / (r.u[1]); 

        double t0_z = (boundingbox.B_min[2] - r.O[2]) / (r.u[2]); 
        double t1_z = (boundingbox.B_max[2] - r.O[2]) / (r.u[2]); 

        double maximum = std::max(std::max(std::min(t0_x, t1_x), std::min(t0_y, t1_y)), std::min(t0_z, t1_z));
        double minimum = std::min(std::min(std::max(t0_x, t1_x), std::max(t0_y, t1_y)), std::max(t0_z, t1_z));

        if(maximum < minimum){
            *inter_distance = maximum;
            return true;     
        }    
        return false;

    }

    Vector compute_barycenter(int i){
        Vector A = scaling_factor * vertices[indices[i].vtxi] + translation;
        Vector B = scaling_factor * vertices[indices[i].vtxj] + translation;
        Vector C = scaling_factor * vertices[indices[i].vtxk] + translation;
        return Vector((A + B + C) / 3.);
    }

    void BVH(Node* node, int starting_triangle, int ending_triangle){
        node->boundingbox = boundingBox(starting_triangle, ending_triangle);
        node->starting_triangle = starting_triangle;
        node->ending_triangle = ending_triangle;
     
        Vector diag = node->boundingbox.B_max - node->boundingbox.B_min;
        Vector middle_diag = node->boundingbox.B_min + 0.5 * diag;
        int longest_axis = 0;
        double longest_distance = diag[0];
        if (diag[1] > longest_distance) {
            longest_distance = diag[1];
            longest_axis = 1;
        }
        if (diag[2] > longest_distance){
            longest_distance = diag[2];
            longest_axis = 2;
        }

        int pivot_index = starting_triangle;
        for (int i = starting_triangle; i < ending_triangle; i++){
            Vector barycenter = compute_barycenter(i);
            if (barycenter[longest_axis] < middle_diag[longest_axis]){
                std::swap(indices[i], indices[pivot_index]);
                pivot_index++;
            }
        }

        if (pivot_index <= starting_triangle || pivot_index >= ending_triangle - 1 || ending_triangle - starting_triangle < 5){
            return;
        }

        node->left = new Node();
        node->right = new Node();
        this->BVH(node->left, starting_triangle, pivot_index);
        this->BVH(node->right, pivot_index, ending_triangle);
    }

    Intersection intersect(const Ray &r) override{
        Intersection intersection;
        intersection.albedo = this->albedo;
        double inter_distance;

        if (!intersect_bounding_box(r, this->root->boundingbox, &inter_distance)){ return Intersection();}
        std::list<Node*> nodes_to_visit;
        nodes_to_visit.push_front(root);
        double best_inter_distance = std::numeric_limits<double>::infinity();
        while (!nodes_to_visit.empty()){
            Node* curNode = nodes_to_visit.back();
            nodes_to_visit.pop_back();
            if (curNode->left){
                if (intersect_bounding_box(r, curNode->left->boundingbox, &inter_distance)){
                   if (inter_distance < best_inter_distance){
                        nodes_to_visit.push_back(curNode->left);
                    }
                }
                if (intersect_bounding_box(r, curNode->right->boundingbox, &inter_distance)){
                    if (inter_distance < best_inter_distance){
                        nodes_to_visit.push_back(curNode->right);
                    }
                }

            }
                else {
                    for (int i = curNode->starting_triangle; i < curNode->ending_triangle; i++){
                        Vector A = scaling_factor * vertices[indices[i].vtxi] + translation;
                        Vector B = scaling_factor * vertices[indices[i].vtxj] + translation;
                        Vector C = scaling_factor * vertices[indices[i].vtxk] + translation;
                        Vector e1 = B - A;
                        Vector e2 = C - A;
                        Vector N = cross(e1, e2);
                        Vector dst_to_center = A - r.O;
                        double dot_product = dot(r.u, N);
                        double beta = dot(e2, cross(dst_to_center, r.u)) / dot_product;
                        double gamma = (-1.) * dot(e1, cross(dst_to_center, r.u)) / dot_product;
                        double alpha = 1.0 - beta - gamma;
                        double t = dot(dst_to_center, N) / dot_product;
                        if (alpha >= 0. && beta >= 0. && gamma>= 0. && t < best_inter_distance && t >= 0.){
                            best_inter_distance = t;
                            intersection.distance = t;
                            intersection.hit = true;
                            intersection.point = A + beta * e1 + gamma * e2;
                            intersection.normal = N;
                            // Vector N_A = normals[indices[i].vtxi] + translation;
                            // Vector N_B = normals[indices[i].vtxj] + translation;
                            // Vector N_C = normals[indices[i].vtxk] + translation;

                            // Vector N_phong = alpha * N_A + beta * N_B + gamma * N_C;
                            // N_phong.normalize();
                            // intersection.normal_phong = N_phong;
                        }
                    }
                }
        }
        return intersection;
    }
};


class Sphere: public Geometry {
private:
    Vector C;
    double R;
    Vector albedo; 
    double refractive_index;
    bool is_mirror;
    bool is_transparent;
    bool is_inside;
    
public:
    Sphere(const Vector& C, double R, const Vector& albedo, double refractive_index = 1., bool is_mirror = false, bool is_transparent = false, bool is_inside = false): 
        C(C), R(R), albedo(albedo), refractive_index(refractive_index), is_mirror(is_mirror), is_transparent(is_transparent), is_inside(is_inside) {};
    
    Intersection intersect(const Ray& r) override{
        Intersection intersection(this->albedo, this->is_mirror, this->is_transparent, this->refractive_index);
        Vector distance = r.O - C;
        double dot_product = dot(r.u, distance);
        double delta = std::pow(dot_product, 2) - (distance.norm2() - R * R);
        if (delta < 0){
            intersection.hit = false;
            return intersection;
        }
        double t1 = - dot_product - sqrt(delta);
        double t2 = - dot_product + sqrt(delta);
            if (t2 < 0) {
                intersection.hit = false;
                return intersection;
            }
            double t = t2;
            if (t1 > 0){
                t = t1;
            }
            Vector P = r.O + t * r.u;
            Vector N = P - C;
            N.normalize();
            intersection.hit = true;
            intersection.distance = t;
            intersection.point = P;
            intersection.normal = N;
            if (this->is_inside){
                intersection.normal = (-1.) * intersection.normal;
            }
            return intersection;
    }

};

class Scene {
public:
    Scene(const Vector& light_source) : light_source(light_source) {};
    void addGeometry(Geometry* geometry) {
        objects.push_back(geometry);
    }
    Intersection intersect(const Ray& r) const { 
        Intersection closestIntersection;
        double closestT = std::numeric_limits<double>::infinity();
        for (const auto& object : objects) {
            Intersection intersection = object->intersect(r);
            if (intersection.hit) {
                double t = intersection.distance; 
                if (t < closestT) {
                    closestT = t;
                    closestIntersection = intersection;
                }
            }
        }
        return closestIntersection;
    }

    Vector random_cos(const Vector &N){
        int min_component = 3;
        double min_component_value = std::numeric_limits<double>::infinity();
        for (int i = 0; i < 3; i++){
            if (std::abs(N[i]) < min_component_value){
                min_component_value = std::abs(N[i]);
                min_component = i;
            }
        }

        Vector T1(0.,0.,0.);
        if (min_component == 0){
            T1 = Vector(0., (-1) * N[2], N[1]);
        }
        else if(min_component == 1){
            T1 = Vector((-1) * N[2], 0., N[0]);
        }
        else{
            T1 = Vector((-1) * N[1], N[0], 0.);
        }
        T1.normalize();

        Vector T2 = cross(N, T1);

        double r1 = uniform(engine);
        double r2 = uniform(engine);
        double x = std::cos(2 * M_PI * r1) * std::sqrt(1 - r2);
        double y = std::sin(2 * M_PI * r1) * std::sqrt(1 - r2);
        double z = std::sqrt(r2);

        Vector V = x * T1 + y * T2 + z * N;
        return V;
    }

    Vector getColor(const Ray& r, int ray_depth){
        if(ray_depth < 0) {return Vector(0.,0.,0.);}

        Intersection intersection = this->intersect(r);
        Vector P = intersection.point;
        Vector N = intersection.normal;
        double eps = 1e-6; 
        Vector Lo(0.,0.,0.);
        double dot_product = dot(r.u, N);
        if (intersection.hit){
            if (intersection.reflection){
                P = P + eps * N;
                Ray reflected_ray(P, r.u - 2. * dot_product * N);
                return getColor(reflected_ray, ray_depth - 1);
            }
            else if(intersection.refraction && intersection.refractive_index != 1){
                double refractive_index_object = intersection.refractive_index;
                double refraction_index = 1./refractive_index_object;
                Vector tmp_N = N;
                if (dot_product > 0){
                    tmp_N = (-1.0) * N;
                    refraction_index = 1./refraction_index;
                }
                P = P - eps * tmp_N;
                double new_dot_product = dot(r.u, tmp_N);
                double new_dot_product_squared = std::pow(new_dot_product, 2);
                if (1 - std::pow(refraction_index, 2) * (1 - new_dot_product_squared) > 0.){
                    double k_0 = std::pow((refractive_index_object - 1.), 2) / std::pow((refractive_index_object + 1.), 2);
                    double R = k_0 + (1. - k_0) * std::pow((1 - std::abs(new_dot_product)), 5);
                    double random_value = uniform(engine);
                    if (random_value < R){
                        Ray reflected_ray = Ray(P, r.u - 2. * dot_product * N);
                        return getColor(reflected_ray, ray_depth - 1);
                    }
                    else{
                        Vector refracted_direction_tangential =  refraction_index * (r.u - new_dot_product * tmp_N);
                        Vector refracted_direction_normal = (-1.0) * tmp_N * std::sqrt(1. - std::pow(refraction_index, 2) * (1. - new_dot_product_squared));
                        Ray refracted_ray(P, refracted_direction_normal + refracted_direction_tangential);
                        return getColor(refracted_ray, ray_depth - 1);
                    }
                }
                else{
                    Ray reflected_ray = Ray(P, r.u - 2. * dot_product * N);
                    return getColor(reflected_ray, ray_depth - 1);
                }

            }
            else {
                P = P + eps * N;
                Vector albedo = intersection.albedo;
                Vector direction = light_source - P;
                double direction_norm = direction.norm();
                direction.normalize();
                    
                Ray shadow_ray(P, direction);
                Intersection intersection_shadow = this->intersect(shadow_ray);

                double visibility = 1.;

                if (intersection_shadow.hit){
                    Vector tmp_P = intersection_shadow.point;
                    if ((tmp_P - P).norm() < direction_norm) {
                        visibility = 0.;
                    }
                }
                Lo = 1. / (4. * M_PI * std::pow(direction_norm, 2)) * (1./ M_PI) * std::max(0.,dot(N, direction)) * albedo * visibility;
                Ray randomRay = Ray(P, random_cos(N));
                Lo = Lo + albedo * getColor(randomRay, ray_depth - 1);
                return Lo;
        }
        }

        else{
            return Vector(0.,0.,0.);
        }

    }

    Vector light_source;
    std::vector<Geometry*> objects;

};

void boxMuller(double stdev , double &x, double &y) { 
    double r1 = uniform(engine);
    double r2 = uniform(engine);
    x = std::sqrt((-2.) * std::log(r1)) * std::cos(2 * M_PI * r2) * stdev;
    y = std::sqrt((-2.) * std::log(r1)) * std::sin(2 * M_PI * r2) * stdev;
};


int main() {
    auto start = std::chrono::high_resolution_clock::now();

    int W = 512;
    int H = 512;
    double fov = 60 * M_PI / 180.;
    double z = - W / (2 * tan(fov / 2.));

    Vector albedo(1.,1.,1.);

    // Sphere S(Vector(0.,0.,0.), 10, albedo, 1.5);

    // Sphere mirror(Vector(-20,0,0), 10, albedo, 1.5, true);
    // Sphere refracted(Vector(0,0,0), 10, albedo, 1.5, false, true);
    // Sphere hollow_ext(Vector(20,0,0), 10, albedo, 1.5, false, true);
    // Sphere hollow_int(Vector(20,0,0), 9.8, albedo, 1.5, false, true, true);

    Vector Q(0,0,55);
    Vector light_source(-10,20,40);
    Scene scene(light_source);
    Sphere sphere1(Vector(0,1000,0), 940, Vector(1,0,0));
    Sphere sphere2(Vector(0,0,-1000), 940, Vector(0,1,0));
    Sphere sphere3(Vector(0,-1000,0), 990, Vector(0,0,1));
    Sphere sphere4(Vector(0,0,1000), 940, Vector(1,0,1));
    Sphere sphere5(Vector(1000, 0, 0), 940, Vector(1, 1, 0));
    Sphere sphere6(Vector(-1000, 0, 0), 940, Vector(0, 1, 1));
  
    // scene.addGeometry(&mirror);
    // scene.addGeometry(&refracted);
    // scene.addGeometry(&hollow_ext);
    // scene.addGeometry(&hollow_int);
    // scene.addGeometry(&S);

    scene.addGeometry(&sphere1);
    scene.addGeometry(&sphere2);
    scene.addGeometry(&sphere3);
    scene.addGeometry(&sphere4);
    scene.addGeometry(&sphere5);
    scene.addGeometry(&sphere6);

    TriangleMesh cat = TriangleMesh(Vector(0, -10, 0), 0.6, Vector(1., 1., 1.));
    cat.readOBJ("cat.obj");
    cat.BVH(cat.root, 0, cat.indices.size());
    scene.addGeometry(&cat);
    
    double light_intensity = 2e10;
    double gamma_correction = 2.2;
    double num_rays = 1024;

    std::vector<unsigned char> image(W * H * 3, 0);


    #pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {
            Vector color_pixel = Vector(0.,0.,0.);
            double x, y;

            for (int k = 0; k < num_rays; k++){          
                boxMuller(0.5, x, y);
                Vector random_dir = Vector(j + y - W / 2 + 0.5, H / 2 - i - x - 0.5, z);
                random_dir.normalize();
                Ray r(Q, random_dir);
                color_pixel = color_pixel + light_intensity * scene.getColor(r,5);
            }
            color_pixel = color_pixel / num_rays;
            image[(i * W + j) * 3 + 0] = std::min(std::pow(color_pixel[0], 1./gamma_correction), 255.);
            image[(i * W + j) * 3 + 1] = std::min(std::pow(color_pixel[1], 1./gamma_correction), 255.);
            image[(i * W + j) * 3 + 2] = std::min(std::pow(color_pixel[2], 1./gamma_correction), 255.);

        }
    }
    
    auto end = std::chrono::high_resolution_clock::now();
    stbi_write_png("image13.png", W, H, 3, &image[0], 0);

   
    std::chrono::duration<double> duration = end - start;
    double seconds = duration.count();
    std::cout << "Time taken: " << seconds << " seconds" << std::endl;

    return 0;
}