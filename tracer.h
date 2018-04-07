#include "colorlib.h"
#define WALL_SIDE 499

int raster[WALL_SIDE+1][WALL_SIDE+1][3];
double hsvraster[WALL_SIDE+1][WALL_SIDE+1][3];
double zBuffer[WALL_SIDE+1][WALL_SIDE+1];
double Ka = 0.03, Ia = 0.5, Kd = 1, Ip = 1000, Ks = 1, specN = 40,c1 = 1, c2 = 1, c3 = 1;
map<pair<int,int>,double> hsImap;
bool CrossShadows = false;

struct color{
    int R, G, B;
};

inline int maximum(int a, int b, int c, int d){
    return max(max(max(a,b),c),d);
}

inline int minimum(int a, int b, int c, int d){
    return min(min(min(a,b),c),d);
}

template<typename T>
class Vec3
{
public:
    T x, y, z;
    Vec3() : x(T(0)), y(T(0)), z(T(0)) {}
    Vec3(T xx) : x(xx), y(xx), z(xx) {}
    Vec3(T xx, T yy, T zz) : x(xx), y(yy), z(zz) {}
    Vec3& normalize()
    {
        T nor2 = length2();
        if (nor2 > 0) {
            T invNor = 1 / sqrt(nor2);
            x *= invNor, y *= invNor, z *= invNor;
        }
        return *this;
    }
    Vec3<T> operator * (const T &f) const { return Vec3<T>(x * f, y * f, z * f); }
    Vec3<T> operator * (const Vec3<T> &v) const { return Vec3<T>(x * v.x, y * v.y, z * v.z); }
    Vec3<T> operator / (const T &f) const { return Vec3<T>(x /f, y /f, z/f); }
    T dot(const Vec3<T> &v) const { return x * v.x + y * v.y + z * v.z; }
    Vec3<T> operator - (const Vec3<T> &v) const { return Vec3<T>(x - v.x, y - v.y, z - v.z); }
    Vec3<T> operator + (const Vec3<T> &v) const { return Vec3<T>(x + v.x, y + v.y, z + v.z); }
    Vec3<T>& operator += (const Vec3<T> &v) { x += v.x, y += v.y, z += v.z; return *this; }
    Vec3<T>& operator *= (const Vec3<T> &v) { x *= v.x, y *= v.y, z *= v.z; return *this; }
    Vec3<T> operator - () const { return Vec3<T>(-x, -y, -z); }
    T length2() const { return x * x + y * y + z * z; }
    T length() const { return sqrt(length2()); }
    Vec3<double> itof(){
        Vec3<double> ans(this->x,this->y,this->z);
        return ans;
    }
    void print(){
        cout << x << "  "<< y << "  "<< z << endl;
    }
};

typedef Vec3<double> Vec3f;
typedef Vec3<int> Vec3i;

Vec3f rotateX(Vec3f p, double theta);
Vec3f rotateY(Vec3f p, double theta);
Vec3f rotateZ(Vec3f p, double theta);
Vec3f unit_cross_product(Vec3f v1, Vec3f v2);

class Wall{
    public:
        Vec3f a, b, c, d;
        Vec3f n;
        color col;
    Wall(){}
    Wall(Vec3f aa, Vec3f bb, Vec3f cc, Vec3f dd, int RR, int GG, int BB) : a(aa), b(bb), c(cc), d(dd){
        n = unit_cross_product(a-b,b-c);
        col.R = RR;
        col.G = GG;
        col.B = BB;
    }
};

Wall back_wall(Vec3f(0,0,-WALL_SIDE), Vec3f(WALL_SIDE,0,-WALL_SIDE), Vec3f(WALL_SIDE,WALL_SIDE,-WALL_SIDE), Vec3f(0,WALL_SIDE,-WALL_SIDE),128,128,128);
Wall up_wall(Vec3f(0,WALL_SIDE,0), Vec3f(0,WALL_SIDE,-WALL_SIDE), Vec3f(WALL_SIDE,WALL_SIDE,-WALL_SIDE), Vec3f(WALL_SIDE,WALL_SIDE,0),0,0,200);
Wall down_wall(Vec3f(0,0,0), Vec3f(WALL_SIDE,0,0), Vec3f(WALL_SIDE,0,-WALL_SIDE), Vec3f(0,0,-WALL_SIDE),0,0,200);
Wall left_wall(Vec3f(0,0,0), Vec3f(0,0,-WALL_SIDE), Vec3f(0,WALL_SIDE,-WALL_SIDE), Vec3f(0,WALL_SIDE,0),128,128,128);
Wall right_wall(Vec3f(WALL_SIDE,0,0), Vec3f(WALL_SIDE,WALL_SIDE,0), Vec3f(WALL_SIDE,WALL_SIDE,-WALL_SIDE), Vec3f(WALL_SIDE,0,-WALL_SIDE),128,128,128);

Wall walls[5] = {back_wall,up_wall,down_wall,left_wall,right_wall};

class meshPlane
{
    public:
        int a, b, c;
        Vec3f n;
        color col;

        meshPlane(){}
        meshPlane(int aa, int bb, int cc, int RR, int GG, int BB) : a(aa), b(bb), c(cc){
            col.R = RR;
            col.G = GG;
            col.B = BB;
        }    
};

class mesh{
    public:
        vector<Vec3f> vertices;
        vector<meshPlane> faces;
        Vec3f centre;
        int num_vertices,num_faces;
        double maxX = -1e5,maxY = -1e5,maxZ = -1e5;
        double minX = 1e5,minY = 1e5,minZ = 1e5;

        mesh(string objfile, color col, Vec3f cen, double scale_factor);
        mesh(mesh &originalMesh,double &thetaX, double &thetaY, double &thetaZ);
};

bool isOnFace(Vec3f &pt, meshPlane &face,mesh &parent_mesh);
bool isOnFace(Vec3f &pt, Wall &w);
bool getIntersection(Vec3f p1, Vec3f p2, Vec3f &pt ,meshPlane &face, mesh &parent_mesh);
bool getIntersection(Vec3f p1, Vec3f p2, Vec3f &pt, Wall &w);
void make_ppm(string &filename);
void NormaliseV();
bool ShadowCheckIntersection(Vec3f &light, Vec3f pt, int faceInd,mesh &parent_mesh);
bool ShadowCheckIntersection(Vec3f &light, Vec3f pt, mesh &parent_mesh);
double illuminatePoint(Vec3f &pt, Vec3f &cam, Vec3f &light,int faceInd, mesh &parent_mesh);
double illuminatePoint(Vec3f &pt, Vec3f &cam, Vec3f &light,Wall &w, mesh &parent_mesh, int wallInd);
double illuminatePoint(Vec3f &pt, Vec3f &cam, Vec3f &light,int faceInd, int meshIndex, vector<mesh> &meshes);
double illuminatePoint(Vec3f &pt, Vec3f &cam, Vec3f &light,Wall &w, vector<mesh> &meshes, int wallInd);
void raycast(mesh &mymesh,Vec3f cam, Vec3f light, string &filename);
void raycast(vector<mesh> &meshes,Vec3f cam, Vec3f light, string &filename);



