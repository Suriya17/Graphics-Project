#include "colorlib.h"
#define WALL_SIDE 499

int raster[WALL_SIDE+1][WALL_SIDE+1][3];
double hsvraster[WALL_SIDE+1][WALL_SIDE+1][3];
double zBuffer[WALL_SIDE+1][WALL_SIDE+1];
double Ka = 0.03, Ia = 0.5, Kd = 1, Ip = 1000, Ks = 1, specN = 4,c1 = 1, c2 = 1, c3 = 1;
int cubeR, cubeG, cubeB;
int sideR, sideG, sideB;
int topdowR, topdowG, topdowB;
int backR, backG, backB;
double prevIntensity = -1;
int r;


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

Vec3f light;

Vec3f rotateX(Vec3f p, double theta){
    Vec3f ans;
	ans.y = cos(theta)*p.y - sin(theta)*p.z;
	ans.z= cos(theta)*p.z + sin(theta)*p.y;
	ans.x = p.x;
	return ans;
}

Vec3f rotateY(Vec3f p, double theta){
	Vec3f ans;
	ans.x = cos(theta)*p.x + sin(theta)*p.z;
	ans.z = cos(theta)*p.z - sin(theta)*p.x;
	ans.y = p.y;
	return ans;
}

Vec3f rotateZ(Vec3f p, double theta){
	Vec3f ans;
	ans.x = cos(theta)*p.x - sin(theta)*p.y;
	ans.y = cos(theta)*p.y + sin(theta)*p.x;
	ans.z = p.z;
	return ans;
}


inline bool isObtuse(Vec3f a, Vec3f v){
    return ((a.x * v.x + a.y * v.y + a.z * v.z) < 0);
}

Vec3f cross_product(Vec3f v1, Vec3f v2){
    Vec3f ans;
    ans.x = v1.y*v2.z - v1.z*v2.y;
    ans.y = v1.z*v2.x - v1.x*v2.z;
    ans.z = v1.x*v2.y - v1.y*v2.x;
    
    ans.normalize();
    // ans.print();
    return ans;
}



class Plane
{
public:
    Vec3f a,b,c,d;
    Vec3f normal;
    int maxX,maxY,maxZ;
    int minX,minY,minZ;
    int R,G,B;
    Plane() {}
    Plane(Vec3f aa, Vec3f bb, Vec3f cc, Vec3f dd, int RR, int GG, int BB) : a(aa), b(bb), c(cc), d(dd), R(RR), G(GG) , B(BB){
        normal = cross_product(a-b,b-c);
        maxX = maximum(a.x,b.x,c.x,d.x);
        minX = minimum(a.x,b.x,c.x,d.x);
        maxY = maximum(a.y,b.y,c.y,d.y);
        minY = minimum(a.y,b.y,c.y,d.y);
        maxZ = maximum(a.z,b.z,c.z,d.z);
        minZ = minimum(a.z,b.z,c.z,d.z);
        // cout << maxX << " " << minX <<"  "<< maxY << " " << minY <<"  "<< maxZ << " " << minZ << endl; 
    }

    bool isOnFace(Vec3f pt){
        Vec3f P1_P4 = a - d;
        Vec3f P3_P4 = c - d;
        Vec3f TWO_P_C = pt*2 - a - c;
        return (P3_P4.dot(TWO_P_C - P3_P4) <= 1e-5 && P3_P4.dot(TWO_P_C + P3_P4) >= -1e-5) &&
         (P1_P4.dot(TWO_P_C - P1_P4) <= 1e-5 && P1_P4.dot(TWO_P_C + P1_P4) >= -1e-5);
        // return (P3_P4.dot(TWO_P_C - P3_P4) <= 0 && P3_P4.dot(TWO_P_C + P3_P4) >= 0) &&
        // (P1_P4.dot(TWO_P_C - P1_P4) <= 0 && P1_P4.dot(TWO_P_C + P1_P4) >= 0);
    }
};


Plane back_wall(Vec3f(0,0,-WALL_SIDE), Vec3f(WALL_SIDE,0,-WALL_SIDE), Vec3f(WALL_SIDE,WALL_SIDE,-WALL_SIDE), Vec3f(0,WALL_SIDE,-WALL_SIDE),128,0,255);
Plane up_wall(Vec3f(0,WALL_SIDE,0), Vec3f(0,WALL_SIDE,-WALL_SIDE), Vec3f(WALL_SIDE,WALL_SIDE,-WALL_SIDE), Vec3f(WALL_SIDE,WALL_SIDE,0),255,0,128);
Plane down_wall(Vec3f(0,0,0), Vec3f(WALL_SIDE,0,0), Vec3f(WALL_SIDE,0,-WALL_SIDE), Vec3f(0,0,-WALL_SIDE),255,128,255);
Plane left_wall(Vec3f(0,0,0), Vec3f(0,0,-WALL_SIDE), Vec3f(0,WALL_SIDE,-WALL_SIDE), Vec3f(0,WALL_SIDE,0),255,90,90);
Plane right_wall(Vec3f(WALL_SIDE,0,0), Vec3f(WALL_SIDE,WALL_SIDE,0), Vec3f(WALL_SIDE,WALL_SIDE,-WALL_SIDE), Vec3f(WALL_SIDE,0,-WALL_SIDE),255,0,255);

Plane walls[5] = {back_wall,up_wall,down_wall,left_wall,right_wall};

class Cube
{
    public:
        Vec3f center;
        double side;
        Vec3f vertices[12];
        Plane faces[6];

        Cube() {}
        Cube(Vec3f cc, double ss) : center(cc), side(ss){
            vertices[0] = Vec3f(center.x - side/2,center.y - side/2,center.z + side/2);
            vertices[1] = Vec3f(center.x + side/2,center.y - side/2,center.z + side/2);
            vertices[2] = Vec3f(center.x + side/2,center.y + side/2,center.z + side/2);
            vertices[3] = Vec3f(center.x - side/2,center.y + side/2,center.z + side/2);
            vertices[4] = Vec3f(center.x + side/2,center.y + side/2,center.z - side/2);
            vertices[5] = Vec3f(center.x + side/2,center.y - side/2,center.z - side/2);
            vertices[6] = Vec3f(center.x - side/2,center.y - side/2,center.z - side/2);
            vertices[7] = Vec3f(center.x - side/2,center.y + side/2,center.z - side/2);

            faces[0] = Plane(vertices[0], vertices[1], vertices[2], vertices[3],cubeR,cubeG,cubeB);
            faces[1] = Plane(vertices[1], vertices[5], vertices[4], vertices[2],cubeR,cubeG,cubeB);
            faces[2] = Plane(vertices[2], vertices[4], vertices[7], vertices[3],cubeR,cubeG,cubeB);
            faces[3] = Plane(vertices[0], vertices[3], vertices[7], vertices[6],cubeR,cubeG,cubeB);
            faces[4] = Plane(vertices[0], vertices[6], vertices[5], vertices[1],cubeR,cubeG,cubeB);
            faces[5] = Plane(vertices[4], vertices[5], vertices[6], vertices[7],cubeR,cubeG,cubeB);
        }

        Cube(Cube originalCube,double thetaX, double thetaY, double thetaZ){
            side = originalCube.side;
            center = originalCube.center;
            for (int i = 0; i < 8; ++i){
                this->vertices[i] = rotateX(rotateY(rotateZ(\
                originalCube.vertices[i] - originalCube.center,thetaZ),thetaY),thetaX) + originalCube.center;
                // this->vertices[i].print();
            }
            this->faces[0] = Plane(this->vertices[0], this->vertices[1], this->vertices[2], this->vertices[3],cubeR,cubeG,cubeB);
            this->faces[1] = Plane(this->vertices[1], this->vertices[5], this->vertices[4], this->vertices[2],cubeR,cubeG,cubeB);
            this->faces[2] = Plane(this->vertices[2], this->vertices[4], this->vertices[7], this->vertices[3],cubeR,cubeG,cubeB);
            this->faces[3] = Plane(this->vertices[0], this->vertices[3], this->vertices[7], this->vertices[6],cubeR,cubeG,cubeB);
            this->faces[4] = Plane(this->vertices[0], this->vertices[6], this->vertices[5], this->vertices[1],cubeR,cubeG,cubeB);
            this->faces[5] = Plane(this->vertices[4], this->vertices[5], this->vertices[6], this->vertices[7],cubeR,cubeG,cubeB);
        }
};


class face3{
public:
    int x,y,z;
    Vec3f normal;
    int R,G,B;
    face3(int aa, int bb, int cc, Vec3i col) : R(col.x), G(col.y) , B(col.z){
        x = aa;
        y = bb;
        z = cc;
    }

};

class mesh{
public:
    vector<Vec3f> vertices;
    vector<face3> faces;
    Vec3f centre;
    int num_vertices;
    int num_faces;
    double maxX = -1e5,maxY = -1e5,maxZ = -1e5;
    double minX = 1e5,minY = 1e5,minZ = 1e5;

    mesh(string objfile, Vec3i col,Vec3f pos,int scale_factor){
        ifstream objstream(objfile.c_str());
        string line;
        num_vertices = 0;
        num_faces = 0;
        while(getline(objstream,line)){
            if(line.c_str()[0] == '#')
                continue;

            if(line.c_str()[0] == 'v' && line.c_str()[1] == ' '){
                stringstream s(line);
                Vec3f vertex;
                char tmp;
                s >> tmp >> vertex.x >> vertex.y >> vertex.z;
                vertex.x *= scale_factor;
                vertex.x += pos.x;
                vertex.y *= scale_factor;
                vertex.y += pos.y;
                vertex.z *= scale_factor;
                vertex.z += pos.z;
                maxX = max(maxX,vertex.x);
                maxY = max(maxY,vertex.y);
                maxZ = max(maxZ,vertex.z);
                minX = min(minX,vertex.x);
                minY = min(minY,vertex.y);
                minZ = min(minZ,vertex.z);
                vertices.push_back(vertex);
                cout << vertex.x << " " << vertex.y << " " << vertex.z << endl;
            }

            centre.x = (maxX + minX)/2 + pos.x;
            centre.y = (maxY + minY)/2 + pos.y;
            centre.z = (maxZ + minZ)/2 + pos.x;

            if(line.c_str()[0] == 'f' && line.c_str()[1] == ' '){
                stringstream s(line);
                int a, b, c;
                char tmp;
                s >> tmp >> a >> b >> c;
                // Plane p(vertices[a-1],vertices[b-1],vertices[c-1],col);
                face3 p(a-1,b-1,c-1,col);
                
                p.normal = cross_product(vertices[a-1] - vertices[b-1],vertices[b-1] - vertices[c-1]);
                faces.push_back(p);

            }
        }
        cout << maxX << " " << maxY << " " << maxZ << endl;
        cout << minX << " " << minY << " " << minZ << endl;
        objstream.close();
        
        num_faces = faces.size();
        num_vertices = vertices.size();
        // cout << num_faces << "  " << num_vertices << endl;
        
    }

    mesh(mesh originalMesh,double thetaX, double thetaY, double thetaZ){
        centre = originalMesh.centre;
        for (int i = 0; i < originalMesh.num_vertices; ++i){
            this->vertices.push_back(rotateX(rotateY(rotateZ(\
            originalMesh.vertices[i] - originalMesh.centre,thetaZ),thetaY),thetaX) + originalMesh.centre);
            // this->vertices[i].print();
        }

        for (int i = 0; i < originalMesh.num_faces; ++i){
            faces.push_back(originalMesh.faces[i]);
            Vec3f a = vertices[faces[i].x];
            Vec3f b = vertices[faces[i].y];
            Vec3f c = vertices[faces[i].z];
            faces[i].normal = cross_product(a-b,b-c);
        }
    }
};

bool isOnface(mesh mymesh, face3 triface, Vec3f pt){
    Vec3f v0 = mymesh.vertices[triface.z] - mymesh.vertices[triface.x];
    Vec3f v1 = mymesh.vertices[triface.y] - mymesh.vertices[triface.x];
    Vec3f v2 = pt - mymesh.vertices[triface.x];

    double dot00 = v0.dot(v0);
    double dot01 = v0.dot(v1);
    double dot02 = v0.dot(v2);
    double dot11 = v1.dot(v1);
    double dot12 = v1.dot(v2);

    double invDenom = 1/(dot11*dot00 -dot01*dot01);

    double u = (dot11 * dot02 - dot01 * dot12) * invDenom;
    double v = (dot00 * dot12 - dot01 * dot02) * invDenom;

    return (u > -1e-5) and (v > -1e-5) and (u + v < 1+1e-5);
}

double illumPoint(Vec3f p, Vec3f cam, Vec3f L, mesh mymesh, face3 triface, int facid){
    int param = 1;
    double ambient = Ka * Ia;
    Vec3f n = triface.normal;
    Vec3f l = (L - p).normalize();
    Vec3f v = (cam - p).normalize();
    double dist = (cam - p).length();
    double fatt = min(1.0,1/(c1 + c2 * dist + c3 *dist*dist));
    double diffusion = Kd * Ip * fatt * l.dot(n);
    double specular = Ks * Ip * fatt * pow(((2 * l.dot(n) * v.dot(n)) - l.dot(v)),specN);

    for(int i = 0; i < mymesh.num_faces; i++){
        if(i != facid){
            double t1 = (L-p).dot(mymesh.faces[i].normal);
            double t2 = (mymesh.vertices[mymesh.faces[i].x] - p).dot(mymesh.faces[i].normal);

            if(t1 == 0)
                continue;

            double t = t2/t1;
            Vec3f pt = p + (L-p)*t;
            // cout << t << endl;
            if(t > 0 && t < 1 && isOnface(mymesh,mymesh.faces[i],pt)){
                param = 0;
                break;
            }
        }
    }
    double ans =  ambient + param * max(diffusion + specular,0.0);
    if( ans < 0)
        cout << "YO"<< endl;
    return ans;
}


double illumPoint(Vec3f p, Vec3f cam, Vec3f L, mesh mymesh, Plane wall, int facid){
    int param = 1;
    double ambient = Ka * Ia;
    Vec3f n = wall.normal;
    Vec3f l = (L - p).normalize();
    Vec3f v = (cam - p).normalize();
    double dist = (cam - p).length();
    double fatt = min(1.0,1/(c1 + c2 * dist + c3 *dist*dist));
    double diffusion = Kd * Ip * fatt * l.dot(n);
    double specular = Ks * Ip * fatt * pow(((2 * l.dot(n) * v.dot(n)) - l.dot(v)),specN);

    for(int i = 0; i < mymesh.num_faces; i++){
        if(i != facid){
            double t1 = (L-p).dot(mymesh.faces[i].normal);
            double t2 = (mymesh.vertices[mymesh.faces[i].x] - p).dot(mymesh.faces[i].normal);

            if(t1 == 0)
                continue;

            double t = t2/t1;
            Vec3f pt = p + (L-p)*t;
            // cout << t << endl;

            // if(isOnface(mymesh,mymesh.faces[i],pt))
            //     cout << isOnface(mymesh,mymesh.faces[i],pt) << endl;
            if(t > 0 && t < 1 && isOnface(mymesh,mymesh.faces[i],pt)){
                // cout << "This ah" << endl;
                param = 0;
                break;
            }
        }
    }
    double ans =  ambient + param * max(diffusion + specular,0.0);
    if( ans < 0)
        cout << "YO"<< endl;
    return ans;

}


double illumPoint(Vec3f p, Vec3f cam, Vec3f L, Cube cube, Plane myface, int facid){
    int param = 1;
    double ambient = Ka * Ia;
    Vec3f n = myface.normal;
    Vec3f l = (L - p).normalize();
    Vec3f v = (cam - p).normalize();
    double dist = (cam - p).length();
    double fatt = min(1.0,1/(c1 + c2 * dist + c3 *dist*dist));
    double diffusion = Kd * Ip * fatt * l.dot(n);
    double specular = Ks * Ip * fatt * pow(((2 * l.dot(n) * v.dot(n)) - l.dot(v)),specN);
    
    // cout << "ambient - " << ambient << " diff - " << diffusion << " spec -" << specular << endl;
    for(int i = 0; i < 6; i++){
        if(i != facid){
            double t1 = (L-p).dot(cube.faces[i].normal);
            double t2 = (cube.faces[i].a - p).dot(cube.faces[i].normal);

            if(t1 == 0)
                continue;

            double t = t2/t1;
            Vec3f pt = p + (L-p)*t;
            // cout << t << endl;
            if(t > 0 && t < 1 && cube.faces[i].isOnFace(pt)){
                param = 0;
                break;
            }
        }
    }
    double ans =  ambient + param * max(diffusion + specular,0.0);
    if( ans < 0)
        cout << "YO"<< endl;
    return ans;
}

void make_ppm(string filename){
    ofstream file;
	file.open(filename.c_str());
	file << "P3\n";
	file << WALL_SIDE+1 << " " << WALL_SIDE+1 << "\n";
	file << "255\n";
	// cout << "Hello" << endl;
	
	for(int i = WALL_SIDE; i >= 0; i--){

		for(int j = 0; j <= WALL_SIDE; j++){
			file << raster[j][i][0] << " " << raster[j][i][1] << " " << raster[j][i][2] << " ";
		}
		
		file << "\n";
	}
	
	file << "\n";
	file.close();
}

bool find(Vec3f cam, Vec3f pix, Plane myface, Cube cube,int facid){
    double t1 = (pix - cam).dot(myface.normal);
    int i = pix.x;
    int j = pix.y;
    
    if(t1 == 0)
        return false;
    Vec3f arbpt = (myface.a + myface.b)/2.0;
    double t2 = (arbpt - cam).dot(myface.normal);
    double t = t2/t1;
    if(t < 1)
        return false;

    Vec3f pt = cam + (pix - cam)*t;

    if(pt.z <= zBuffer[i][j] && myface.isOnFace(pt)){
        double fR = myface.R/255.0, fG = myface.G/255.0, fB = myface.B/255.0;
        RGBtoHSV(fR,fG,fB,hsvraster[i][j][0],hsvraster[i][j][1],hsvraster[i][j][2]);
        Vec3f tempvec;
        if(r == 0)
            hsvraster[i][j][2] = illumPoint(pt, cam, light, cube, myface,facid);
        else{
            hsvraster[i][j][2] = 0;
            
            for(int d = 0; d <= r; d++){
                for(int e = 0; e <=r; e++){
                    for(int f = 0; f <=r; f++)
                    tempvec = light - Vec3f(d,e,f);
                    hsvraster[i][j][2] += illumPoint(pt, cam, tempvec, cube, myface,facid);
                }
            }
        }
        // if(i == 0 and j == 0)
        //     cout << hsvraster[i][j][2] << endl;
        // cout << hsvraster[i][j][2] << endl;
        // RGBtoHSV(raster[i][j][0],raster[i][j][1],raster[i][j][2],hsvraster[i][j][0],hsvraster[i][j][1],hsvraster[i][j][2]);
        return true;
    }
    return false;
}


bool find(){

}


bool find(Vec3f cam, Vec3f pix, face3 triface, mesh mymesh ,int facid){
    double t1 = (pix - cam).dot(triface.normal);
    int i = pix.x;
    int j = pix.y;
    
    if(t1 == 0)
        return false;
    Vec3f arbpt = (mymesh.vertices[triface.x] + mymesh.vertices[triface.y] + mymesh.vertices[triface.z])/3.0;
    double t2 = (arbpt - cam).dot(triface.normal);
    double t = t2/t1;
    if(t < 1)
        return false;

    Vec3f pt = cam + (pix - cam)*t;

    if(pt.z <= zBuffer[i][j] && isOnface(mymesh,triface,pt)){
        double fR = triface.R/255.0, fG = triface.G/255.0, fB = triface.B/255.0;
        RGBtoHSV(fR,fG,fB,hsvraster[i][j][0],hsvraster[i][j][1],hsvraster[i][j][2]);
        Vec3f tempvec;
        if(r == 0)
            hsvraster[i][j][2] = illumPoint(pt, cam, light, mymesh, triface,facid);
        else{
            hsvraster[i][j][2] = 0;
            
            for(int d = 0; d <= r; d++){
                for(int e = 0; e <=r; e++){
                    for(int f = 0; f <=r; f++)
                    tempvec = light - Vec3f(d,e,f);
                    hsvraster[i][j][2] += illumPoint(pt, cam, tempvec, mymesh, triface,facid);
                }
            }
        }
        
        return true;
    }
    return false;

}

bool find(Vec3f cam, Vec3f pix, mesh mymesh, Plane wall ,int facid){
    double t1 = (pix - cam).dot(wall.normal);
    int i = pix.x;
    int j = pix.y;
    
    if(t1 == 0)
        return false;
    Vec3f arbpt = (wall.a + wall.b + wall.c)/3.0;
    double t2 = (arbpt - cam).dot(wall.normal);
    double t = t2/t1;
    if(t < 1)
        return false;

    Vec3f pt = cam + (pix - cam)*t;

    if(pt.z <= zBuffer[i][j] && wall.isOnFace(pt)){
        double fR = wall.R/255.0, fG = wall.G/255.0, fB = wall.B/255.0;
        RGBtoHSV(fR,fG,fB,hsvraster[i][j][0],hsvraster[i][j][1],hsvraster[i][j][2]);
        Vec3f tempvec;
        if(r == 0)
            hsvraster[i][j][2] = illumPoint(pt, cam, light, mymesh, wall,facid);
        else{
            hsvraster[i][j][2] = 0;
            
            for(int d = 0; d <= r; d++){
                for(int e = 0; e <=r; e++){
                    for(int f = 0; f <=r; f++)
                    tempvec = light - Vec3f(d,e,f);
                    hsvraster[i][j][2] += illumPoint(pt, cam, tempvec, mymesh, wall,facid);
                }
            }
        }
        
        return true;
    }
    return false;

}

map<pair<int,int>,double> hsImap;




void raycast(Cube cube, Vec3f cam, string filename){
    bool visible_faces[6];
    std::memset(raster, 0, sizeof(raster));
    // cout << raster[20][20][2] << endl;
    for(int i = 0; i <= WALL_SIDE; i++){
        for(int j = 0; j <= WALL_SIDE; j++){
            zBuffer[i][j] = 1e8;
        }
    }
    for(int i = 0; i < 6; i++){
        Vec3f shootVec = cube.faces[i].a - cam;
        
        if(isObtuse(shootVec,cube.faces[i].normal)){
            visible_faces[i] = true;
        }
        else
            visible_faces[i] = false;
    }
    
    for(int i = 0; i <= WALL_SIDE; i++){
        for(int j = 0; j <= WALL_SIDE; j++){
            Vec3f pix(i,j,0);
            bool isColored = false;
            for(int k = 0; k < 6;k++){
                if(visible_faces[k]){
                    bool tmp = find(cam,pix,cube.faces[k],cube,k);
                    isColored = tmp or isColored;
                }
            }
            if(isColored)
                continue;
            else{
                for(int k = 0; k < 5; k++){
                    find(cam,pix,walls[k],cube,-1);
                }
            }
        }
    }
    // double Imax = -1e8;
    
    pair<int,int> key;
    map<pair<int,int>,double>::const_iterator it;

    for(int i = 0; i <= WALL_SIDE; i++){
        for(int j = 0; j <= WALL_SIDE; j++){
            key.first = (int)(hsvraster[i][j][0] * 1000) ;
            key.second = (int)(hsvraster[i][j][1] * 100000);
            it = hsImap.find(key);
            if( it == hsImap.end()){
                hsImap[key] = hsvraster[i][j][2];
            }
            else if(it->second < hsvraster[i][j][2]){
                hsImap[key] = hsvraster[i][j][2];
            }
            // Imax = max(hsvraster[i][j][2],Imax);
            
        }
    }
    double fR,fG,fB;
    for(int i = 0; i <= WALL_SIDE; i++){
        for(int j = 0; j <= WALL_SIDE; j++){
            key.first = (int)(hsvraster[i][j][0] * 1000) ;
            key.second = (int)(hsvraster[i][j][1] * 100000);
            it = hsImap.find(key);
            // if((i == 50 && j == 50) or (i== 200 && j == 200)){
            //     cout << it->second << endl;
            // }
            hsvraster[i][j][2] = hsvraster[i][j][2]/it->second;

            HSVtoRGB(fR,fG,fB,hsvraster[i][j][0],hsvraster[i][j][1],hsvraster[i][j][2]);
            raster[i][j][0] = ((int)(fR * 255));
            raster[i][j][1] = ((int)(fG * 255));
            raster[i][j][2] = ((int)(fB * 255));
        }
    }

    make_ppm(filename);
}

void raycast(mesh mymesh, Vec3f cam, string filename){
    bool visible_faces[mymesh.num_faces];
    std::memset(raster, 0, sizeof(raster));

    for(int i = 0; i <= WALL_SIDE; i++){
        for(int j = 0; j <= WALL_SIDE; j++){
            zBuffer[i][j] = 1e8;
        }
    }

    for(int i = 0; i < mymesh.num_faces; i++){
        Vec3f shootVec = mymesh.vertices[mymesh.faces[i].x] - cam;
        
        if(isObtuse(shootVec,mymesh.faces[i].normal)){
            visible_faces[i] = true;
        }
        else
            visible_faces[i] = false;
    }


    for(int i = 0; i <= WALL_SIDE; i++){
        for(int j = 0; j <= WALL_SIDE; j++){
            Vec3f pix(i,j,0);
            bool isColored = false;
            for(int k = 0; k < mymesh.num_faces;k++){
                if(visible_faces[k]){
                    bool tmp = find(cam,pix,mymesh.faces[k],mymesh,k);
                    isColored = tmp or isColored;
                }
            }
            if(isColored)
                continue;
            else{
                for(int k = 0; k < 5; k++){
                    find(cam,pix,mymesh,walls[k],-1);
                }
            }
        }
    }

    pair<int,int> key;
    map<pair<int,int>,double>::const_iterator it;

    for(int i = 0; i <= WALL_SIDE; i++){
        for(int j = 0; j <= WALL_SIDE; j++){
            key.first = (int)(hsvraster[i][j][0] * 1000) ;
            key.second = (int)(hsvraster[i][j][1] * 100000);
            it = hsImap.find(key);
            if( it == hsImap.end()){
                hsImap[key] = hsvraster[i][j][2];
            }
            else if(it->second < hsvraster[i][j][2]){
                hsImap[key] = hsvraster[i][j][2];
            }
            // Imax = max(hsvraster[i][j][2],Imax);
            
        }
    }
    double fR,fG,fB;
    for(int i = 0; i <= WALL_SIDE; i++){
        for(int j = 0; j <= WALL_SIDE; j++){
            key.first = (int)(hsvraster[i][j][0] * 1000) ;
            key.second = (int)(hsvraster[i][j][1] * 100000);
            it = hsImap.find(key);
            // if((i == 50 && j == 50) or (i== 200 && j == 200)){
            //     cout << it->second << endl;
            // }
            hsvraster[i][j][2] = hsvraster[i][j][2]/it->second;

            HSVtoRGB(fR,fG,fB,hsvraster[i][j][0],hsvraster[i][j][1],hsvraster[i][j][2]);
            raster[i][j][0] = ((int)(fR * 255));
            raster[i][j][1] = ((int)(fG * 255));
            raster[i][j][2] = ((int)(fB * 255));
        }
    }

    make_ppm(filename);

}



// int main(int argc, char** argv){

//     int whiteLight;

//     cout << "1 for white light and 0 for colored light and white objs" << endl;
//     cin >> whiteLight;

//     if(whiteLight){
//         cout << "RGB of cube" << endl;
//         cin >> cubeR >> cubeG >> cubeB;
//         cout << "RGB of side walls" << endl;
//         cin >> sideR >> sideG >> sideB;
//         cout << "RGB of top and bottom walls" << endl;
//         cin >> topdowR >> topdowG >> topdowB;
//         cout << "RGB of back wall" << endl;
//         cin >> backR >> backG >> backB;

//         for(int i = 0; i < 6; i++){
//             switch(i){
//                 case 0:
//                     walls[i].R = backR;
//                     walls[i].G = backG;
//                     walls[i].B = backB;
//                 break;
//                 case 1:
//                 case 2:
//                     walls[i].R = topdowR;
//                     walls[i].G = topdowG;
//                     walls[i].B = topdowB;
//                 break;
//                 case 3:
//                 case 4:
//                     walls[i].R = sideR;
//                     walls[i].G = sideG;
//                     walls[i].B = sideB;
//                 break;
//             }
//         }
//     }
//     else{
//         cin >> cubeR >> cubeG >> cubeB;
//         sideR = cubeR;
//         topdowR = cubeR;
//         backR = cubeR;

//         sideG = cubeG;
//         topdowG = cubeG;
//         backG = cubeG;

//         sideB = cubeB;
//         topdowB = cubeB;
//         backB = cubeB;

//         for(int i = 0; i < 5; i++){
//             switch(i){
//                 case 0:
//                     walls[i].R = backR;
//                     walls[i].G = backG;
//                     walls[i].B = backB;
//                 break;
//                 case 1:
//                 case 2:
//                     walls[i].R = topdowR;
//                     walls[i].G = topdowG;
//                     walls[i].B = topdowB;
//                 break;
//                 case 3:
//                 case 4:
//                     walls[i].R = sideR;
//                     walls[i].G = sideG;
//                     walls[i].B = sideB;
//                 break;
//             }
//         }
//         cout << walls[0].R << " " << walls[0].G << " "<< walls[0].B << endl;
//     }
//     int a,b,c,s;

//     cout << "Centre of cube :" << endl;
//     cin >> a >> b >> c;
//     Vec3f centre(a,b,c);
//     cout << "Location of camera :" << endl;
//     cin >> a >> b >> c;
//     Vec3f cam(a,b,c);
//     cout << "length of the side of cube :" << endl;
//     cin >> s;
//     Cube cube(centre,s);
//     cout << "Location of light :" << endl;
//     cin >> a >> b >> c;
//     light = Vec3f(a,b,c);
//     cout << "Radius of light source(0 for point source. for big r, too slow):" << endl;
//     cin >> r;

    
//     double x1,x2,x3;
// 	double thetaX = 0,thetaY = 0,thetaZ =0;
// 	double magnitude;
//     double dtheta = 0.05;
//     double theta = 0;
// 	for (int i = 0; i < 100; ++i)
// 	{
//         char name[20];
// 		sprintf(name,"cube%03d.ppm",i+1);
// 		string cppname(name);
// 		cout << "Plotting " << cppname << endl;
// 		Cube temp = Cube(cube,thetaX,thetaY,thetaZ);
// 		raycast(temp,cam,cppname);
// 		theta += dtheta;
// 		x1 =  cos(theta);
// 		x2 =  sin(2*theta);
// 		x3 =  sin(3*theta);
// 		magnitude = sqrt(x1*x1 + x2*x2 + x3*x3);
// 		thetaX = acos(x1/magnitude);
// 		thetaY = acos(x2/magnitude);
// 		thetaZ = acos(x3/magnitude);
// 	}
// }



int main(){
    string meshname("apple2.obj");
    Vec3i col(255,255,0);
    Vec3f pos(300,250,-150);
    mesh mymesh(meshname,col,pos,1000);
    Cube mycube(pos,200);
    Vec3f cam(250,250,250);
    light.x = 250;
    light.y = 250;
    light.z = -50;
    string output("out.ppm");

    // raycast(mymesh,cam,output);

    double x1,x2,x3;
    double thetaX = 0,thetaY = 0,thetaZ =0;
    double magnitude;
    double dtheta = 0.05;
    double theta = 0;

    for (int i = 0; i < 100; ++i)
	{
        char name[20];
		sprintf(name,"cube%03d.ppm",i+1);
		string cppname(name);
		cout << "Plotting " << cppname << endl;
        
		mesh temp = mesh(mymesh,thetaX,thetaY,thetaZ);
        if(i > 40 and i <= 60)
		    raycast(temp,cam,cppname);
		theta += dtheta;
		x1 =  cos(theta);
		x2 =  sin(2*theta);
		x3 =  sin(3*theta);
		magnitude = sqrt(x1*x1 + x2*x2 + x3*x3);
		thetaX = acos(x1/magnitude);
		thetaY = acos(x2/magnitude);
		thetaZ = acos(x3/magnitude);
	}
    // raycast(mycube,cam,output);
}