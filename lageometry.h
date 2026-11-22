#pragma once
//#define DEBUG

#ifdef DEBUG
    #include <stdio.h>
#endif
#include <stdio.h>

#include <math.h>

namespace lag{

    typedef float mat4[4][4];
    typedef float mat3[3][3];
    struct vec3f;
    struct vec2f;

    struct vec2f{
        float x;
        float y;

        inline vec2f operator+(const float imfuckingdrunk){
            return {this->x + imfuckingdrunk, this->y + imfuckingdrunk};
        }

        inline vec2f operator+(const vec2f&other) const {
            return {this->x+other.x, this->y+other.y};
        }

        inline vec2f operator*(const float shitass) const {
            return {this->x*shitass, this->y*shitass};
        }  

        inline vec2f& operator+=(const vec2f&other){
            this->x += other.x;
            this->y += other.y;
            return *this;
        }

        inline void set(float shit){
            this->x = shit;
            this->y = shit;
        }

        static inline vec2f findLinCoeff(const vec2f & p1, const vec2f & p2){
            lag::vec2f coeff = {};
                coeff.x = (p1.y - p2.y) / (p1.x - p2.x);
                coeff.y = p1.y - coeff.x * p1.x; 
            return coeff;
        }

        inline float length() const{
            return sqrt(this->x * this->x + this->y * this->y);
        }

        inline float squareLength() const{
            return this->x * this->x + this->y * this->y;
        }
    };

    struct vec3f{
        float x;
        float y;
        float z;

        inline vec2f toVec2() const{
            return {this->x, this->y};
        }

        inline void rotareAlongX(float angle){
            float ca = cos(angle); float sa = sin(angle);
            float z = this->z;
            float y = this->y;

            this->y = y * ca + z * sa ;
            this->z = z * ca - y * sa ;
        }

        inline void rotareAlongY(float angle){
            float ca = cos(angle); float sa = sin(angle);
            float z = this->z;
            float x = this->x;

            this->x = x * ca - z * sa;
            this->z = z * ca + x * sa;
        }

        inline void rotareAlongZ(float angle){
            float ca = cos(angle); float sa = sin(angle);
            float x = this->x;
            float y = this->y;

            this->y = y * ca + x * sa;
            this->x = x * ca - y * sa;
        }

        inline float length() const{
            return sqrt(this->x * this->x + this->y * this->y + this->z * this->z);
        }

        inline float squareLength() const{
            return this->x * this->x + this->y * this->y + this->z * this->z;
        }

        inline vec3f operator+(const vec3f&other) const {
            return {this->x+other.x, this->y+other.y, this->z+other.z};
        }

        inline vec3f operator+(const vec2f&other) const {
            return {this->x+other.x, this->y+other.y, this->z};
        }

        inline vec3f operator*(const vec3f&other) const {
            return {this->x*other.x, this->y*other.y, this->z*other.z};
        }

        inline vec3f operator*(float sclr) const {
            return {this->x*sclr, this->y*sclr, this->z*sclr};
        }

        inline void set(float v){
            this->x = v; this->y = v; this->z = v; 
        }

        inline void set(const vec2f & other){
            this->x = other.x; this->y = other.y; 
        }
//        inline vec3f operator=(float shit) const {
//            return {shit, shit, shit};
//        }
//
//        inline vec3f operator=(const vec3f & other) const {
//            return {other.x, other.y, other.z};
//        }

        inline vec3f operator-(const vec3f&other) const {
            return {this->x-other.x, this->y-other.y, this->z-other.z};
        }

        inline float getYAngle() const{
            return atan2(this->z, this->x);
        }

        inline void normalize(){
            float length = this->length();
            if(length > 0.0f){
                this->x /= length;
                this->y /= length;
                this->z /= length;    
            }else{
                *this={0.0f, 0.0f, 0.0f};
            }
        } 

        inline void removeNaN(){
            if(isnan(this->x)) this->x = 0.0f;
            if(isnan(this->y)) this->y = 0.0f;
            if(isnan(this->z)) this->z = 0.0f;
        }

        inline void makeAbsolute(){
            this->x *= -1.0f + 2.0f *(this->x >= 0.0f);
            this->y *= -1.0f + 2.0f *(this->y >= 0.0f);
            this->z *= -1.0f + 2.0f *(this->z >= 0.0f);
        }

        inline vec3f& operator+=(const vec3f&other){
            this->x += other.x;
            this->y += other.y;
            this->z += other.z;
            return *this;
        }

        inline vec3f& operator*=(const vec3f&other){
            this->x *= other.x;
            this->y *= other.y;
            this->z *= other.z;
            return *this;
        }

        inline vec3f& operator+=(float sclr){
            this->x += sclr;
            this->y += sclr;
            this->z += sclr;
            return *this;
        }

        inline vec3f& operator+=(const vec2f&other){
            this->x += other.x;
            this->y += other.y;
            return *this;
        }

        inline vec3f& operator-=(const vec3f&other){
            this->x -= other.x;
            this->y -= other.y;
            this->z -= other.z;
            return *this;
        }

        inline vec3f& operator*=(const float (&mul)){
            this->x *= mul;
            this->y *= mul;
            this->z *= mul;
            return *this;
        }
        
        inline vec3f& operator*=(const float (&trns)[3][3]){
            vec3f tmp = {};
            tmp.x = this->x * trns[0][0] + this->y * trns[0][1] + this->z * trns[0][2];
            tmp.y = this->x * trns[1][0] + this->y * trns[1][1] + this->z * trns[1][2];
            tmp.z = this->x * trns[2][0] + this->y * trns[2][1] + this->z * trns[2][2];
            *this = tmp;
            return *this;      
        }

        inline vec3f& operator*=(const float (&trns)[4][4]){  // igneres last row and collumn
            vec3f tmp = {};
            tmp.x = this->x * trns[0][0] + this->y * trns[0][1] + this->z * trns[0][2];
            tmp.y = this->x * trns[1][0] + this->y * trns[1][1] + this->z * trns[1][2];
            tmp.z = this->x * trns[2][0] + this->y * trns[2][1] + this->z * trns[2][2];
            *this = tmp;
            return *this;      
        }
    };

    struct vec4f{
        float x;
        float y;
        float z;
        float w;
    };

    struct vect2df{
        vec2f position;
        vec3f color;
    };

    struct vector3df{
        vec3f position;
        vec3f color;
        vec2f texCoord;
    };

    struct vector3dPer{
        vec3f base;
        vec3f end;

        inline float length() const{return (end - base).length();}
        
        inline vec3f getNormalized() const {
            vec3f v = (end - base);
                v.normalize();
            return v;
        }
        
        inline void loadRand() {
            base.x = drand48();end.x = drand48();
            base.y = drand48();end.y = drand48();
            base.z = drand48();end.z = drand48();
        }

        inline void randomLength(bool norm = false){
            end.x = 2.0f * drand48() - 1.0f;
            end.y = 0.0f;
            end.z = 2.0f * drand48() - 1.0f;//drand48();

            if(norm)end.normalize();

            end += base;
        }

        inline float getDistance(const vec3f & point) const {
            return (point - this->base).length();
        }
    };

    struct PerlinNoice
    {
        vector3dPer * points;
        unsigned int width;
        unsigned int height;

        // xy_res number of points in axiesess
        // normalized returns normalized values 
        // reand positions - positnions are random not in the gread        
        PerlinNoice(unsigned int x_res, unsigned int y_res, bool normalized = true, bool randPositions = false);
        ~PerlinNoice();

        vector3dPer & getPoint(unsigned int x, unsigned int y);

        float getValuef(float x, float y);
    };

    struct Dmatrixf
    {
        size_t w = 0;
        size_t h = 0;

        float * mxData = nullptr;

        ~Dmatrixf(){delete[]  mxData;}

        void newMatrix(size_t width, size_t height);
        void reSize(size_t width, size_t height);

        float getValue(unsigned int x, unsigned int y)const;
        void setValue(unsigned int x, unsigned int y, float v);

        inline float fgetValue(unsigned int x, unsigned int y)const {    return *(this->mxData + (this->w * y + x));}
        inline void fsetValue(unsigned int x, unsigned int y, float v){    *(this->mxData + (this->w * y + x)) = v;}
    
        void loadTestValues();
    };
    
    

    //typedef mat3f float[3][3];
    // A*B=C
    void fastMultiply(float (&mA)[4][4], float (&matB)[4][4], float (&mC)[4][4]);
    void fastMultiply(mat3 & mA, mat3 & mB, mat3 & mC);

    void multiplyMx(const Dmatrixf & a, const Dmatrixf & b, Dmatrixf & c);

    void prontMatrix(const float (&matrix)[4][4]);
    void prontMatrix(const mat3 & matrix);
    void prontMatrix(const Dmatrixf & matrix);

    void makeIdentityMatrix(float(&matrix)[4][4]);
    void makeIdentityMatrix(mat3 & matrix);

    void project(float(&matrix)[4][4], float fov, float aspect, float near, float far);
    void lookAt( float(&matrix)[4][4], vec3f &eye, vec3f &lat, vec3f &up);
    
    void move(float(&matrix)[4][4], float x_axis, float y_axis, float z_axis);
    void rotate(float (&matrix)[4][4], float x_axis, float y_axis, float z_axis, float angle);
    void rotate(float (&matrix)[4][4], float x_angle, float y_angle, float z_angle);
    void rotate(mat3 &matrix, float x_angle, float y_angle, float z_angle);
    void scale(float(&matrix)[4][4], float x_axis, float y_axis, float z_axis);

    float angleDistancef(float a, float b);

    float getHeight(vec3f &trA, vec3f &trB, vec3f &trC, vec2f xyPos);
    vec3f cross(vec3f &a, vec3f &b); 
    inline float dot(const vec3f &a, const vec3f &b){ return a.x*b.x + a.y*b.y + a.z*b.z;} 
    vec3f normal3f(vec3f &trA, vec3f &trB, vec3f &trC);
    vec3f interpol(const vec3f &a, const vec3f &b, const vec3f &t);
    inline float interpol(float a, float b, float t) {return a * (1.0f - t) + b * t;}
    inline float interpolSquare(float a, float b, float t){return interpol(a, b, 3 * t * t - 2 * t * t * t);}
}