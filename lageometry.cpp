#include "lageometry.h"
#include <math.h>
#include <string.h>

void lag::fastMultiply(float (&mA)[4][4], float (&matB)[4][4], float (&mC)[4][4]){
    mC[0][0] = (mA[0][0] * matB[0][0]) + (mA[0][1] * matB[1][0]) + (mA[0][2] * matB[2][0]) + (mA[0][3] * matB[3][0]);
    mC[0][1] = (mA[0][0] * matB[0][1]) + (mA[0][1] * matB[1][1]) + (mA[0][2] * matB[2][1]) + (mA[0][3] * matB[3][1]);
    mC[0][2] = (mA[0][0] * matB[0][2]) + (mA[0][1] * matB[1][2]) + (mA[0][2] * matB[2][2]) + (mA[0][3] * matB[3][2]);
    mC[0][3] = (mA[0][0] * matB[0][3]) + (mA[0][1] * matB[1][3]) + (mA[0][2] * matB[2][3]) + (mA[0][3] * matB[3][3]);

    mC[1][0] = (mA[1][0] * matB[0][0]) + (mA[1][1] * matB[1][0]) + (mA[1][2] * matB[2][0]) + (mA[1][3] * matB[3][0]);
    mC[1][1] = (mA[1][0] * matB[0][1]) + (mA[1][1] * matB[1][1]) + (mA[1][2] * matB[2][1]) + (mA[1][3] * matB[3][1]);
    mC[1][2] = (mA[1][0] * matB[0][2]) + (mA[1][1] * matB[1][2]) + (mA[1][2] * matB[2][2]) + (mA[1][3] * matB[3][2]);
    mC[1][3] = (mA[1][0] * matB[0][3]) + (mA[1][1] * matB[1][3]) + (mA[1][2] * matB[2][3]) + (mA[1][3] * matB[3][3]);

    mC[2][0] = (mA[2][0] * matB[0][0]) + (mA[2][1] * matB[1][0]) + (mA[2][2] * matB[2][0]) + (mA[2][3] * matB[3][0]);
    mC[2][1] = (mA[2][0] * matB[0][1]) + (mA[2][1] * matB[1][1]) + (mA[2][2] * matB[2][1]) + (mA[2][3] * matB[3][1]);
    mC[2][2] = (mA[2][0] * matB[0][2]) + (mA[2][1] * matB[1][2]) + (mA[2][2] * matB[2][2]) + (mA[2][3] * matB[3][2]);
    mC[2][3] = (mA[2][0] * matB[0][3]) + (mA[2][1] * matB[1][3]) + (mA[2][2] * matB[2][3]) + (mA[2][3] * matB[3][3]);

    mC[3][0] = (mA[3][0] * matB[0][0]) + (mA[3][1] * matB[1][0]) + (mA[3][2] * matB[2][0]) + (mA[3][3] * matB[3][0]);
    mC[3][1] = (mA[3][0] * matB[0][1]) + (mA[3][1] * matB[1][1]) + (mA[3][2] * matB[2][1]) + (mA[3][3] * matB[3][1]);
    mC[3][2] = (mA[3][0] * matB[0][2]) + (mA[3][1] * matB[1][2]) + (mA[3][2] * matB[2][2]) + (mA[3][3] * matB[3][2]);
    mC[3][3] = (mA[3][0] * matB[0][3]) + (mA[3][1] * matB[1][3]) + (mA[3][2] * matB[2][3]) + (mA[3][3] * matB[3][3]);
}

void lag::fastMultiply(mat3 &mA, mat3 &mB, mat3 &mC)
{
    mC[0][0] = (mA[0][0] * mB[0][0]) + (mA[0][1] * mB[1][0]) + (mA[0][2] * mB[2][0]);
    mC[0][1] = (mA[0][0] * mB[0][1]) + (mA[0][1] * mB[1][1]) + (mA[0][2] * mB[2][1]);
    mC[0][2] = (mA[0][0] * mB[0][2]) + (mA[0][1] * mB[1][2]) + (mA[0][2] * mB[2][2]);

    mC[1][0] = (mA[1][0] * mB[0][0]) + (mA[1][1] * mB[1][0]) + (mA[1][2] * mB[2][0]);
    mC[1][1] = (mA[1][0] * mB[0][1]) + (mA[1][1] * mB[1][1]) + (mA[1][2] * mB[2][1]);
    mC[1][2] = (mA[1][0] * mB[0][2]) + (mA[1][1] * mB[1][2]) + (mA[1][2] * mB[2][2]);

    mC[2][0] = (mA[2][0] * mB[0][0]) + (mA[2][1] * mB[1][0]) + (mA[2][2] * mB[2][0]);
    mC[2][1] = (mA[2][0] * mB[0][1]) + (mA[2][1] * mB[1][1]) + (mA[2][2] * mB[2][1]);
    mC[2][2] = (mA[2][0] * mB[0][2]) + (mA[2][1] * mB[1][2]) + (mA[2][2] * mB[2][2]);
}

void lag::multiplyMx(const Dmatrixf &a, const Dmatrixf &b, Dmatrixf &c)
{
    if(b.h != a.w) return;

    if(c.w != b.w || c.h != a.h){
        c.reSize(b.w, a.h);
    }

    for(int i = 0; i < c.w; ++i){
        for(int j = 0; j < c.h; ++j){
            float sum = 0.0f;

            for(int k = 0; k < a.w; ++k){
                float av = a.getValue(k, i);
                float bv = b.getValue(j, k);
                sum += av * bv;
            }

            c.setValue(i, j, sum);
        }
    }
}

void lag::multiplyMxPartial(const Dmatrixf &a, const Dmatrixf &b, Dmatrixf &c, const MxChunk &chunk)
{
    if(b.h != a.w) return;

    if(c.w != b.w || c.h != a.h){
        c.reSize(b.w, a.h);
    }

    if(chunk.x_end >= c.w || chunk.y_end >= c.h) return;

    for(int i = chunk.x_start; i < chunk.x_end; ++i){
        for(int j = chunk.y_start; j < chunk.y_end; ++j){
            float sum = 0.0f;

            for(int k = 0; k < a.w; ++k){
                float av = a.getValue(k, i);
                float bv = b.getValue(j, k);
                sum += av * bv;
            }

            c.setValue(i, j, sum);
        }
    }
}

void lag::prontMatrix(const float (&matrix)[4][4])
{
    #ifdef DEBUG
        for(int i = 0; i < 4; i ++){
            printf("|");
            for(int j = 0; j < 4; j++){
                printf(" %f ", matrix[j][i]);
            }
            printf("|\n");
        }
    #endif
}

void lag::prontMatrix(const mat3 & matrix)
{
    //#ifdef DEBUG
        for(int i = 0; i < 3; i ++){
            printf("|");
            for(int j = 0; j < 3; j++){
                printf(" %f ", matrix[j][i]);
            }
            printf("|\n");
        }
    //#endif
}

void lag::prontMatrix(const Dmatrixf &matrix)
{
    for(int i = 0; i < matrix.w; ++i){
        printf("| ");
        for(int j = 0; j < matrix.h; ++j){
            printf(" %.1f ", matrix.getValue(i, j));
        }
        printf("|\n");
    }
}

void lag::makeIdentityMatrix(float (&matrix)[4][4])
{
    int k = 0; int j = 0;
    for(int i = 0; i < 16; i++){
        matrix[j][k] = 1.0 * (j == k);
        j += 1 * (k == 3);
        k = (k + 1) * (k != 3);
    }
}

void lag::makeIdentityMatrix(mat3 &matrix)
{
    matrix[0][0] = 1.0f; matrix[1][0] = 0.0f; matrix[2][0] = 0.0f;
    matrix[0][1] = 0.0f; matrix[1][1] = 1.0f; matrix[2][1] = 0.0f;
    matrix[0][2] = 0.0f; matrix[1][2] = 0.0f; matrix[2][2] = 1.0f;
}

void lag::project(float (&matrix)[4][4], float fov, float aspect, float near, float far)
{
    //#define LH_ZO 
    //#define LH_NO 
    //#define RH_ZO 
    #define RH_NO 
    float s = (tan(fov*0.5));

    #ifdef LH_ZO
    matrix[0][0] = 1.0/(aspect * s);  matrix[0][1] = 1;         matrix[0][2] = 0;                       matrix[0][3] = 0;
    matrix[1][0] = 0;                 matrix[1][1] = 1.0/s;     matrix[1][2] = 0;                       matrix[1][3] = 0;
    matrix[2][0] = 0;                 matrix[2][1] = 1;         matrix[2][2] = far/(far-near);          matrix[2][3] = 1.0;
    matrix[3][0] = 0;                 matrix[3][1] = 0;         matrix[3][2] = -(far*near)/(far-near);  matrix[3][3] = 0;
    #endif

    #ifdef LH_NO
    matrix[0][0] = 1.0/(aspect * s);  matrix[0][1] = 1;         matrix[0][2] = 0;                       matrix[0][3] = 0;
    matrix[1][0] = 0;                 matrix[1][1] = 1.0/s;     matrix[1][2] = 0;                       matrix[1][3] = 0;
    matrix[2][0] = 0;                 matrix[2][1] = 1;         matrix[2][2] = (far + near)/(far-near); matrix[2][3] = 1.0;
    matrix[3][0] = 0;                 matrix[3][1] = 0;         matrix[3][2] = -(2.0*far*near)/(far-near);  matrix[3][3] = 0;
    #endif

    #ifdef RH_ZO
    matrix[0][0] = 1.0/(aspect * s);  matrix[0][1] = 1;         matrix[0][2] = 0;                       matrix[0][3] = 0;
    matrix[1][0] = 0;                 matrix[1][1] = 1.0/s;     matrix[1][2] = 0;                       matrix[1][3] = 0;
    matrix[2][0] = 0;                 matrix[2][1] = 1;         matrix[2][2] = far/(near-far);          matrix[2][3] = -1.0;
    matrix[3][0] = 0;                 matrix[3][1] = 0;         matrix[3][2] = -(far*near)/(far-near);  matrix[3][3] = 0;
    #endif

    #ifdef RH_NO
    matrix[0][0] = 1.0/(aspect * s);  matrix[0][1] = 0;          matrix[0][2] = 0;                          matrix[0][3] = 0;
    matrix[1][0] = 0;                 matrix[1][1] = -1.0/s;     matrix[1][2] = 0;                          matrix[1][3] = 0;
    matrix[2][0] = 0;                 matrix[2][1] = 1;          matrix[2][2] = -(far + near)/(far-near);   matrix[2][3] = -1.0;
    matrix[3][0] = 0;                 matrix[3][1] = 0;          matrix[3][2] = -(2.0*far*near)/(far-near); matrix[3][3] = 0.0;
    #endif
}

void lag::move(float (&matrix)[4][4], float x_axis, float y_axis, float z_axis)
{
    matrix[3][0] += x_axis;
    matrix[3][1] += y_axis;
    matrix[3][2] += z_axis;

}

void lag::rotate(float (&matrix)[4][4], float x_angle, float y_angle, float z_angle){
    float sA = sin(x_angle); float cA = cos(x_angle);
    float sB = sin(y_angle); float cB = cos(y_angle);
    float sJ = sin(z_angle); float cJ = cos(z_angle);

    float resoult [4][4] = {};
    float tmp [4][4] = {};

    resoult[0][0] = cB*cJ;  resoult[1][0] = sA*sB*cJ - cA*sJ;   resoult[2][0] = cA*sB*cJ + sA*sJ;   resoult[3][0] = 0.0;
////////////////
    resoult[0][1] = cB*sJ;  resoult[1][1] = sA*sB*sJ + cA*cJ;   resoult[2][1] = cA*sB*sJ - sA*cJ;   resoult[3][1] = 0.0;
////////////////
    resoult[0][2] = -sB;    resoult[1][2] = sA*cB;              resoult[2][2] = cA*cB;              resoult[3][2] = 0.0;

    resoult[3][3] = 1.0;

    lag::fastMultiply(matrix, resoult, tmp);
    memcpy(matrix, tmp, sizeof(tmp));
}

void lag::rotate(mat3 &matrix, float x_angle, float y_angle, float z_angle)
{
    float sA = sin(x_angle); float cA = cos(x_angle);
    float sB = sin(y_angle); float cB = cos(y_angle);
    float sJ = sin(z_angle); float cJ = cos(z_angle);

    float resoult [3][3] = {};
    float tmp [3][3] = {};

    resoult[0][0] = cB*cJ;  resoult[1][0] = sA*sB*cJ - cA*sJ;   resoult[2][0] = cA*sB*cJ + sA*sJ;   //resoult[3][0] = 0.0;

    resoult[0][1] = cB*sJ;  resoult[1][1] = sA*sB*sJ + cA*cJ;   resoult[2][1] = cA*sB*sJ - sA*cJ;   //resoult[3][1] = 0.0;

    resoult[0][2] = -sB;    resoult[1][2] = sA*cB;              resoult[2][2] = cA*cB;              //resoult[3][2] = 0.0;

    //resoult[3][3] = 1.0;

    lag::fastMultiply(matrix, resoult, tmp);
    memcpy(matrix, tmp, sizeof(tmp));
}

void lag::rotate(float (&matrix)[4][4], float x_axis, float y_axis, float z_axis, float angle)
{
    //x rotation
    //static float x_multiplier[4][4] = {
    //    {1.0, 0.0,          0.0,         0.0},
    //    {0.0, cos(x_axis), -sin(x_axis), 0.0},
    //    {0.0, sin(x_axis), cos(x_axis),  0.0},   
    //    {0.0, 0.0,         0.0,          1.0}
    //};
//
    //static float y_multiplier[4][4] = {
    //    {cos(y_axis),   0.0,          sin(y_axis),  0.0},
    //    {0.0,           1.0,          0.0,          0.0},
    //    {-sin(y_axis),  0.0,          cos(y_axis),  0.0},       
    //    {0.0,           0.0,           0.0,          1.0}    
    //};
//
    //static float z_multiplier[4][4] = {
    //    {cos(z_axis),   -sin(z_axis),  0.0,          0.0},
    //    {sin(z_axis),   cos(z_axis),   0.0,          0.0},  
    //    {0.0,           0.0,           1.0,          0.0},      
    //    {0.0,           0.0,           0.0,          1.0}   
    //};
//
    //// X=X*Y*Z
    //float out_multiplier[4][4];
    //lag::fastMultiply(x_multiplier, y_multiplier, out_multiplier);
    //lag::fastMultiply(out_multiplier, z_multiplier, x_multiplier);
//
    //matrix[0][0] += x_multiplier[0][0];    matrix[1][0] += x_multiplier[1][0];    matrix[2][0] += x_multiplier[2][0];    matrix[3][0] += x_multiplier[3][0];
    //matrix[0][1] += x_multiplier[0][1];    matrix[1][1] += x_multiplier[1][1];    matrix[2][1] += x_multiplier[2][1];    matrix[3][1] += x_multiplier[3][1];
    //matrix[0][2] += x_multiplier[0][2];    matrix[1][2] += x_multiplier[1][2];    matrix[2][2] += x_multiplier[2][2];    matrix[3][2] += x_multiplier[3][2];

    float resoult [4][4] = {};
    float tmp [4][4] = {};
    double c_r = cos(angle);
    double d_r = sin(angle);

    resoult[0][0] = c_r + x_axis * x_axis * (1.0-c_r); //g
    resoult[1][0] =       y_axis * x_axis * (1.0-c_r) - z_axis * d_r; //g
    resoult[2][0] =       z_axis * x_axis * (1-c_r)   + y_axis * d_r; //g
    resoult[3][0] = 0.0;
////////////////
    resoult[0][1] = x_axis * y_axis * (1.0-c_r) + z_axis * d_r; //g
    resoult[1][1] = y_axis * y_axis * (1.0-c_r)          + c_r; //g
    resoult[2][1] = z_axis * y_axis * (1.0-c_r) - x_axis * d_r; //g
    resoult[3][1] = 0.0;
////////////////
    resoult[0][2] = x_axis * z_axis * (1.0-c_r) - y_axis * d_r; //g
    resoult[1][2] = y_axis * z_axis * (1.0-c_r) + x_axis * d_r; //g
    resoult[2][2] = z_axis * z_axis * (1.0-c_r)          + c_r;
    resoult[3][2] = 0.0;

    resoult[3][3] = 1.0;

    lag::fastMultiply(matrix, resoult, tmp);
    memcpy(matrix, tmp, sizeof(tmp));
}

lag::vec3f normalize(lag::vec3f &v){
    float mag = sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
    return {v.x/mag, v.y/mag, v.z/mag};
}

lag::vec3f cross(lag::vec3f &v, lag::vec3f &deez){
    lag::vec3f res = {};
    res.x = (v.y*deez.z - v.z*deez.y);
    res.y = (v.z*deez.x - v.x*deez.z);
    res.z = (v.x*deez.y - v.y*deez.x);
    return res;
}

void lag::lookAt( float(&matrix)[4][4], lag::vec3f &eye, lag::vec3f &lat, lag::vec3f &up){
    vec3f f = {lat.x-eye.x, lat.z-eye.z, lat.z-eye.z};
    f = normalize(f);

    vec3f s = cross(f, up);
    s = normalize(s);

    vec3f u = cross(s, f);

    matrix[0][0] = s.x;
	matrix[1][0] = s.y;
	matrix[2][0] = s.z;

	matrix[0][1] = u.x;
	matrix[1][1] = u.y;
	matrix[2][1] = u.z;
	
	matrix[0][2] =-f.x;
	matrix[1][2] =-f.y;
	matrix[2][2] =-f.z;

    matrix[3][0] = -(s.x*eye.x + s.y*eye.y + s.z*eye.z);
    matrix[3][1] = -(u.x*eye.x + u.y*eye.y + u.z*eye.z);
    matrix[3][2] =  (f.x*eye.x + f.y*eye.y + f.z*eye.z);
}

void lag::scale(float (&matrix)[4][4], float x_axis, float y_axis, float z_axis)
{
    matrix[0][0] *= x_axis;
    matrix[1][1] *= y_axis;
    matrix[2][2] *= z_axis;
}

float lag::angleDistancef(float a, float b)
{
    float retangle = b - a;
    retangle += M_PI*2.0f*(retangle < -M_PI);
    retangle -= M_PI*2.0f*(retangle > M_PI);
    return retangle;
}

float lag::getHeight(lag::vec3f &trA, lag::vec3f &trB, lag::vec3f &trC, lag::vec2f xyPos)
{
    vec3f n = normal3f(trA, trB, trC);

    float z = trA.z -(n.x * (xyPos.x - trA.x) + n.y * (xyPos.y - trA.y)) / n.z;
    return z;
}

lag::vec3f lag::cross(vec3f &a, vec3f &b)
{
    lag::vec3f s;
    s.x = a.y * b.z - a.z * b.y;
    s.y = a.z * b.x - a.x * b.z;
    s.z = a.x * b.y - a.y * b.x;
    return s;
}

lag::vec3f lag::normal3f(vec3f &trA, vec3f &trB, vec3f &trC)
{
    // points -> vectors tangent to the plane
    vec3f ab = trB - trA;
    vec3f ac = trC - trA;
 
    // now normal is cross cb x ac
    vec3f normal = cross(ac, ab);
    normal.normalize();
    return normal;
}

lag::vec3f lag::interpol(const vec3f &a, const vec3f &b, const vec3f &t)
{
   vec3f value;

   value.x = a.x * (1.0f - t.x) + b.x * (t.x);
   value.y = a.y * (1.0f - t.y) + b.y * (t.y);
   value.z = a.z * (1.0f - t.z) + b.z * (t.z);

   return value;
}

lag::PerlinNoice::PerlinNoice(unsigned int x_res, unsigned int y_res, bool normalized, bool randPositions)
{
    unsigned int pointArraySize = x_res * y_res;
    this->points = new vector3dPer [pointArraySize];

    this->width = x_res;
    this->height = y_res;

    float x_size = 1.0f / (float)x_res;
    float y_size = 1.0f / (float)y_res;

    for(int x = 0; x < x_res; x++){
        for(int y = 0; y < y_res; y++){
            vector3dPer & point = this->getPoint(x, y);

            if(randPositions){
                point.base.x = x_size * (float)x + x_size * drand48();
                point.base.y = y_size * (float)y + y_size * drand48();
            }else{
                point.base.x = x_size * (float)x;
                point.base.y = y_size * (float)y;
            }
            point.base.z = 0.0f;
            point.randomLength(normalized);
        }
    }
}

lag::PerlinNoice::~PerlinNoice()
{
    delete[] this->points;
}

lag::vector3dPer &lag::PerlinNoice::getPoint(unsigned int x, unsigned int y)
{
    x *= x < this->width;
    y *= y < this->height;
    
    unsigned int index = y * this->width + x;

    return this->points[index];
}

float lag::PerlinNoice::getValuef(float x, float y)
{
    unsigned int posx = x * (float)this->width;
    unsigned int posy = y * (float)this->height;

    vec3f point = {x, y, 0};

    vector3dPer & a = this->getPoint(posx,     posy);
    vector3dPer & b = this->getPoint(posx + 1, posy);
    vector3dPer & c = this->getPoint(posx,     posy + 1);
    vector3dPer & d = this->getPoint(posx + 1, posy + 1);

    vec3f pointNorm = {1.0f, 0.0f, 0.0f};
        pointNorm.rotareAlongY(20.0f * M_PI * 2.0f * x * y);

    float va = dot(a.getNormalized(), pointNorm);
    float vb = dot(b.getNormalized(), pointNorm);
    float vc = dot(c.getNormalized(), pointNorm);
    float vd = dot(d.getNormalized(), pointNorm);


    float cd_value = interpolSquare(vc, vd, (point.x - c.base.x) / (d.base.x - c.base.x));
    float ab_value = interpolSquare(va, vb, (point.x - a.base.x) / (b.base.x - a.base.x));

    float y_value = interpolSquare(ab_value, cd_value, (point.y - a.base.y) / (c.base.y - a.base.y));

    return y_value;
}

void lag::Dmatrixf::newMatrix(size_t width, size_t height)
{
    this->mxData = new float[width * height];

    this->w = width;
    this->h = height;
}

void lag::Dmatrixf::reSize(size_t width, size_t height)
{
    delete[] this->mxData;
    this->newMatrix(width, height);
}

float lag::Dmatrixf::getValue(unsigned int x, unsigned int y) const 
{
    if(x >= this->w || y >= this->h){
        return NAN;
    }
    
    return *(this->mxData + (this->w * y + x));
}

void lag::Dmatrixf::setValue(unsigned int x, unsigned int y, float v)
{
    if(x >= this->w || y >= this->h){
        return;
    }

    *(this->mxData + (this->w * y + x)) = v;
}

void lag::Dmatrixf::loadTestValues()
{
    float *beg = this->mxData;
    float *end = beg + this->w * this->h;

    float l = 1.0f, n = 1.0f;

    for(; beg < end; ++beg){
        *beg = l + n;
        l = n;
        n = *beg;

        if(n > 1e10){
            n -= 0.9f * n;
            l = n;
        }
    }
}
