
#include <windows.h>
#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include<bits/stdc++.h>
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<iostream>
#include <vector>
#include "bitmap_image.hpp"

#define pi (2*acos(0.0))
#define REFRACTING_INDEX 1.5
#define RAY_OFFSET 0.0001

///---------------init global vars-------------//
int sphere_slice = 30;
int sphere_stack = 30;
using namespace std;


struct point{
    double x;
    double y;
    double z;
    point(double a,double b,double c){
        x=a;
        y=b;
        z=c;
	}
	point(){}
	point operator-() const {
		return { -x, -y, -z };
	}

	point operator+(const point& p) const {
		return { x + p.x, y + p.y, z + p.z };
	}

	point operator-(const point& p) const {
		return { x - p.x, y - p.y, z - p.z };
	}

	point operator*(const point& p) const {
		return { x * p.x, y * p.y, z * p.z };
	}

	point operator*(const double constant) const {
		return { x * constant, y * constant, z * constant };
	}

	point operator/(const double constant) const {
		return { x / constant, y / constant, z / constant };
	}

	point operator+=(const point& plusplus) {
		return *this = *this + plusplus;
	}

	point operator-=(const point& minusminus) {
		return *this = *this - minusminus;
	}

	friend point operator+=(double constant, point const& p);
	friend point operator-=(double constant,point const& p);
	friend point operator*(double constant,point const& p);
	friend point operator/(double constant,point const& p);
};

inline point operator*(double constant, point const& p) {
	return p * constant;
}
inline point operator/(double constant, point const& p) {
	return p / constant;
}

///-------------------matrix operations--------------///
point unit_vector(point a)
{
    double VALUE=sqrt(a.x*a.x+a.y*a.y+a.z*a.z);
    struct point temp;
    temp.x=(a.x/VALUE);
    temp.y=(a.y/VALUE);
    temp.z=(a.z/VALUE);

    return temp;

}
double length_vector(point a)
{
    double VALUE=sqrt(a.x*a.x+a.y*a.y+a.z*a.z);


    return VALUE;

}
double dot_prod(point a,point b){
    double prod;
    prod=(a.x*b.x)+(a.y*b.y)+(a.z*b.z);
    return prod;
}
point cross_prod(point a,point b){
    struct point prod;
    prod.x=(a.y*b.z)-(b.y*a.z);
    prod.y=(a.z*b.x)-(a.x*b.z);
    prod.z=(a.x*b.y)-(a.y*b.x);
    return prod;
}
void matrix_mult(int temp,double mtx1[4][4],double mtx2[4][4]){
    double sum=0,mult_mtx[4][4];

    for(int i=0;i<4;i++){
        for(int j=0;j<temp;j++){
                mult_mtx[i][j]=0;
           for(int k=0;k<4;k++){
                mult_mtx[i][j]+=mtx1[i][k]*mtx2[k][j];
           }
        }
    }
    if(temp == 4){
        for(int i=0;i<4;i++){
            for(int j=0;j<4;j++){
                ///copy array
                //initial_mtx[i][j]=mult_mtx[i][j];
            }
        }

    }
    else{

        //pointVect.clear();
        ///push matrix
        //pointVect.push_back(mult_mtx);
    }


}
struct point negVector(struct point a){
    struct point temp;
    temp.x=-a.x;
    temp.y=-a.y;
    temp.z=-a.z;
    return temp;
}


///------------------------clip color----------------------///
void clipColor(double red,double green,double blue)
{
    if(red<0){
        red=0;
    }
    else if(red>1){
        red=1;
    }
    if(green<0){
        green=0;
    }
    else if(green>1){
        green=1;
    }
    if(blue<0){
        blue=0;
    }
    else if(blue>1){
        blue=1;
    }

}
///------------------------clip color----------------------///

///-------------------------ALL CLASSES-------------------///
class Ray
{
public:
    struct point start;
    struct point dir;
    Ray(){}
    Ray(struct point a,struct point b)
    {
        start=a;
        dir=b;
    }
};



class PointLight
{
public:
    struct point light_pos;
    double color[3];
    PointLight(){}
    PointLight(struct point light_pos){
        this->light_pos=light_pos;
    }
    void setColor(double r,double g,double b)
    {
        this->color[0]=r;
        this->color[1]=g;
        this->color[2]=b;
    }
    void draw()
    {
        glPointSize(5);
        glColor3f(color[0],color[1],color[2]);
        glBegin(GL_POINTS);
        {
            glVertex3f(light_pos.x,light_pos.y,light_pos.z);
        }glEnd();

    }
};
class SpotLight
{
public:
    PointLight point_light;
    struct point light_dir;
    double cutoff_angle;
    SpotLight(){}
    SpotLight(PointLight *point_light,struct point light_dir,double cutoff_angle){
        struct point pos(point_light->light_pos.x,point_light->light_pos.y,point_light->light_pos.z);
        this->point_light=PointLight(pos);
        this->point_light.setColor(point_light->color[0],point_light->color[1],point_light->color[2]);
        this->light_dir = light_dir;
        this->cutoff_angle = cutoff_angle;
    }
    void setLightDirection(double x,double y,double z)
    {
        this->light_dir.x=x;
        this->light_dir.y=y;
        this->light_dir.z=z;
    }
    void draw()
    {
        glPointSize(5);
        glColor3f(this->point_light.color[0],this->point_light.color[1],this->point_light.color[2]);
        glBegin(GL_POINTS);
        {
            glVertex3f(this->point_light.light_pos.x,this->point_light.light_pos.y,this->point_light.light_pos.z);
        }glEnd();

    }
};
struct point getRefractingPoint(Ray *r,struct point N)
{
    double product = dot_prod(r->dir,N);

    double k = 1 - (pow(REFRACTING_INDEX,2) * (1 - pow(product,2)));

    struct point temp_refracted_pt;
    if(k < 0){
        temp_refracted_pt.x = 0;
        temp_refracted_pt.y = 0;
        temp_refracted_pt.z = 0;
    }
    else{
        temp_refracted_pt.x = (REFRACTING_INDEX * r->dir.x) - (REFRACTING_INDEX * product + sqrt(k)) * N.x;
        temp_refracted_pt.y = (REFRACTING_INDEX * r->dir.y) - (REFRACTING_INDEX * product + sqrt(k)) * N.y;
        temp_refracted_pt.z = (REFRACTING_INDEX * r->dir.z) - (REFRACTING_INDEX * product + sqrt(k)) * N.z;
    }
    return unit_vector(temp_refracted_pt);
}
struct point getReflectingPoint(Ray *r,struct point N)
{
    struct point temp_reflected_pt;
    temp_reflected_pt.x = r->dir.x-N.x * 2.0 * dot_prod(r->dir,N);
    temp_reflected_pt.y = r->dir.y-N.y * 2.0 * dot_prod(r->dir,N);
    temp_reflected_pt.z = r->dir.z-N.z * 2.0 * dot_prod(r->dir,N);
    return unit_vector(temp_reflected_pt);
}
class Object
{
public:
    struct point reference_point;
    double height, width, length;

    double color[3];
    double coEfficients[4];// reflection coefficients
    int shine;// exponent term of specular component
    Object(){}
    virtual void draw(){}

    void setColor(double red,double green,double blue)
    {
        this->color[0]=red;
        this->color[1]=green;
        this->color[2]=blue;
    }
    void setShine(int shine)
    {
        this->shine=shine;
    }
    void setCoEfficients(double a,double b,double c,double d)
    {
        this->coEfficients[0]=a;
        this->coEfficients[1]=b;
        this->coEfficients[2]=c;
        this->coEfficients[3]=d;
    }

    void setCurrentColor(double *current_color)
    {
        current_color[0]=this->color[0]*this->coEfficients[0];
        current_color[1]=this->color[1]*this->coEfficients[0];
        current_color[2]=this->color[2]*this->coEfficients[0];
    }


    virtual struct point getUnitVect(struct point a)
	{
        return unit_vector(a);
	}
	struct point getRefractedRay(Ray *r, struct point N)
    {
        struct point refracted_pt;
        refracted_pt = getRefractingPoint(r,N);

        return refracted_pt;

    }
    struct point getReflectedRay(Ray *r,struct point N)
	{
	    struct point reflected_pt;
        reflected_pt = getReflectingPoint(r,N);
        return reflected_pt;
	}
	Ray getRay(PointLight *light,struct point intersectionPoint)
	{
	    struct point RAY1Start,RAY1Dir;
        RAY1Dir.x=light->light_pos.x-intersectionPoint.x;
        RAY1Dir.y=light->light_pos.y-intersectionPoint.y;
        RAY1Dir.z=light->light_pos.z-intersectionPoint.z;

        double lengthOfLDir=length_vector(RAY1Dir);
        RAY1Dir=unit_vector(RAY1Dir);

        RAY1Start.x=intersectionPoint.x+(RAY1Dir.x*RAY_OFFSET);
        RAY1Start.y=intersectionPoint.y+(RAY1Dir.y*RAY_OFFSET);
        RAY1Start.z=intersectionPoint.z+(RAY1Dir.z*RAY_OFFSET);
        Ray RAY1(RAY1Start,RAY1Dir);
        return RAY1;
	}
	Ray getRay2(PointLight *light,struct point intersectionPoint)
	{
	    struct point RAY1Start,RAY1Dir;
        RAY1Dir.x=intersectionPoint.x-light->light_pos.x;
        RAY1Dir.y=intersectionPoint.y-light->light_pos.y;
        RAY1Dir.z=intersectionPoint.z-light->light_pos.z;

        double lengthOfLDir=length_vector(RAY1Dir);
        RAY1Dir=unit_vector(RAY1Dir);

        RAY1Start.x=intersectionPoint.x+(RAY1Dir.x*RAY_OFFSET);
        RAY1Start.y=intersectionPoint.y+(RAY1Dir.y*RAY_OFFSET);
        RAY1Start.z=intersectionPoint.z+(RAY1Dir.z*RAY_OFFSET);
        Ray RAY1(RAY1Start,RAY1Dir);
        return RAY1;
	}
	Ray getRecursiveRay(Ray *r,struct point normal_at_intersect,struct point intersectionPoint)
	{
	    struct point ReflectionRAY1Dir,ReflectionRAY1Start;
        ReflectionRAY1Dir=getReflectedRay(r,normal_at_intersect);


        ReflectionRAY1Start.x=intersectionPoint.x+ReflectionRAY1Dir.x*RAY_OFFSET;
        ReflectionRAY1Start.y=intersectionPoint.y+ReflectionRAY1Dir.y*RAY_OFFSET;
        ReflectionRAY1Start.z=intersectionPoint.z+ReflectionRAY1Dir.z*RAY_OFFSET;

        Ray ReflectedRay(ReflectionRAY1Start,ReflectionRAY1Dir);
        return ReflectedRay;
	}
	double getLenOfRAYDir(PointLight *light,struct point intersectionPoint)
	{
	    struct point RAY1Dir;
        RAY1Dir.x=light->light_pos.x-intersectionPoint.x;
        RAY1Dir.y=light->light_pos.y-intersectionPoint.y;
        RAY1Dir.z=light->light_pos.z-intersectionPoint.z;

        double lengthOfDir=length_vector(RAY1Dir);
        return lengthOfDir;
	}

    virtual double get_tMIN(Ray *r)
    {
        return -1;
	}
	virtual double intersect(Ray *r, double *current_color, int level){
        return -1;
	}



};



extern string PATH;
extern int recursion_level;
extern vector<PointLight*> point_lights;
extern vector<SpotLight*> spot_lights;
extern vector<Object*> objects;


///---------------------HELPING FUNCTIONS-------------------.///



int calculateNearest(Ray *ReflectedRay)
{
    int nearest=-1;
    int tMin=999999;



    double *dummyColor=new double[3];
    for(int k=0;k<objects.size();k++){
        double temp=objects[k]->intersect(ReflectedRay,dummyColor,0);
        if(temp<=0){
            continue;
        }
        else if(temp<tMin){
            tMin=temp;
            nearest=k;


        }
    }
    delete[] dummyColor;
    return nearest;

}
int calculatetMIN(Ray *ReflectedRay)
{
    int tMin=999999;



    double *dummyColor=new double[3];
    for(int k=0;k<objects.size();k++){
        double temp=objects[k]->intersect(ReflectedRay,dummyColor,0);
        if(temp<=0){
            continue;
        }
        else if(temp<tMin){
            tMin=temp;



        }
    }
    delete[] dummyColor;
    return tMin;
}
///---------------------HELPING FUNCTIONS-------------------.///

class Sphere:public Object
{
public:
    Sphere(struct point centre,double radius){
        reference_point=centre;
        length=radius;
        height=radius;
    }

    double get_tMIN(Ray *r)
    {
        struct point p;
        p.x = r->start.x-reference_point.x;
        p.y = r->start.y-reference_point.y;
        p.z = r->start.z-reference_point.z;

        double AX_2=dot_prod(r->dir,r->dir);
        double BX=2.0*(dot_prod(r->dir,p));
        double C_=(dot_prod(p,p)-(pow(length,2)));
        double DET=pow(BX,2)-(4*AX_2*C_);
        if(DET<0)
        {
            return -1;
        }
        ///ELSE TWO CASES
        ///CASE1:
        double t1=(-BX+sqrt(DET))/(2.0*AX_2);
        ///CASE2:
        double t2=(-BX-sqrt(DET))/(2.0*AX_2);
        if(t1<=0 && t2<=0)
        {
            return -1;
        }
        else if(t1<=0 && t2>0)
        {
            return t2;
        }
        else if(t1>0 && t2<=0)
        {
            return t1;
        }
        else
        {
            if(t1<t2) return t1;
            else return t2;

        }

    }
    bool isObscured(Ray *r,Ray *RAY1,double Dir_len)
    {
        bool obscured=false;

        ///itr through all objects
        for(int j=0;j<objects.size();j++){
            struct point temp_point2;
            temp_point2.x = r->start.x - reference_point.x;
            temp_point2.y = r->start.y - reference_point.y;
            temp_point2.z = r->start.z - reference_point.z;
            double tr=objects[j]->get_tMIN(RAY1);
            if(tr>0 && abs(tr)<Dir_len){
                obscured=true;
                break;
            }
        }
        ///here write codes
        return obscured;
    }


    double intersect(Ray *r, double *current_color, int level){
        struct point temp_point;
        temp_point.x = r->start.x - reference_point.x;
        temp_point.y = r->start.y - reference_point.y;
        temp_point.z = r->start.z - reference_point.z;

        ///CALCULATE TMIN
        double tMIN=get_tMIN(r);
        if(tMIN<=0){
            return -1;
        }
        ///if level is 0, return tmin
        if(level==0){
            return tMIN;
        }
        struct point intersectionPoint(r->start.x+(r->dir.x*tMIN),r->start.y+(r->dir.y*tMIN),r->start.z+(r->dir.z*tMIN));
        ///double intersectionPointColor[3];
        ///setCurrentColor(current_color);


        current_color[0]=this->color[0]*this->coEfficients[0];
        current_color[1]=this->color[1]*this->coEfficients[0];
        current_color[2]=this->color[2]*this->coEfficients[0];

        ///calculate normal at intersect
        struct point normal_at_intersect,temp;

        temp.x=intersectionPoint.x-reference_point.x;
        temp.y=intersectionPoint.y-reference_point.y;
        temp.z=intersectionPoint.z-reference_point.z;
        normal_at_intersect=unit_vector(temp);


        for(int i=0;i<point_lights.size();i++){

            ///cast rayl from pl.light_pos to intersectionPoint

            Ray RAY1;
            RAY1 = getRay(point_lights[i],intersectionPoint);
            double lengthOfLDir = getLenOfRAYDir(point_lights[i],intersectionPoint);

            ///check if obscured
            bool obscured;

            ///IF NOT OBSCURED
            if(isObscured(r,&RAY1,lengthOfLDir)==false){
                double lambertValue=dot_prod(RAY1.dir,normal_at_intersect);

                double phongValue=dot_prod(getReflectedRay(&RAY1,normal_at_intersect),r->dir);
                if(lambertValue <= 0) lambertValue = 0;
                if(phongValue <= 0) phongValue = 0;


                current_color[0]+=(point_lights[i]->color[0]*lambertValue*this->color[0]*this->coEfficients[1])+(point_lights[i]->color[0]*pow(phongValue,this->shine)*this->coEfficients[2]);//*this->color[0]);
                current_color[1]+=(point_lights[i]->color[1]*lambertValue*this->color[1]*this->coEfficients[1])+(point_lights[i]->color[1]*pow(phongValue,this->shine)*this->coEfficients[2]);//*this->color[1]);
                current_color[2]+=(point_lights[i]->color[2]*lambertValue*this->color[2]*this->coEfficients[1])+(point_lights[i]->color[2]*pow(phongValue,this->shine)*this->coEfficients[2]);//*this->color[2]);
            }
        }

        ///code for spotlight
        for(int i=0;i<spot_lights.size();i++){
            PointLight pl = spot_lights[i]->point_light;

            ///cast rayl from pl.light_pos to intersectionPoint

            Ray RAY1,RAY2;
            RAY2 = getRay(&spot_lights[i]->point_light,intersectionPoint);
            RAY1 = getRay2(&spot_lights[i]->point_light,intersectionPoint);  ///get line from intersect to spotlight
            double lengthOfLDir = getLenOfRAYDir(&spot_lights[i]->point_light,intersectionPoint);
            struct point LineFactor = unit_vector(RAY1.dir);   ///get unitvector of line
            double spotFactor = dot_prod(LineFactor,spot_lights[i]->light_dir);  ///dotprod with line dir to get cos(theta)
            double OBJ_angle = (acos(spotFactor)*180)/pi;  ///get theta in degree




            ///check if obscured
            bool obscured;

            ///IF NOT OBSCURED
            if(isObscured(r,&RAY2,lengthOfLDir)==false && OBJ_angle < spot_lights[i]->cutoff_angle){   ///compare with cutoff_angle
                double lambertValue=dot_prod(RAY2.dir,normal_at_intersect);

                double phongValue=dot_prod(getReflectedRay(&RAY2,normal_at_intersect),r->dir);
                if(lambertValue <= 0) lambertValue = 0;
                if(phongValue <= 0) phongValue = 0;


                current_color[0]+=(spot_lights[i]->point_light.color[0]*lambertValue*this->color[0]*this->coEfficients[1])+(spot_lights[i]->point_light.color[0]*pow(phongValue,this->shine)*this->coEfficients[2]);//*this->color[0]);
                current_color[1]+=(spot_lights[i]->point_light.color[1]*lambertValue*this->color[1]*this->coEfficients[1])+(spot_lights[i]->point_light.color[1]*pow(phongValue,this->shine)*this->coEfficients[2]);//*this->color[1]);
                current_color[2]+=(spot_lights[i]->point_light.color[2]*lambertValue*this->color[2]*this->coEfficients[1])+(spot_lights[i]->point_light.color[2]*pow(phongValue,this->shine)*this->coEfficients[2]);//*this->color[2]);
            }
        }



        ///reflection
        if(level>=recursion_level) return tMIN;
        ///ELSE RECURSION CALL
        if(level<recursion_level){



            Ray ReflectedRay;
            ReflectedRay = getRecursiveRay(r,normal_at_intersect,intersectionPoint);


            int nearest;
            double tMin;


            nearest = calculateNearest(&ReflectedRay);
            tMin = calculatetMIN(&ReflectedRay);

            //double tMin=999999;
            double *reflectedColor=new double[3];

            if(nearest!=-1){
                 tMin=objects[nearest]->intersect(&ReflectedRay,reflectedColor,level+1);
                 if(tMin!=-1){
                    current_color[0]+=reflectedColor[0]*this->coEfficients[3];
                    current_color[1]+=reflectedColor[1]*this->coEfficients[3];
                    current_color[2]+=reflectedColor[2]*this->coEfficients[3];
                }

            }

            clipColor(current_color[0],current_color[1],current_color[2]);

            delete[] reflectedColor;

            return tMin;
        }


    }

    void draw(){
        glColor3f(color[0], color[1], color[2]);

        vector<vector<point>> points(sphere_stack + 2, vector<point>(sphere_slice + 2));
        int i,j;
        double h,r;

        /// points
        for(i=0;i<=sphere_stack;i++)
        {
            h=length*sin(((double)i/(double)sphere_stack)*(pi/2.0));
            r=length*cos(((double)i/(double)sphere_stack)*(pi/2.0));
            for(j=0;j<=sphere_slice;j++)
            {
                points[i][j].x=r*cos(((double)j/(double)sphere_slice)*2.0*pi);
                points[i][j].y=r*sin(((double)j/(double)sphere_slice)*2.0*pi);
                points[i][j].z=h;
            }
        }
        ///draw quads
        for(i=0;i<sphere_stack;i++)
        {
            for(j=0;j<sphere_slice;j++)
            {
                glBegin(GL_QUADS);{
                    ///upper hemisphere
                    glVertex3f(points[i][j].x + reference_point.x, points[i][j].y + reference_point.y, points[i][j].z + reference_point.z);
                    glVertex3f(points[i][j+1].x + reference_point.x, points[i][j+1].y + reference_point.y, points[i][j+1].z + reference_point.z);
                    glVertex3f(points[i+1][j+1].x + reference_point.x, points[i+1][j+1].y + reference_point.y, points[i+1][j+1].z + reference_point.z);
                    glVertex3f(points[i+1][j].x + reference_point.x, points[i+1][j].y + reference_point.y, points[i+1][j].z + reference_point.z);
                    ///lower hemisphere
                    glVertex3f(points[i][j].x + reference_point.x,points[i][j].y + reference_point.y,-points[i][j].z + reference_point.z);
                    glVertex3f(points[i][j+1].x + reference_point.x,points[i][j+1].y + reference_point.y,-points[i][j+1].z + reference_point.z);
                    glVertex3f(points[i+1][j+1].x + reference_point.x,points[i+1][j+1].y + reference_point.y,-points[i+1][j+1].z + reference_point.z);
                    glVertex3f(points[i+1][j].x + reference_point.x,points[i+1][j].y + reference_point.y,-points[i+1][j].z + reference_point.z);
                }glEnd();
            }
        }
    }



};

///---------------------floor class-----------------///
class Floor:public Object
{
public:
    Floor(double floorWidth,double tileWidth){
        reference_point.x=(-1.0)*(floorWidth/2);
        reference_point.y=(-1.0)*(floorWidth/2);
        reference_point.z=0;
        length=tileWidth;

    }
    double get_tMIN(Ray *r)
    {
        struct point normal;
        normal.x=0;
        normal.y=0;
        normal.z=1;
        if(r->dir.z==0){
            return -1;
        }
        return (-1.0)*(dot_prod(normal, r->start) / (dot_prod(normal, r->dir)));

    }
    bool isObscured(Ray *RAY1,double Dir_len)
    {
        bool obscured=false;

        ///itr through all objects
        for(int j=0;j<objects.size();j++){
            struct point normal;
            normal.x=0;
            normal.y=0;
            normal.z=1;

            double tr=objects[j]->get_tMIN(RAY1);
            if(tr>0 && abs(tr)<Dir_len){
                //cout<<obscured<<endl;
                obscured=true;
                return obscured;
                //break;
            }
        }
        //cout<<obscured<<endl;
        ///here write codes
        return obscured;
    }



    double intersect(Ray *r, double *current_color, int level){
        struct point p_N(0,0,1);

        double tMIN=get_tMIN(r);
        if(tMIN<=0){
            return -1;
        }
        ///next day start from here---------
        struct point intersectionPoint(r->start.x+(r->dir.x*tMIN),r->start.y+(r->dir.y*tMIN),r->start.z+(r->dir.z*tMIN));


        if(intersectionPoint.x < reference_point.x || intersectionPoint.x > -reference_point.x || intersectionPoint.y < reference_point.y || intersectionPoint.y > -reference_point.y){
            return -1;
        }
        if(level==0){
            return tMIN;
        }

        int difX = (int)((intersectionPoint.x - reference_point.x) / length);
        int dirY = (int)((intersectionPoint.y - reference_point.y) / length);



        if((difX + dirY) % 2 == 0){
            color[0] = color[1] = color[2] = 1;
        }
        else {
            color[0] = color[1] = color[2] = 0;

        }
        ///setCurrentColor(current_color);
        current_color[0]=this->color[0]*this->coEfficients[0];
        current_color[1]=this->color[1]*this->coEfficients[0];
        current_color[2]=this->color[2]*this->coEfficients[0];

        struct point normal_at_intersect(0,0,1);

        for(int i=0;i<point_lights.size();i++){
            Ray RAY1;
            RAY1 = getRay(point_lights[i],intersectionPoint);
            double lengthOfLDir = getLenOfRAYDir(point_lights[i],intersectionPoint);


            bool obscured;


            if(isObscured(&RAY1,lengthOfLDir)==false){

                double lambertValue=dot_prod(RAY1.dir,normal_at_intersect);

                double phongValue=dot_prod(getReflectedRay(&RAY1,normal_at_intersect),r->dir);

                if(lambertValue <= 0) lambertValue = 0;
                if(phongValue <= 0) phongValue = 0;
                current_color[0]+=(point_lights[i]->color[0]*lambertValue*this->color[0]*this->coEfficients[1])+(point_lights[i]->color[0]*pow(phongValue,this->shine)*this->coEfficients[2]);
                current_color[1]+=(point_lights[i]->color[1]*lambertValue*this->color[1]*this->coEfficients[1])+(point_lights[i]->color[1]*pow(phongValue,this->shine)*this->coEfficients[2]);
                current_color[2]+=(point_lights[i]->color[2]*lambertValue*this->color[2]*this->coEfficients[1])+(point_lights[i]->color[2]*pow(phongValue,this->shine)*this->coEfficients[2]);
            }
        }
        ///code for spotlight
        for(int i=0;i<spot_lights.size();i++){
            PointLight pl = spot_lights[i]->point_light;

            ///cast rayl from pl.light_pos to intersectionPoint

            Ray RAY1,RAY2;
            RAY2 = getRay(&spot_lights[i]->point_light,intersectionPoint);
            RAY1 = getRay2(&spot_lights[i]->point_light,intersectionPoint);  ///get line from intersect to spotlight
            double lengthOfLDir = getLenOfRAYDir(&spot_lights[i]->point_light,intersectionPoint);
            struct point LineFactor = unit_vector(RAY1.dir);   ///get unitvector of line
            double spotFactor = dot_prod(LineFactor,spot_lights[i]->light_dir);  ///dotprod with line dir to get cos(theta)
            double OBJ_angle = (acos(spotFactor)*180)/pi;  ///get theta in degree




            ///check if obscured
            bool obscured;

            ///IF NOT OBSCURED
            if(isObscured(&RAY2,lengthOfLDir)==false && OBJ_angle < spot_lights[i]->cutoff_angle){   ///compare with cutoff_angle
                double lambertValue=dot_prod(RAY2.dir,normal_at_intersect);

                double phongValue=dot_prod(getReflectedRay(&RAY2,normal_at_intersect),r->dir);
                if(lambertValue <= 0) lambertValue = 0;
                if(phongValue <= 0) phongValue = 0;


                current_color[0]+=(spot_lights[i]->point_light.color[0]*lambertValue*this->color[0]*this->coEfficients[1])+(spot_lights[i]->point_light.color[0]*pow(phongValue,this->shine)*this->coEfficients[2]);//*this->color[0]);
                current_color[1]+=(spot_lights[i]->point_light.color[1]*lambertValue*this->color[1]*this->coEfficients[1])+(spot_lights[i]->point_light.color[1]*pow(phongValue,this->shine)*this->coEfficients[2]);//*this->color[1]);
                current_color[2]+=(spot_lights[i]->point_light.color[2]*lambertValue*this->color[2]*this->coEfficients[1])+(spot_lights[i]->point_light.color[2]*pow(phongValue,this->shine)*this->coEfficients[2]);//*this->color[2]);
            }
        }

        if(level >= recursion_level) return tMIN;
        ///ELSE RECURSION CALL
        if(level < recursion_level){
            Ray ReflectedRay;
            ReflectedRay = getRecursiveRay(r,normal_at_intersect,intersectionPoint);

            int nearest;
            double tMin;


            nearest = calculateNearest(&ReflectedRay);
            tMin = calculatetMIN(&ReflectedRay);
            double *reflectedColor=new double[3];

            if(nearest!=-1){
                 tMin=objects[nearest]->intersect(&ReflectedRay,reflectedColor,level+1);
                 if(tMin!=-1){
                    current_color[0]+=reflectedColor[0]*this->coEfficients[3];
                    current_color[1]+=reflectedColor[1]*this->coEfficients[3];
                    current_color[2]+=reflectedColor[2]*this->coEfficients[3];
                }

            }

            clipColor(current_color[0],current_color[1],current_color[2]);

            delete[] reflectedColor;

            return tMin;


        }

    }
    void draw(){

        double refX,refY,refZ;
        refX = reference_point.x;
        refY = reference_point.y;
        refZ = reference_point.z;

        int toggle = 0;
        for(refY = reference_point.y; refY < -reference_point.y; refY += length)
        {
            for(refX=reference_point.x; refX < -reference_point.x; refX += length)
            {
                if(toggle == 1)
                {
                    glColor3f(0, 0, 0);
                }
                else
                {
                    glColor3f(1, 1, 1);
                }


                glBegin(GL_QUADS);
                {
                    glVertex3f(refX, refY, refZ);
                    glVertex3f(refX, refY + length, refZ);
                    glVertex3f(refX + length, refY + length, refZ);
                    glVertex3f(refX + length, refY, refZ);
                }
                glEnd();
                toggle =  toggle ^ 1;

            }

            toggle = toggle ^ 1;
        }
    }
};
///---------------------floor class-----------------///




///---------------------------triangle class----------------------///
class Triangle:public Object
{
public:
    struct point a,b,c;

    Triangle(){}
    Triangle(struct point a,struct point b,struct point c){
        this->a=a;
        this->b=b;
        this->c=c;
    }
    void draw()
    {
        glColor3f(color[0],color[1],color[2]);
        glBegin(GL_TRIANGLES);
        {
            glVertex3f(a.x,a.y,a.z);
			glVertex3f(b.x,b.y,b.z);
			glVertex3f(c.x,c.y,c.z);
        }
        glEnd();
    }
    bool isObscured(Ray *RAY1,double Dir_len)
    {
        bool obscured=false;

        ///itr through all objects
        for(int j=0;j<objects.size();j++){
            //struct point normal;
            //normal.x=0;
            //normal.y=0;
            //normal.z=1;

            double tr=objects[j]->get_tMIN(RAY1);
            if(tr>0 && abs(tr)<Dir_len){
                //cout<<obscured<<endl;
                obscured=true;
                return obscured;
                //break;
            }
        }
        //cout<<obscured<<endl;
        ///here write codes
        return obscured;
    }

    double get_tMIN(Ray *r)
    {
        struct point AB(a.x-b.x,a.y-b.y,a.z-b.z);
        struct point AC(a.x-c.x,a.y-c.y,a.z-c.z);
        struct point AR(a.x-r->start.x,a.y-r->start.y,a.z-r->start.z);


        double A=0.0,Gamma=0.0,Beta=0.0,t=0.0;

        ///DETERMINANT VALUE

        A = A + (AB.x*(AC.y*r->dir.z - r->dir.y*AC.z));

        A = A - (AB.y*(AC.x*r->dir.z - r->dir.x*AC.z));

        A = A + (AB.z*(AC.x*r->dir.y - r->dir.x*AC.y));

        if(A==0)
        {
            return -1;
        }



        ///DETERMINANT VALUE

        Gamma = Gamma + (AB.x*(AR.y*r->dir.z - r->dir.y*AR.z));

        Gamma = Gamma - (AB.y*(AR.x*r->dir.z - r->dir.x*AR.z));

        Gamma = Gamma + (AB.z*(AR.x*r->dir.y - r->dir.x*AR.y));
        Gamma = Gamma / A;

        if(Gamma<0 || Gamma>1.0)
        {
            return -1;
        }



        ///DETERMINANT VALUE

        Beta = Beta + (AR.x*(AC.y*r->dir.z - r->dir.y*AC.z));

        Beta = Beta - (AR.y*(AC.x*r->dir.z - r->dir.x*AC.z));

        Beta = Beta + (AR.z*(AC.x*r->dir.y - r->dir.x*AC.y));

        Beta = Beta / A;
        //cout<<Beta<<endl;
        if(Beta < 0 || Beta > 1-Gamma)
        {
            return -1;
        }

        ///determinant

        t = t + (AB.x*(AC.y*AR.z-AR.y*AC.z));
        t = t - (AB.y*(AC.x*AR.z-AR.x*AC.z));
        t = t + (AB.z*(AC.x*AR.y-AR.x*AC.y));
        t = t / A;

        if(t<0){
            return -1;
        }
        return t;

    }
    double intersect(Ray *r, double *current_color, int level){
        double tMIN=get_tMIN(r);
        if(tMIN<=0){
            return -1;
        }
        struct point intersectionPoint(r->start.x+(r->dir.x*tMIN),r->start.y+(r->dir.y*tMIN),r->start.z+(r->dir.z*tMIN));

        if(level==0){
            return tMIN;
        }

        ///setCurrentColor(current_color);
        current_color[0]=this->color[0]*this->coEfficients[0];
        current_color[1]=this->color[1]*this->coEfficients[0];
        current_color[2]=this->color[2]*this->coEfficients[0];
        ///-----------------------NORMAL CALC--------------------///
        ///calc ab,ac
        struct point AB,AC;
        struct point n;
        AB.x = b.x - a.x;
        AB.y = b.y - a.y;
        AB.z = b.z - a.z;

        AC.x = c.x - a.x;
        AC.y = c.y - a.y;
        AC.z = c.z - a.z;

        ///calculate normal at intersect
        struct point normal_at_intersect,temp;
        temp = cross_prod(AB,AC);


        normal_at_intersect=unit_vector(temp);


        if(dot_prod(intersectionPoint, normal_at_intersect) > 0){
            normal_at_intersect.x = -normal_at_intersect.x;
            normal_at_intersect.y = -normal_at_intersect.y;
            normal_at_intersect.z = -normal_at_intersect.z;
        }
        ///------------------------END OF NORMAL CALC----------------------///

        for(int i=0;i<point_lights.size();i++){
            Ray RAY1;
            RAY1 = getRay(point_lights[i],intersectionPoint);
            double lengthOfLDir = getLenOfRAYDir(point_lights[i],intersectionPoint);


            bool obscured;


            if(isObscured(&RAY1,lengthOfLDir)==false){

                double lambertValue=dot_prod(RAY1.dir,normal_at_intersect);

                double phongValue=dot_prod(getReflectedRay(&RAY1,normal_at_intersect),r->dir);

                if(lambertValue <= 0) lambertValue = 0;
                if(phongValue <= 0) phongValue = 0;
                current_color[0]+=(point_lights[i]->color[0]*lambertValue*this->color[0]*this->coEfficients[1])+(point_lights[i]->color[0]*pow(phongValue,this->shine)*this->coEfficients[2]);
                current_color[1]+=(point_lights[i]->color[1]*lambertValue*this->color[1]*this->coEfficients[1])+(point_lights[i]->color[1]*pow(phongValue,this->shine)*this->coEfficients[2]);
                current_color[2]+=(point_lights[i]->color[2]*lambertValue*this->color[2]*this->coEfficients[1])+(point_lights[i]->color[2]*pow(phongValue,this->shine)*this->coEfficients[2]);

            }
        }
        ///code for spotlight
        for(int i=0;i<spot_lights.size();i++){
            PointLight pl = spot_lights[i]->point_light;

            ///cast rayl from pl.light_pos to intersectionPoint

            Ray RAY1,RAY2;
            RAY2 = getRay(&spot_lights[i]->point_light,intersectionPoint);
            RAY1 = getRay2(&spot_lights[i]->point_light,intersectionPoint);  ///get line from intersect to spotlight
            double lengthOfLDir = getLenOfRAYDir(&spot_lights[i]->point_light,intersectionPoint);
            struct point LineFactor = unit_vector(RAY1.dir);   ///get unitvector of line
            double spotFactor = dot_prod(LineFactor,spot_lights[i]->light_dir);  ///dotprod with line dir to get cos(theta)
            double OBJ_angle = (acos(spotFactor)*180)/pi;  ///get theta in degree




            ///check if obscured
            bool obscured;

            ///IF NOT OBSCURED
            if(isObscured(&RAY2,lengthOfLDir)==false && OBJ_angle < spot_lights[i]->cutoff_angle){   ///compare with cutoff_angle
                double lambertValue=dot_prod(RAY2.dir,normal_at_intersect);

                double phongValue=dot_prod(getReflectedRay(&RAY2,normal_at_intersect),r->dir);
                if(lambertValue <= 0) lambertValue = 0;
                if(phongValue <= 0) phongValue = 0;


                current_color[0]+=(spot_lights[i]->point_light.color[0]*lambertValue*this->color[0]*this->coEfficients[1])+(spot_lights[i]->point_light.color[0]*pow(phongValue,this->shine)*this->coEfficients[2]);//*this->color[0]);
                current_color[1]+=(spot_lights[i]->point_light.color[1]*lambertValue*this->color[1]*this->coEfficients[1])+(spot_lights[i]->point_light.color[1]*pow(phongValue,this->shine)*this->coEfficients[2]);//*this->color[1]);
                current_color[2]+=(spot_lights[i]->point_light.color[2]*lambertValue*this->color[2]*this->coEfficients[1])+(spot_lights[i]->point_light.color[2]*pow(phongValue,this->shine)*this->coEfficients[2]);//*this->color[2]);
            }
        }

       //reflection
       if(level >= recursion_level) return tMIN;

        if(level < recursion_level){
            Ray ReflectedRay;
            ReflectedRay = getRecursiveRay(r,normal_at_intersect,intersectionPoint);

            int nearest;
            double tMin;


            nearest = calculateNearest(&ReflectedRay);
            tMin = calculatetMIN(&ReflectedRay);
            double *reflectedColor=new double[3];

            reflectedColor[0]=reflectedColor[1]=reflectedColor[2]=0.0;

            if(nearest!=-1){
                 tMin=objects[nearest]->intersect(&ReflectedRay,reflectedColor,level+1);
                 if(tMin!=-1){
                    current_color[0]+=reflectedColor[0]*this->coEfficients[3];
                    current_color[1]+=reflectedColor[1]*this->coEfficients[3];
                    current_color[2]+=reflectedColor[2]*this->coEfficients[3];
                }

            }

            clipColor(current_color[0],current_color[1],current_color[2]);

            delete[] reflectedColor;
            return tMin;
        }




    }
};

///---------------------------triangle class---------------------///


///--------------------generic class ------------------------------///
class GeneralQuadratic:public Object
{
public:
    double A,B,C,D,E,F,G,H,I,J;

    GeneralQuadratic(){}
    GeneralQuadratic(double A,double B,double C,double D,double E,double F,double G,double H,double I,double J,double length,double width,double height,struct point reference_point){
        this->A=A;
        this->B=B;
        this->C=C;
        this->D=D;
        this->E=E;
        this->F=F;
        this->G=G;
        this->H=H;
        this->I=I;
        this->J=J;
        this->reference_point=reference_point;
        this->length=length;
        this->width=width;
        this->height=height;
    }
    void setParams(double A,double B,double C,double D,double E,double F,double G,double H,double I,double J)
    {
        this->A=A;
        this->B=B;
        this->C=C;
        this->D=D;
        this->E=E;
        this->F=F;
        this->G=G;
        this->H=H;
        this->I=I;
        this->J=J;
    }
    void draw()
    {
        /*glColor3f(color[0],color[1],color[2]);
        glBegin(GL_TRIANGLES);
        {
            glVertex3f(a.x,a.y,a.z);
			glVertex3f(b.x,b.y,b.z);
			glVertex3f(c.x,c.y,c.z);
        }
        glEnd();*/
    }
    bool isObscured(Ray *RAY1,double Dir_len)
    {
        bool obscured=false;

        ///itr through all objects
        for(int j=0;j<objects.size();j++){
            //struct point normal;
            //normal.x=0;
            //normal.y=0;
            //normal.z=1;

            double tr=objects[j]->get_tMIN(RAY1);
            if(tr>0 && abs(tr)<Dir_len){
                //cout<<obscured<<endl;
                obscured=true;
                return obscured;
                //break;
            }
        }
        //cout<<obscured<<endl;
        ///here write codes
        return obscured;
    }
    double calculateAX_2(Ray *r)
    {
        double ax_2 = A * pow(r->dir.x,2) + B * pow(r->dir.y,2) + C * pow(r->dir.z,2) + D * r->dir.x * r->dir.y + E * r->dir.x * r->dir.z + F * r->dir.y * r->dir.z;
        return ax_2;
    }
    double calculateBX(Ray *r)
    {
        double bx = 2 * A * r->start.x * r->dir.x + 2 * B * r->start.y * r->dir.y + 2 * C * r->start.z * r->dir.z + D * (r->start.x * r->dir.y + r->start.y * r->dir.x) + E * (r->start.x * r->dir.z + r->start.z * r->dir.x) + F * (r->start.y * r->dir.z + r->dir.y * r->start.z) + G * r->dir.x + H * r->dir.y + I * r->dir.z;
        return bx;
    }
    double calculateC(Ray *r)
    {
        double c = A * pow(r->start.x,2) + B * pow(r->start.y,2) + C * pow(r->start.z,2) + D * r->start.x * r->start.y + E * r->start.x * r->start.z + F * r->start.y * r->start.z + G * r->start.x + H * r->start.y + I * r->start.z + J;
        return c;

    }
    double calculateDET(double a,double b,double c)
    {
        return pow(b,2) - 4 * a * c;
    }
    double getActualT(struct point intersectionPoint1,struct point intersectionPoint2,double t1,double t2)
    {
        bool check1=true,check2=true;
        ///check length boundary
        if(length!=0)
        {
            ///check out of boundary or not for point1
            if(intersectionPoint1.x<reference_point.x || intersectionPoint1.x > reference_point.x+length){
                check1=false;
            }
            ///check out of boundary or not for point2
            if(intersectionPoint2.x<reference_point.x || intersectionPoint2.x > reference_point.x+length){
                check2=false;
            }
        }
        ///check width boundary
        if(width!=0)
        {
            if(intersectionPoint1.y<reference_point.y || intersectionPoint1.y > reference_point.y+width){
                check1=false;
            }
            if(intersectionPoint2.y<reference_point.y || intersectionPoint2.y > reference_point.y+width){
                check2=false;
            }
        }
        ///check height boundary
        if(height!=0)
        {
            if(intersectionPoint1.z<reference_point.z || intersectionPoint1.z > reference_point.z+height){
                check1=false;
            }
            if(intersectionPoint2.z<reference_point.z || intersectionPoint2.z > reference_point.z+height){
                check2=false;
            }
        }
        if(check1 == true && check2 == true){
            if (t1 < t2) return t1;
            else return t2;

        }
        else if(check1 == true){
            return t1;
        }
        else if(check2 == true){
            return t2;
        }
        else return -1;

    }
    double get_tMIN(Ray *r)
    {
        ///ax^2 + bx + c



        double AX_2 = calculateAX_2(r);
        double BX = calculateBX(r);
        double C_ = calculateC(r);
        double DET = calculateDET(AX_2,BX,C_);
        if(AX_2 == 0){
            return -C_ / BX;
        }

        ///if det < -1 complex root
        if(DET < 0) return -1;

        ///else

        ///case 1:
        struct point intersectionPoint1;
        double t1 = (-BX - sqrt(DET))/(2 * AX_2);
        intersectionPoint1.x = r->start.x + r->dir.x * t1;
        intersectionPoint1.y = r->start.y + r->dir.y * t1;
        intersectionPoint1.z = r->start.z + r->dir.z * t1;

        ///case 2:
        struct point intersectionPoint2;
        double t2 = (-BX + sqrt(DET))/(2 * AX_2);
        intersectionPoint2.x = r->start.x + r->dir.x * t2;
        intersectionPoint2.y = r->start.y + r->dir.y * t2;
        intersectionPoint2.z = r->start.z + r->dir.z * t2;

        double t_actual = getActualT(intersectionPoint1,intersectionPoint2,t1,t2);
        return t_actual;


    }
    double intersect(Ray *r, double *current_color, int level){
        double tMIN=get_tMIN(r);
        if(tMIN<=0){
            return -1;
        }
        struct point intersectionPoint(r->start.x+(r->dir.x*tMIN),r->start.y+(r->dir.y*tMIN),r->start.z+(r->dir.z*tMIN));

        if(level==0){
            return tMIN;
        }
        ///setCurrentColor(current_color);

        current_color[0]=this->color[0]*this->coEfficients[0];
        current_color[1]=this->color[1]*this->coEfficients[0];
        current_color[2]=this->color[2]*this->coEfficients[0];
        ///-----------------------NORMAL CALC--------------------///

        ///calculate normal at intersect
        struct point normal_at_intersect,temp;
        ///       dF/dx
        temp.x = 2 * A * intersectionPoint.x + D * intersectionPoint.y + E * intersectionPoint.z + G;
        ///       dF/dy
        temp.y = 2 * B * intersectionPoint.y + D * intersectionPoint.x + F * intersectionPoint.z + H;
        ///       dF/dz
        temp.z = 2 * C * intersectionPoint.z + E * intersectionPoint.x + F * intersectionPoint.y + I;


        normal_at_intersect = unit_vector(temp);


        if(dot_prod(intersectionPoint, normal_at_intersect) > 0){
            normal_at_intersect.x = -normal_at_intersect.x;
            normal_at_intersect.y = -normal_at_intersect.y;
            normal_at_intersect.z = -normal_at_intersect.z;
        }
        ///------------------------END OF NORMAL CALC----------------------///

        for(int i=0;i<point_lights.size();i++){
            Ray RAY1;
            RAY1 = getRay(point_lights[i],intersectionPoint);
            double lengthOfLDir = getLenOfRAYDir(point_lights[i],intersectionPoint);


            bool obscured;


            if(isObscured(&RAY1,lengthOfLDir)==false){

                double lambertValue=dot_prod(RAY1.dir,normal_at_intersect);

                double phongValue=dot_prod(getReflectedRay(&RAY1,normal_at_intersect),r->dir);

                if(lambertValue <= 0) lambertValue = 0;
                if(phongValue <= 0) phongValue = 0;
                current_color[0]+=(point_lights[i]->color[0]*lambertValue*this->color[0]*this->coEfficients[1])+(point_lights[i]->color[0]*pow(phongValue,this->shine)*this->coEfficients[2]);
                current_color[1]+=(point_lights[i]->color[1]*lambertValue*this->color[1]*this->coEfficients[1])+(point_lights[i]->color[1]*pow(phongValue,this->shine)*this->coEfficients[2]);
                current_color[2]+=(point_lights[i]->color[2]*lambertValue*this->color[2]*this->coEfficients[1])+(point_lights[i]->color[2]*pow(phongValue,this->shine)*this->coEfficients[2]);

            }
        }
        ///code for spotlight
        for(int i=0;i<spot_lights.size();i++){
            PointLight pl = spot_lights[i]->point_light;

            ///cast rayl from pl.light_pos to intersectionPoint

            Ray RAY1,RAY2;
            RAY2 = getRay(&spot_lights[i]->point_light,intersectionPoint);
            RAY1 = getRay2(&spot_lights[i]->point_light,intersectionPoint);  ///get line from intersect to spotlight
            double lengthOfLDir = getLenOfRAYDir(&spot_lights[i]->point_light,intersectionPoint);
            struct point LineFactor = unit_vector(RAY1.dir);   ///get unitvector of line
            double spotFactor = dot_prod(LineFactor,spot_lights[i]->light_dir);  ///dotprod with line dir to get cos(theta)
            double OBJ_angle = (acos(spotFactor)*180)/pi;  ///get theta in degree




            ///check if obscured
            bool obscured;

            ///IF NOT OBSCURED
            if(isObscured(&RAY2,lengthOfLDir)==false && OBJ_angle < spot_lights[i]->cutoff_angle){   ///compare with cutoff_angle
                double lambertValue=dot_prod(RAY2.dir,normal_at_intersect);

                double phongValue=dot_prod(getReflectedRay(&RAY2,normal_at_intersect),r->dir);
                if(lambertValue <= 0) lambertValue = 0;
                if(phongValue <= 0) phongValue = 0;


                current_color[0]+=(spot_lights[i]->point_light.color[0]*lambertValue*this->color[0]*this->coEfficients[1])+(spot_lights[i]->point_light.color[0]*pow(phongValue,this->shine)*this->coEfficients[2]);//*this->color[0]);
                current_color[1]+=(spot_lights[i]->point_light.color[1]*lambertValue*this->color[1]*this->coEfficients[1])+(spot_lights[i]->point_light.color[1]*pow(phongValue,this->shine)*this->coEfficients[2]);//*this->color[1]);
                current_color[2]+=(spot_lights[i]->point_light.color[2]*lambertValue*this->color[2]*this->coEfficients[1])+(spot_lights[i]->point_light.color[2]*pow(phongValue,this->shine)*this->coEfficients[2]);//*this->color[2]);
            }
        }

       ///reflection
       if(level >= recursion_level) return tMIN;

        if(level < recursion_level){
            Ray ReflectedRay;
            ReflectedRay = getRecursiveRay(r,normal_at_intersect,intersectionPoint);

            int nearest;
            double tMin;


            nearest = calculateNearest(&ReflectedRay);
            tMin = calculatetMIN(&ReflectedRay);
            double *reflectedColor=new double[3];

            reflectedColor[0]=reflectedColor[1]=reflectedColor[2]=0.0;

            if(nearest!=-1){
                 tMin=objects[nearest]->intersect(&ReflectedRay,reflectedColor,level+1);
                 if(tMin!=-1){
                    current_color[0]+=reflectedColor[0]*this->coEfficients[3];
                    current_color[1]+=reflectedColor[1]*this->coEfficients[3];
                    current_color[2]+=reflectedColor[2]*this->coEfficients[3];
                }

            }

            clipColor(current_color[0],current_color[1],current_color[2]);

            delete[] reflectedColor;
            return tMin;
        }




    }
};
///--------------------generic class------------------------------///
