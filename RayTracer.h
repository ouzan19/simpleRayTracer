
#include <vector>
#include <string>
#include <iostream>
#include <math.h>
#include <limits>
#include "boost/unordered_map.hpp"
#define constant 4*M_PI
using namespace std;

class RGB {
	public:
	float r;
	float g;
	float b;
	RGB(){
		this->r = 0.0;
		this->g = 0.0;
		this->b = 0.0;
		}
	RGB(float r,float g, float b){
		this->r = r;
		this->g = g;
		this->b = b;
		}
	RGB operator+ (RGB rgb){
			this->r += rgb.r ; 
			this->g += rgb.g ; 
			this->b += rgb.b ; 
			return  *this;
		}
	RGB operator += (RGB rgb){
			this->r += rgb.r ; 
			this->g += rgb.g ; 
			this->b += rgb.b ; 
			return  *this;
		}
	};


class Vector{
	public:
	float x;
	float y;
	float z;
	Vector(){
		this->x = 0.0;
		this->y = 0.0;
		this->z = 0.0;
		}
	Vector(float x,float y, float z){
		this->x = x;
		this->y = y;
		this->z = z;
		}
	
	Vector operator+(const Vector p){
		
		Vector result;
		result.x = this->x + p.x ;
		result.y = this->y + p.y ;
		result.z = this->z + p.z ;
		return result;
		
		}
		
		Vector operator-(const Vector p){
		
		Vector result;
		result.x = this->x - p.x ;
		result.y = this->y - p.y ;
		result.z = this->z - p.z ;
		return result;
		
		}
		
		float operator*(const Vector p){     /*dot product*/
			
			return this->x*p.x + this->y*p.y + this->z*p.z ;
			
			}
		Vector operator%(const Vector p){     /*cross product*/
				float _x = this->y*p.z  - this->z*p.y; 
				float _y = this->z*p.x - this->x*p.z; 
				float _z = this->x*p.y  - this->y *p.x;
				 return Vector(_x,_y,_z);
			}
		Vector operator*(const float p){
			
			Vector result;
		result.x = this->x *p ;
		result.y = this->y *p ;
		result.z = this->z *p ;
		return result;
			
			}
		Vector operator/(const float p){
			
			Vector result;
		result.x = this->x /p ;
		result.y = this->y /p ;
		result.z = this->z /p ;
		return result;
			
			}
			
		Vector operator!(){
				float length = sqrt ( x*x + y*y + z*z );
				return (*this)/length ; 
			
			}
	};
	
	
class Material{
	public:
	Material(RGB ambient,RGB diffuse,RGB specular,RGB reflectance,float specExp){
		
		this->ambient = ambient;
		this->diffuse=diffuse;
		this->specular =specular;
		this->reflectance=reflectance;
		this->specExp=specExp;
		
		}
	RGB ambient;
	RGB diffuse;
	RGB specular;
	RGB reflectance;
	float specExp;
	};
	
	
class Vertex{
	public:
	Vertex(Vector p,int id){
		this->coor = p;
		this->id=id;
		}
	Vector coor;
	int id;
	
	};
	
	
class ImagePlane{
	public:
	ImagePlane(){
		this->l = 0;
		this->r = 0;
		this->b = 0;
		this->t = 0;
		this->d = 0;
		this->HR = 0;
		this->VR = 0;
		}
	ImagePlane(float l,float r,float b,float t,float d,int HR,int VR){
		this->l = l;
		this->r = r;
		this->b = b;
		this->t = t;
		this->d = d;
		this->HR = HR;
		this->VR = VR;
		}
	float l;
	float r;
	float b;
	float t;
	float d;
	int HR;
	int VR;
	};


class Camera{
	public:
	Camera(){}
	Camera(Vector position,Vector gaze,Vector up,ImagePlane imPlane,string outName){
		this->position = position;
		this->w = gaze;
		this->v = up;
		this->u= w%v;
		this->imPlane = imPlane;
		this->outName = outName	;
		pixelLengthH = (imPlane.r - imPlane.l)/ imPlane.HR;
		pixelLengthV = (imPlane.t - imPlane.b)/ imPlane.VR;
		pixelStartX = imPlane.l;
		pixelStartY = imPlane.b;
		this->focal = imPlane.d;
		}	
	Vector position;
	Vector w;
	Vector v;
	Vector u;
	ImagePlane imPlane;
	string outName;
	float pixelLengthH;
	float pixelLengthV;
	float pixelStartX;
	float pixelStartY;
	float focal;
	};


class PLight{
	public:
	PLight(Vector p, RGB intensity){
	   this->position = p;
	   this->intensity = intensity; 
		}
	Vector position;
	RGB intensity;	
	};




class Shape {
	public:
	int Mid;
	virtual Vector getNormal(Vector p){};
	};
	
class Triangle: public Shape{
	public:
	int Vid1;
	int Vid2;
	int Vid3;
	Vector p1;
	Vector p2;
	Vector p3;
	Vector normal;
	Triangle(int Mid,int Vid1,int Vid2,int Vid3);
	Vector getNormal(Vector p){return this->normal;}
	void initialize();

	};


class Sphere:public Shape{
	public:
	int cid;
	float r;
	Vector c;
	Sphere(int cid,float r,int Mid);
	Vector getNormal(Vector p){return (p-c)*2;} 
	void initialize();
	
	};
enum Type { TRIANGLE,SPHERE};


class hitObject{
	public: 
		Type type;
		float t;
		int Sid;
		Vector center;
		float radius;
		hitObject(Type type, float t,int Sid){
			this->type = type;
			this->t =t ;
			this->Sid = Sid;
			}
		Shape* getPtr();
		
	};

class Ray{
	public:
	Vector direction;
	Vector origin;
	int reflections;
	Ray(){
		this->direction = Vector(0,0,0);
		this->origin = Vector(0,0,0);
		reflections=0;
		
		}
	Ray(Vector direction, Vector origin){
		this->direction =direction;
		this->origin = origin;
		reflections=0;
		}
	float intersect(Sphere *s){
				
				Vector d= this->direction;
				Vector e= this->origin;
				Vector c= s->c;
				Vector ec= e-c;
				float R = s->r;
				float A,B,C;    // Ax^2 + Bx + C = 0
				//cout << d.x<< " "<< d.y<< " "<< d.z<< endl;
				//cout << c.x<< " "<< c.y<< " "<< c.z<< endl;
				//cout << R << endl;
				
				A= d*d;
				B=2*(d*(ec));
				C=ec*ec - R*R;
				//cout << A<< " "<< B<< " "<< C<< endl;
				//cout << "*************"<<endl;
 				float disc = B*B - 4*A*C;
				//cout << "dis: "<< disc << endl;
				if(disc < 0) return (-1.0);
				else {
					disc = sqrt(disc);
					float root1 = (-B-disc)/(2*A);
					float root2 = (-B+disc)/(2*A);
					float tt= (root1 < root2) ? (root1):(root2);
						//cout << tt << endl;
					return tt;
				}
			}
				
	float intersect(Triangle *t){
	
	//cout <<"****************************************"<<endl;
	float a= t->p1.x - t->p2.x;
	float b= t->p1.y - t->p2.y;
	float c= t->p1.z - t->p2.z;
	float d= t->p1.x - t->p3.x;
	float e= t->p1.y - t->p3.y;
	float f= t->p1.z - t->p3.z;
	float g= this->direction.x;
	float h= this->direction.y;
	float i= this->direction.z; 
	//cout << this->direction.x<<" "<< this->direction.y<<" "<< this->direction.z<<endl;
	float j=t->p1.x - this->origin.x ;
	float k=t->p1.y - this->origin.y ;
	float l=t->p1.z - this->origin.z ;
	float ei_hf = e * i - h * f;
	float gf_di = g * f - d * i;
	float dh_eg = d * h - e * g;
	float ak_jb = a * k - j * b;
	float jc_al = j * c - a * l;
	float bl_kc = b * l - k * c;

	float M = a * ei_hf + b * gf_di + c * dh_eg;

	
	float gamma = (i * ak_jb + h * jc_al + g * bl_kc) / M;
	
	if(gamma < 0 || gamma > 1)
		return -1;
	
	float beta = (j * ei_hf + k * gf_di + l * dh_eg) / M;
		
	if((beta < 0 || (beta > 1 - gamma)) ) return -1;
	//cout << beta<<endl;
	return -(f * ak_jb + e * jc_al + d * bl_kc) / M;
	
		
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	/*
	
	float ei_hf = e*i - h*f;
	float gf_di = g*f - d*i;
	float dh_eg = d*h - e*g;
	float ak_jb = a*k - j*b;
	float jc_al = j*c - a*l;
	float bl_kc = b*l - k*c;
	float M = -(g*(b*f-c*e)+h*(d*c-a*f)+i*(a*e-d*b));//a*ei_hf + b*gf_di + c*dh_eg ;
	//cout <<"M: "<<M<<endl;
	if(M==0) return -1;
	float beta = -(j*ei_hf + k*gf_di + l*dh_eg) / M ;
	//cout <<"beta: "<<beta<<endl;
	if(beta < 0) return (-1);
	float gama = -(a*(k*i-h*l) + j*(h*c-b*i)+g*(b*l-c*k))/M;//(i*ak_jb + h*jc_al + g*bl_kc) / M ;
	//cout <<"gama: "<<gama<<endl;
	if(gama < 0 || beta + gama > 1) return (-1.0);
	float tt = -(j*(b*f-c*e) + k*(d*c-a*f)+l*(a*e-d*b))/M;  //- (f*ak_jb + e*jc_al + d*bl_kc) / M;
	//cout <<"t: "<<tt<<endl;
	// cout <<"t: "<< tt << endl;
	return tt;*/
	 
	 
}

	hitObject hit();
	};



class ImagePixel{
	public:
		Ray ray;
		float u;
		float v;
		
		/*ImagePixel(){
			this->u =0;
			this->v=0;
			Vector vv(0,0,0);
			Vector dd(0,0,0);
			this->ray(vv,dd);
			}*/
		void generateRay(Camera *c,int i,int j){
			//cout << "generating ray"<<endl;
				u = c->pixelStartX + (i+0.5)* c->pixelLengthH;
				v = c->pixelStartY + (j+0.5)* c->pixelLengthV;
				this->ray.direction.x = (-c->focal*c->w.x - u*c->u.x + v*c->v.x);
				this->ray.direction.y = (-c->focal*c->w.y - u*c->u.y + v*c->v.y);
				this->ray.direction.z = (-c->focal*c->w.z - u*c->u.z + v*c->v.z);
				this->ray.direction = !(this->ray.direction);
				//cout <<"ray:: "<<ray.direction.x<<" "<<ray.direction.y<<" "<<ray.direction.z<<" "<<endl;
				this->ray.origin = c->position;	
				this->ray.reflections =0;
				
				//cout << c->w.x << " " << c->w.y << " " << c->w.z << endl;
				//cout << c->u.x << " " << c->u.y << " " << c->u.z << endl;
				//cout << c->v.x << " " << c->v.y << " " << c->v.z << endl;
				//cout << c->focal << endl;
				//cout << "u: "<< u << "v: "<< v << endl;
				//cout << this->ray.direction.x<< " "<< this->ray.direction.y<< " "<< this->ray.direction.z<< endl;
				//cout << this->ray.origin.x <<" "<< this->ray.origin.y <<" "<< this->ray.origin.z<< endl;
			}
		
		RGB shading();
		//RGB calculateIntensity();
		RGB assignIntensity(Camera *c,int i,int j);
	};





