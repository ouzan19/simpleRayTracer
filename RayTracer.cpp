
#include "parseScene.h"

boost::unordered_map<int,Vertex*> verticess;
boost::unordered_map<int,Sphere*> spheress;
boost::unordered_map<int,Triangle*> triangless;
boost::unordered_map<int,Material*> materialss;
vector<PLight*> lights;
RGB amLight; 
int maxReflections;
RGB background;
float s;
vector<Camera*> cameras;
float eps=0;
 bool isShadow(Ray ray){
	
		/*search in triangles*/
			boost::unordered_map<int,Triangle*>::iterator  end = triangless.end(); 
			for(boost::unordered_map<int,Triangle*>::iterator  it = triangless.begin(); it !=end; ++it){
				float temp = ray.intersect(it->second);				
				if( temp > 0.00065 && temp<1+0.005 ) return true;
			
				}
	 
	 
			/*search in spheres*/
			boost::unordered_map<int,Sphere*>::iterator  end2= spheress.end();
			for(boost::unordered_map<int,Sphere*>::iterator  it = spheress.begin(); it !=end2; ++it){
				float temp = ray.intersect(it->second);	
				if(temp > 0.00065 && temp<1+0.005 ) return true;
				
				}
	 
		return false;
	 
	 
	 
	 
	 }
 Triangle::Triangle(int Mid,int Vid1,int Vid2,int Vid3){
		this->Mid = Mid;
		this->Vid1 = Vid1;
		this->Vid2 = Vid2;
		this->Vid3 = Vid3;
		
		
		}
void Triangle::initialize(){
		this->p1 = verticess.at(Vid1)->coor;
		this->p2 = verticess.at(Vid2)->coor;
		this->p3 = verticess.at(Vid3)->coor;
		this->normal= (p2-p1)%(p3-p1);

}
Sphere::Sphere(int cid,float r,int Mid){
		
		this->cid = cid;
		this->r =r;
		this->Mid = Mid;
		
		}

void Sphere::initialize(){this->c = verticess.at(cid)->coor;}



Shape* hitObject::getPtr(){
			if(this->type == TRIANGLE)
				return triangless.at(this->Sid);
			else return spheress.at(this->Sid);
			
			}
hitObject Ray::hit(){
		
			hitObject o(TRIANGLE,numeric_limits<float>::max(),-1);

			
			/*search in triangles*/
			boost::unordered_map<int,Triangle*>::iterator  end = triangless.end(); 
			for(boost::unordered_map<int,Triangle*>::iterator  it = triangless.begin(); it !=end; ++it){
				float temp_t = this->intersect(it->second);
				if((temp_t >eps) && (temp_t < o.t)){
					o.t= temp_t;
					o.Sid = it->first;	
					}
				}
				/*search in spheres*/
			boost::unordered_map<int,Sphere*>::iterator  end2= spheress.end();
			for(boost::unordered_map<int,Sphere*>::iterator  it = spheress.begin(); it !=end2; ++it){
				float temp_t = this->intersect(it->second);
				if((temp_t >eps) && (temp_t < o.t)){
					o.t= temp_t;
					o.type = SPHERE;
					o.Sid = it->first;				
					}
				}
				
				return o;
		}
RGB ImagePixel::shading(){
				
				hitObject o = this->ray.hit();
				if(o.Sid == -1){	// if no hitting object, return background
					//cout << "**********"<<endl;
					return background;//RGB(0,0,0);
					}
					eps=0.0001;
					//return RGB(0,0,255);
				RGB intens;
				Shape* s = o.getPtr();
				int index = s->Mid;
				Material* m = materialss.at(index);
				Vector intersectPoint = this->ray.origin + this->ray.direction*o.t;
				float kd_r = m->diffuse.r;
				float kd_g = m->diffuse.g;
				float kd_b = m->diffuse.b;
				Vector n = !(s->getNormal(intersectPoint));
				Vector v = this->ray.direction*(-1);
				float ks_r = m->specular.r;
				float ks_g = m->specular.g;
				float ks_b = m->specular.b;
				float p = m->specExp;
				
		
				intens.r+= amLight.r* m->ambient.r;
				intens.g+= amLight.g* m->ambient.g;
				intens.b+= amLight.b* m->ambient.b;
				int  numOfLights = lights.size();
				
				for(unsigned int i=0;i<numOfLights;i++){
						
						Vector l= lights[i]->position - intersectPoint;
						Ray tempRay(l,intersectPoint+this->ray.direction*(0.0001));
						
						if(isShadow(tempRay)) continue;					
						float distance = l*l;
						l = !l;
						float n_l = n*l; 
						n_l = (n_l>0 ? n_l : 0);
						Vector h = !(v+l);
						float n_h = n*h;
						n_h = (n_h>0 ? n_h : 0);
						float sik = constant*distance;
						//cout << sik << endl;
						float I_r = (lights[i]->intensity.r)/(sik);
						float I_g = (lights[i]->intensity.g)/(sik);
						float I_b = (lights[i]->intensity.b)/(sik);
						intens.r += I_r*(kd_r*n_l + ks_r*pow(n_h,p));
						intens.g += I_g*(kd_g*n_l + ks_g*pow(n_h,p));
						intens.b += I_b*(kd_b*n_l + ks_b*pow(n_h,p));
					}
				
						if((m->reflectance.r >0 || m->reflectance.g >0 || m->reflectance.b >0) && ray.reflections < maxReflections){
							this->ray.direction = !( n*((v*n)*2) - v );
							this->ray.origin = intersectPoint;//+this->ray.direction*(0.001);
							this->ray.reflections++;
							RGB temp = shading();
							intens.r += temp.r * m->reflectance.r ; 
							intens.g += temp.g * m->reflectance.g ; 
							intens.b += temp.b * m->reflectance.b ;
						
						}
						
						return intens;
			}

RGB ImagePixel::assignIntensity(Camera *c,int i,int j){ 
	 eps=1;
	 //cout << "*************************"<<endl;
	this->generateRay(c,i,j);
	 //cout << "*************************"<<endl;
	//cout << "i: " << i<< "j: "<< j <<endl;
	//cout << "origin: "<< this-> ray.origin.x << " "<< this-> ray.origin.y <<" "<< this-> ray.origin.z <<endl; 
	//cout << "direksiyon: "<< this-> ray.direction.x << " "<< this-> ray.direction.y <<" "<< this-> ray.direction.z <<endl;
	return this->shading();
}


	
	
	
	
	
	
	

