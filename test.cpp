#include <iostream>
#include "parseScene.h"
#include <cstdlib>
#include "writePPM.h"
#include <time.h>
void generateImage(Camera* c);
ImagePixel *pixel = new ImagePixel();

int main(int argc, char* argv[])
{
	clock_t start = clock();
    //
    // Test XML parsing
    //
    bool result = parseSceneXML(argv[1]);
	cout << "******************************"<<endl;
    if (result)
    {
        std::cout << "Scene file parse successfully" << std::endl; 
    }
    else
    {
        std::cout << "ERROR parsing the scene file" << std::endl; 
    }
    cout << "******************************"<<endl;
	//RGB currentIntensity = pixel->assignIntensity(cameras[0],682,218);cout << currentIntensity.r<<" ";cout << currentIntensity.g<<" ";cout << currentIntensity.b<<endl;
	
			boost::unordered_map<int,Triangle*>::iterator  end = triangless.end(); 
			for(boost::unordered_map<int,Triangle*>::iterator  it = triangless.begin(); it !=end; ++it)
				it->second->initialize();
	 
	 
			/*search in spheres*/
			boost::unordered_map<int,Sphere*>::iterator  end2= spheress.end();
			for(boost::unordered_map<int,Sphere*>::iterator  it = spheress.begin(); it !=end2; ++it)
				it->second->initialize();
				





	int size = cameras.size();

	for(int i=0;i<size;i++)
		generateImage(cameras[i]);

	
   
	cout << (clock() - start)/ ((double) CLOCKS_PER_SEC) << endl;
    return 0;
}


void generateImage(Camera* c){
	
	int jj = c->imPlane.VR;
	int ii = c->imPlane.HR;
	float* data = (float*)malloc(sizeof(float)*3*ii*jj);
	for(int j=0;j<jj;j++){
		for(int i=0;i<ii;i++){
			
		
				//cout << i << " " << j <<"       "<<endl;
				int start = 3*((jj-j-1)*ii + i);
				RGB currentIntensity = pixel->assignIntensity(c,i,j);
				//cout << currentIntensity.r<<" ";cout << currentIntensity.g<<" ";cout << currentIntensity.b<<endl;
				data[start] = currentIntensity.r;
				data[start+1] = currentIntensity.g;
				data[start+2] = currentIntensity.b;
			}
		}
	 writePPM(c->outName.c_str(),ii, jj, data);
	}
