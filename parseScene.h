#ifndef __parseScene_h__
#define __parseScene_h__
#include <vector>
#include "RayTracer.h"
using namespace std;

extern float s;
extern boost::unordered_map<int,Vertex*> verticess;
extern boost::unordered_map<int,Sphere*> spheress;
extern boost::unordered_map<int,Triangle*> triangless;
extern boost::unordered_map<int,Material*> materialss;
extern vector<PLight*> lights;
extern RGB amLight;
extern int maxReflections;
extern RGB background;
extern vector<Camera*> cameras;
bool parseSceneXML(const char* filename);
extern float eps;

#endif //__parseScene_h__
