#include <iostream>
#include <cstdio>
#include "tinyXML/tinyxml.h"
#include "parseScene.h"
bool parseSceneXML(const char* filename)
{
	TiXmlDocument doc(filename);
	bool loadOkay = doc.LoadFile();

	if (!loadOkay)
	{
		std::cout << "Could not load file: " << filename << "Error = " << doc.ErrorDesc() << std::endl;
		return false;
	}
	
	TiXmlNode* pRoot = doc.FirstChild("Scene");
	for (TiXmlNode* pNode = pRoot->FirstChild(); pNode; pNode = pNode->NextSibling())
	{	
        if (pNode->Value() == std::string("Material"))
        {
			TiXmlAttribute* pAtt = pNode->ToElement()->FirstAttribute();
			int index = pAtt->IntValue(); // get material index

            //
            // read reflectance coefficients
            //
            float amb[3], dif[3], spe[3], mir[3];
			float phongExp;
			for (TiXmlNode* pChild = pNode->FirstChild(); pChild; pChild = pChild->NextSibling())
			{
				if (pChild->Value() == std::string("Ambient"))
				{
					sscanf(pChild->FirstChild()->Value(), "%f %f %f", &amb[0], &amb[1], &amb[2]);
				}
				else if (pChild->Value() == std::string("Diffuse"))
				{
					sscanf(pChild->FirstChild()->Value(), "%f %f %f", &dif[0], &dif[1], &dif[2]);
				}
				else if (pChild->Value() == std::string("Specular"))
				{
					sscanf(pChild->FirstChild()->Value(), "%f %f %f %f", &spe[0], &spe[1], &spe[2], &phongExp);
				}
				else if (pChild->Value() == std::string("Reflectance"))
				{
					sscanf(pChild->FirstChild()->Value(), "%f %f %f", &mir[0], &mir[1], &mir[2]);
				}
			}

            //
            // TODO: Save the scanned values for the current material in your data structures.
            //
            RGB ambi(amb[0],amb[1],amb[2]);
            RGB diff(dif[0],dif[1],dif[2]);
            RGB specular(spe[0],spe[1],spe[2]);
            RGB ref(mir[0],mir[1],mir[2]);
            Material* m= new Material(ambi,diff,specular,ref,phongExp);
            materialss[index] = m ;
            //cout << "material read" << endl;

        }
		else if (pNode->Value() == std::string("Vertex"))
        {
			TiXmlAttribute* pAtt = pNode->ToElement()->FirstAttribute();
			int index = pAtt->IntValue(); // get vertex index

			float coords[3];
			TiXmlNode* pChild = pNode->FirstChild();
			sscanf(pChild->FirstChild()->Value(), "%f %f %f", &coords[0], &coords[1], &coords[2]);

            //
            // TODO: Save the scanned values for the current vertex in your data structures.
            //
            Vector p(coords[0],coords[1],coords[2]);
            Vertex *v = new Vertex(p,index); 
            verticess[index] = v;
            //cout << "vertex read"<<endl;

		}
		else if (pNode->Value() == std::string("Triangle"))
        {
			//cout << "***************"<<endl;
			TiXmlAttribute* pAtt = pNode->ToElement()->FirstAttribute();
			int index = pAtt->IntValue(); // get triangle index
				//cout << "***************"<<endl;
			int vIndex[3], mIndex;
			for (TiXmlNode* pChild = pNode->FirstChild(); pChild; pChild = pChild->NextSibling())
			{	
				if (pChild->Value() == std::string("Vertices"))
				{
					sscanf(pChild->FirstChild()->Value(), "%d %d %d", &vIndex[0], &vIndex[1], &vIndex[2]);
				}
				else if (pChild->Value() == std::string("MaterialId"))
				{
					sscanf(pChild->FirstChild()->Value(), "%d", &mIndex);
				}
			}
	//cout << "***************"<<endl;
            //
            // TODO: Save the scanned values for the current triangle in your data structures.
            //
			Triangle *t = new Triangle(mIndex,vIndex[0],vIndex[1],vIndex[2]);
				//cout << "***************"<<endl;
			triangless[index] = t;
			//cout << "tri"<<endl;

		}
		else if (pNode->Value() == std::string("Sphere"))
        {
			TiXmlAttribute* pAtt = pNode->ToElement()->FirstAttribute();
			int index = pAtt->IntValue(); // Sphere index

			int vIndex, mIndex;
			float rad;
			for (TiXmlNode* pChild = pNode->FirstChild(); pChild; pChild = pChild->NextSibling())
			{
				if (pChild->Value() == std::string("Center"))
				{
					sscanf(pChild->FirstChild()->Value(), "%d", &vIndex);
				}
				else if (pChild->Value() == std::string("MaterialId"))
				{
					sscanf(pChild->FirstChild()->Value(), "%d", &mIndex);
				}
				else if (pChild->Value() == std::string("Radius"))
				{
					sscanf(pChild->FirstChild()->Value(), "%f", &rad);
				}
			}

            //
            // TODO: Save the scanned values for the current sphere in your data structures.
            //
            Sphere *s = new Sphere(vIndex,rad,mIndex);
            spheress[index] = s;
            //cout << "sphere"<<endl;

		}
		else if (pNode->Value() == std::string("PointLight"))
        {
			TiXmlAttribute* pAtt = pNode->ToElement()->FirstAttribute();
			int index = pAtt->IntValue(); // Light index

			float pos[3], intensity[3];
			for (TiXmlNode* pChild = pNode->FirstChild(); pChild; pChild = pChild->NextSibling())
			{
				if (pChild->Value() == std::string("Position"))
				{
					sscanf(pChild->FirstChild()->Value(), "%f %f %f", &pos[0], &pos[1], &pos[2]);
				}
				else if (pChild->Value() == std::string("Intensity"))
				{
					sscanf(pChild->FirstChild()->Value(), "%f %f %f", &intensity[0], &intensity[1], &intensity[2]);
				}
			}

            //
            // TODO: Save the scanned values for the current point light in your data structures.
            //
            Vector position(pos[0],pos[1],pos[2]);
            RGB intens(intensity[0],intensity[1],intensity[2]);
            PLight *pl = new PLight(position,intens);
            lights.push_back(pl);
            //cout << "plight read"<<endl;
		}
		else if (pNode->Value() == std::string("AmbientLight"))
        {
			float intensity[3];
			TiXmlNode* pChild = pNode->FirstChild();
			sscanf(pChild->Value(), "%f %f %f", &intensity[0], &intensity[1], &intensity[2]);

            //
            // TODO: Save the scanned values for the ambient light in your data structures.
            //
            amLight.r = intensity[0];
            amLight.g = intensity[1];
            amLight.b = intensity[2];
            //cout << "ambient read"<<endl;
		}
		else if (pNode->Value() == std::string("BackgroundColor"))
        {
            float bgColor[3];
			TiXmlNode* pChild = pNode->FirstChild();
			sscanf(pChild->Value(), "%f %f %f", &bgColor[0], &bgColor[1], &bgColor[2]);

            //
            // TODO: Save the scanned values for the background color in your data structures.
            //
            background.r = bgColor[0];
            background.g = bgColor[1];
            background.b = bgColor[2];
           // cout << "bg read"<<endl;
		}
		else if (pNode->Value() == std::string("RayReflectionCount"))
        {
            int rayReflectCount;
			TiXmlNode* pChild = pNode->FirstChild();
			sscanf(pChild->Value(), "%d", &rayReflectCount);

            //
            // TODO: Save the scanned values for the ray reflection count in your data structures.
            //
            maxReflections = rayReflectCount;
            //cout << "count read"<<endl;
		}
		else if (pNode->Value() == std::string("Camera"))
        {
			TiXmlAttribute* pAtt = pNode->ToElement()->FirstAttribute();
			int index = pAtt->IntValue(); // Camera index

			float gaze[3], up[3], pos[3];
			float left, right, bottom, top, distance;
			int nx, ny;
			std::string imageName;
			for (TiXmlNode* pChild = pNode->FirstChild(); pChild; pChild = pChild->NextSibling())
			{
				if (pChild->Value() == std::string("Position"))
				{
					sscanf(pChild->FirstChild()->Value(), "%f %f %f", &pos[0], &pos[1], &pos[2]);
				}
				else if (pChild->Value() == std::string("Gaze"))
				{
					sscanf(pChild->FirstChild()->Value(), "%f %f %f", &gaze[0], &gaze[1], &gaze[2]);
				}
				else if (pChild->Value() == std::string("Up"))
				{
					sscanf(pChild->FirstChild()->Value(), "%f %f %f", &up[0], &up[1], &up[2]);
				}
				else if (pChild->Value() == std::string("ImagePlane"))
				{
					sscanf(pChild->FirstChild()->Value(), "%f %f %f %f %f %d %d", &left, &right, &bottom, &top, &distance, &nx, &ny);
				}
				else if (pChild->Value() == std::string("OutputName"))
				{
					imageName = pChild->FirstChild()->Value();
				}
			}

            //
            // TODO: Save the scanned values for the current camera in your data structures.
            //
            
            Vector g(-gaze[0],-gaze[1],-gaze[2]);
            Vector u(up[0],up[1],up[2]);
            Vector po(pos[0],pos[1],pos[2]);
            ImagePlane imP(left, right, bottom, top, distance,nx,ny);
            Camera * cam = new Camera(po,g,u,imP,imageName);
            cameras.push_back(cam);
           // cout << "cam read"<<endl;
           
		}
	}

    return true;
}
