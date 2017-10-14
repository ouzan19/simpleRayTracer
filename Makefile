SRC = test.cpp \
      writePPM.cpp \
      parseScene.cpp \
	  tinyXML/tinystr.cpp \
	  tinyXML/tinyxml.cpp \
	  tinyXML/tinyxmlerror.cpp \
	  tinyXML/tinyxmlparser.cpp \
	  RayTracer.cpp

OBJ = $(SRC:.cpp=.o)
CFLAGS = -O2
LIBS = 

.cpp.o:
	g++ -c $< $(CFLAGS)

all: raytracer

raytracer: $(OBJ) 
	g++ $(OBJ) -o raytracer $(CFLAGS) $(LIBS)

clean:
	rm -f *.o raytracer raytracer.ppm
