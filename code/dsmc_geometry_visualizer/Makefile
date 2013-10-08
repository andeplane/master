VPATH	= .


# target might be WINDOWS, OS_X or LINUX

TARGET = OS_X

# Home directory of the acw program

# point to the inclue directory
CINCLUDE = 
SOURCEDIR = src

INCLUDES = -I$(SOURCEDIR)
# library dir
LIBDIR = 

# compiler specific flags
CFLAGS =  -O3 -D$(TARGET)
# CFLAGS = -D$(TARGET)
FFLAGS = -lglew -lGLFW -framework OpenGL

PROJECT = main

_obj 	=  main.o Camera.o DSMCOpenGL.o CUtil.o CVector.o CMath.o mesh.o

obj = $(patsubst %,$(SOURCEDIR)/%,$(_obj))

CC 	= icpc

default: $(PROJECT)

$(PROJECT):  $(obj) 
	$(CC)  $(INCLUDES) -o $(PROJECT) $(obj) $(FFLAGS)  

%.o: %.cpp
	$(CC) -c -o $@ $^ $(INCLUDES) $(CFLAGS)   

%.o: %.cu
	$(CC) -c -o $@ $^ $(INCLUDES) $(CFLAGS)   


clean:	
	rm src/*.o


