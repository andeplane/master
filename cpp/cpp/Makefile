TARGET = LINUX

# Home directory of the acw program

# point to the inclue directory
CINCLUDE =  
SOURCEDIR = .

INCLUDES = -I/sw/include -I$(CINCLUDE)  -I$(SOURCEDIR)
# library dir
LIBDIR =  CiRIOBase/libs

# compiler specific flags
CFLAGS =  -O3 -D$(TARGET)

FFLAGS = -lm

PROJECT = main

_obj 	= cell.o grid.o main.o molecule.o Random.o Sorter.o statisticssampler.o Stdafx.o system.o Wall.o


obj = $(patsubst %,$(SOURCEDIR)/%,$(_obj))

CC 	= g++ -m64

default: $(PROJECT)

$(PROJECT):  $(obj) 
	$(CC)  -L$(LIBDIR)  $(INCLUDES) -o $(PROJECT) $(obj) $(FFLAGS)  

%.o: %.cpp
	$(CC) -c -o $@ $^ $(INCLUDES) $(CFLAGS)   

%.o: %.cu
	$(CC) -c -o $@ $^ $(INCLUDES) $(CFLAGS)   


clean:	
	rm *.o source/*.o work	


