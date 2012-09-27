CXX = g++
# Flags for optimized execs
CXXFLAGS = 
# Flags for debugging
# CXXFLAGS = -ftemplate-depth-30 -g -DBZ_DEBUG
# OPTFLAG = -O3
LIBS = -lm -larmadillo -fopenmp

EXEC = dsmc
SRC = dsmc.cpp Molecule.cpp System.cpp Sorter.cpp lib.cpp Cell.cpp StatisticsSampler.cpp
OBJ = dsmc.o Molecule.o System.o Sorter.o lib.o Cell.o StatisticsSampler.o

$(SCR): 
	$(CXX) $(CXXFLAGS) -c $(SRC)  

$(EXEC): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $(EXEC) $(OBJ) $(LIBS)

clean:
	rm *.o
