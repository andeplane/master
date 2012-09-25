CXX = c++
# Flags for optimized execs
CXXFLAGS = 
# CXXFLAGS = 
# Flags for debugging
#CXXFLAGS = -ftemplate-depth-30 -g -DBZ_DEBUG

LIBS = -lm -larmadillo

EXEC = dsmc
SRC = dsmc.cpp Molecule.cpp System.cpp Sorter.cpp lib.cpp Cell.cpp StatisticsSampler.cpp
OBJ = dsmc.o Molecule.o System.o Sorter.o lib.o Cell.o StatisticsSampler.o

$(SCR): 
	$(CXX) $(CXXFLAGS) -c $(SRC)  

$(EXEC): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $(EXEC) $(OBJ) $(LIBS)

clean:
	rm *.o
