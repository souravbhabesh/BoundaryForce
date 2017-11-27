SRC = $(wildcard *.cpp)
OBJS = $(SRC:.cpp=.o)
DEPS = variables.h
CXX = g++
DEBUG = -g
CXXFLAGS = -Wall -c $(DEBUG) -std=c++11
LFLAGS = $(DEBUG) -O2 -Wall 

FORCE : $(OBJS)
	$(CXX) -o FORCE $(OBJS) $(LFLAGS)

gsd.o : $(DEPS) gsd.h  gsd_fn.cpp
	$(CXX) $(CXXFLAGS) gsd_fn.cpp

force.o : $(DEPS) force.h  force.cpp 
	$(CXX) $(CXXFLAGS) force.cpp

gsd_read.o : $(DEPS) gsd.h gsd_read.h  gsd_read.cpp
	$(CXX) $(CXXFLAGS) gsd_read.cpp

main.o : $(DEPS) gsd.h gsd_read.h force.h main.cpp
	$(CXX) $(CXXFLAGS) main.cpp

clean:
	\rm *.o *~ FORCE

