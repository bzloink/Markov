CXX = g++
CXXFLAGS = -std=c++0x -g -Wall

markovcluster: markovcluster.o
	$(CXX) $(CXXFLAGS) markovcluster.o -o markovcluster.exe

markovcluster.o: markovcluster.cpp
	$(CXX) $(CXXFLAGS) -c markovcluster.cpp
