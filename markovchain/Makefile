CXX = g++
CXXFLAGS = -std=c++0x -g -Wall

markovchain: markovchain.o
	$(CXX) $(CXXFLAGS) markovchain.o -o markovchain.exe

markovchain.o: markovchain.cpp
	$(CXX) $(CXXFLAGS) -c markovchain.cpp
