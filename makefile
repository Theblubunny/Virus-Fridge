CXX = g++
CXXFLAGS = -Wall -g

proj4: vdetect.o mytest.o
	$(CXX) $(CXXFLAGS) vdetect.o mytest.o -o proj4

mytest.o: vdetect.o mytest.cpp
	$(CXX) $(CXXFLAGS) -c mytest.cpp

vdetect.o: vdetect.cpp vdetect.h
	$(CXX) $(CXXFLAGS) -c vdetect.cpp

clean:
	rm *.o*
	rm *~
