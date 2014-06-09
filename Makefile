default: fastqz

fastqz: fastqz15.cpp libzpaq.3.pod libzpaq.cpp libzpaq.h
	g++ -O3 -msse2 -s -lpthread fastqz15.cpp libzpaq.cpp -o $@
