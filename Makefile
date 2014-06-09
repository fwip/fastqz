default: fastqz fapack fapacks

fastqz: fastqz15.cpp libzpaq.3.pod libzpaq.cpp libzpaq.h
	g++ -O3 -msse2 -s -lpthread fastqz15.cpp libzpaq.cpp -o $@

clean:
	- rm -f fastqz fapack fapacks
