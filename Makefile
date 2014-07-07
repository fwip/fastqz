default: fastqz fapack fapacks

fastqz: fastqz15.cpp libzpaq.3.pod libzpaq.cpp libzpaq.h
	g++ -O2 -msse2 -lpthread fastqz15.cpp libzpaq.cpp -o $@

clean:
	- rm -f fastqz fapack fapacks
