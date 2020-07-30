all:kseq.h align.cpp
		g++ -g -std=c++14 -O2 align.cpp /.mounts/labs/simpsonlab/sw/edlib/edlib/src/edlib.cpp -o align -I /.mounts/labs/simpsonlab/sw/edlib/edlib/include -lz -fopenmp
