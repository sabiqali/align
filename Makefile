all:kseq.h align.cpp
		g++ -g -std=c++14  -O2 align.cpp -I /.mounts/labs/simpsonlab/sw/parasail/2.4.2/include -L /.mounts/labs/simpsonlab/sw/parasail/2.4.2/lib -lparasail -o align -lz -fopenmp
