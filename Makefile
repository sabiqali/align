all:kseq.h align.cpp
		g++ -g -std=c++14 -O2 align.cpp /.mounts/labs/simpsonlab/sw/edlib/edlib/src/edlib.cpp -I /.mounts/labs/simpsonlab/sw/edlib/edlib/include -I /.mounts/labs/simpsonlab/sw/parasail/2.4.2/include -L /.mounts/labs/simpsonlab/sw/parasail/2.4.2/lib -lparasail -o align -lz -fopenmp
