all:kseq.h parasail_align_alternate.cpp
		g++ -g -std=c++14 -O2 parasail_align_alternate.cpp /.mounts/labs/simpsonlab/sw/edlib/edlib/src/edlib.cpp -o parasail_align_alternate -I /.mounts/labs/simpsonlab/sw/edlib/edlib/include -lz -fopenmp
