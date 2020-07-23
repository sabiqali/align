all:kseq.h parasail_align_alternate.cpp
		g++ -g -O2 parasail_align_alternate.cpp -I /.mounts/labs/simpsonlab/sw/parasail/2.4.2/include -L /.mounts/labs/simpsonlab/sw/parasail/2.4.2/lib -lparasail -o parasail_align_alternate -lz -fopenmp
