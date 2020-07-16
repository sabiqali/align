all:kseq.h parasail_align.cpp
		g++ -g -O2 parasail_align.cpp -I /.mounts/labs/simpsonlab/sw/parasail/2.4.2/include -L /.mounts/labs/simpsonlab/sw/parasail/2.4.2/lib -lparasail -o parasail_align -lz -fopenmp
