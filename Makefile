all:kseq.h parasail_align_alternate.cpp
		g++ -g -O2 parasail_align_alternate.cpp ../edlib/edlib/src/edlib.cpp -o parasail_align_alternate -I ../edlib/edlib/include -lz -fopenmp