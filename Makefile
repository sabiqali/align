all:kseq.h align.cpp
		gcc -g -std=c++14 -O2 align.cpp /.mounts/labs/simpsonlab/sw/edlib/edlib/src/edlib.cpp -I/.mounts/labs/simpsonlab/sw/WFA -L/.mounts/labs/simpsonlab/sw/WFA/build -I /.mounts/labs/simpsonlab/sw/edlib/edlib/include -I /.mounts/labs/simpsonlab/sw/parasail/2.4.2/include -L /.mounts/labs/simpsonlab/sw/parasail/2.4.2/lib -lstdc++ -lparasail -lwfa -o align -lz -fopenmp -ftree-vectorize
