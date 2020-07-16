#! /usr/env/python

import parasail
import sys
import pysam

reads_file = sys.argv[1]
control_file = sys.argv[2]

# read in the control sequences
control_seqs = list()
control_fh = open(control_file)

# skip header
header = control_fh.readline()
for line in control_fh:
    control_seqs.append(line.rstrip().split()[1:3])
print(control_seqs[0])

scoring_matrix = parasail.matrix_create("ACGT", 5, -1)

align_data_file = open('align_data_file_complete.txt','w')

for read in pysam.FastxFile(reads_file):
    print(read)

    # align this read against all oligos
    best_score = 0
    idx = 0
    best_seq = []
    second_best_seq = []
    second_best = 0
    for (control_name, control_sequence) in control_seqs:
        #result = parasail.sg_trace_striped_16(read.sequence, control_sequence, 5, 4, scoring_matrix)
        result = parasail.sg_trace_scan_16(read.sequence, control_sequence, 5, 4, scoring_matrix)
        #result = parasail.sw_trace(read.sequence, control_sequence, 10, 1, parasail.blosum62)
        if idx % 100000 == 0:
            print("Aligned to %d controls" % (idx))
        idx += 1
        if result.score > best_score:
            traceback = result.traceback
            #print("new best score %d to %s" % (result.score, control_name))
            #print(traceback.query)
            #print(traceback.comp)
            #print(traceback.ref)
            #print(result.cigar.decode)
            #print()
            second_best = best_score
            second_best_seq = best_seq
            best_seq = [control_name, traceback.query, traceback.comp, traceback.ref,result.cigar.decode]
            best_score = result.score
    #print("Second best score %d to %s" % (second_best, second_best_seq[0]))
    align_data_file.write("Second best score %d to %s" % (second_best, second_best_seq[0]))
    #print(second_best_seq[1])
    align_data_file.write("\n")
    align_data_file.write(second_best_seq[1])
    #print(second_best_seq[2])
    align_data_file.write("\n")
    align_data_file.write(second_best_seq[2])
    #print(second_best_seq[3])
    align_data_file.write("\n")
    align_data_file.write(second_best_seq[3])
    align_data_file.write("\n")
    #print(second_best_seq[4])
    #align_data_file.write(second_best_seq[4])
    #print("Best score %d to %s" % (best_score, best_seq[0]))
    align_data_file.write("Best score %d to %s" % (best_score, best_seq[0]))
    #print(best_seq[1])
    align_data_file.write("\n")
    align_data_file.write(best_seq[1])
    #print(best_seq[2])
    align_data_file.write("\n")
    align_data_file.write(best_seq[2])
    #print(best_seq[3])
    align_data_file.write("\n")
    align_data_file.write(best_seq[3])
    align_data_file.write("\n")
    #print(best_seq[4])
    #align_data_file.write(best_seq[4])

align_data_file.close()

