Tool to align a genome sequence to a reference genome sequence using the Parasail library and Edlib Library.

C++ Implementation:

Dependencies: 
1. Edlib
2. Parasail
3. g++
4. OpenMP
5. C++14

Usage:
./align [OPTIONS] --reads reads.fastq --control control_oligos.txt --output_file output_file.txt --parasail

  -p, --print_alignment                display best alignment for every read
      --version                        display version
      --help                           display this help and exit
  -r, --reads=FILE                     the ONT reads are in fastq FILE
  -c, --control=FILE                   the control oligos that are in a text FILE
  -f, --output_file=FILE               the file to which the score will be output to
      --table_out                      the output to file will be in a tabular format including only the score and names of the read and the control
      --align_out                      the output to file will be in an alignment format which will include the read, control oligo and the SAM cigar
      --SAM_out                        the output to file will be a SAM file
  -t, --threads=NUM                    use NUM threads (default: 1)
      --parasail                       use the parasail library to generate the alignment sequence
      --edlib                          use the edlib library to produce the alignment sequence

Mandatory Fields:
1. --reads
2. --control
3. --parasail or --edlib

Compilation:
1. Clone/donwload the repository into your drive.
2. Change the working directory to the one you just cloned.
3. Make sure you have g++ and all dependencies installed.
4. Run "make"
5. Run the tool using the usage statement provided above.



Python Implementation:

Dependencies:
1. Parasail
2. pysam
3. python2.7

Compilation:
1. Run the command python "python align.py <path to reads files with file name> <path to control file with file name>

Coded in Python and C++. 


