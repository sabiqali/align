#include <zlib.h>
#include <stdio.h>
#include "kseq.h"
#include "parasail.h"
#include <cstring>
#include <iostream>
#include <stdlib.h>
#include <stdbool.h>
#include <fstream>
#include <unordered_map>
#include <unordered_map>
#include "edlib.h"
#include <vector>
#include <iomanip>
#include <algorithm>
#include <unordered_set>
#include <chrono>
#include <sstream>
#include <getopt.h>
#include <omp.h>

KSEQ_INIT(gzFile, gzread)

//std::unordered_map<std::string,std::string> control_substrings;
//std::unordered_map<std::string,int> test_substrings;

static const char *ALIGN_MESSAGE = 
"Usage: ./align [OPTIONS] --reads reads.fastq --control control_oligos.txt --output_file output_file.txt --parasail\n"
"Align reads to control oligos provided.\n"
"\n"
"  -p, --print_alignment                display best alignment for every read\n"
"      --version                        display version\n"
"      --help                           display this help and exit\n"
"  -r, --reads=FILE                     the ONT reads are in fastq FILE\n"
"  -c, --control=FILE                   the control oligos that are in a text FILE\n"
"  -f, --output_file=FILE               the file to which the score will be output to\n"
"      --table_out                      the output to file will be in a tabular format including only the score and names of the read and the control\n"
"      --align_out                      the output to file will be in an alignment format which will include the read, control oligo and the SAM cigar\n"
"  -t, --threads=NUM                    use NUM threads (default: 1)\n"
"      --parasail                       watch the sequencing run directory DIR and call methylation as data is generated\n"
"      --edlib                          in watch mode, write the alignments for each fastq\n\n";

namespace opt {
    static unsigned int print_alignment = 0;
    static std::string reads_file;
    static std::string control_file;
    static std::string output_file;
    static int table_out = 0;
    static int align_out = 0;
    static int parasail = 0;
    static int edlib = 0;
    static int num_threads = 1;
}

static const char* shortopts = "r:c:f:t:p";

enum { OPT_HELP = 1, OPT_VERSION, OPT_TABLE_OUT, OPT_ALIGN_OUT, OPT_PARASAIL, OPT_EDLIB };

static const struct option longopts[] = {
    { "print_alignment",      no_argument,       NULL, 'p' },
    { "reads",                required_argument, NULL, 'r' },
    { "control",              required_argument, NULL, 'c' },
    { "output_file",          required_argument, NULL, 'f' },
    { "threads",              required_argument, NULL, 't' },
    { "table_out",            no_argument,       NULL, OPT_TABLE_OUT },
    { "align_out",            no_argument,       NULL, OPT_ALIGN_OUT },
    { "parasail",             no_argument,       NULL, OPT_PARASAIL },
    { "edlib",                no_argument,       NULL, OPT_EDLIB },
    { "help",                 no_argument,       NULL, OPT_HELP },
    { "version",              no_argument,       NULL, OPT_VERSION },
    { NULL, 0, NULL, 0 }
};

void parse_align_options(int argc, char** argv) {
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) {
            case 'p': opt::print_alignment = 1; break;
            case 'r': arg >> opt::reads_file; break;
            case 'c': arg >> opt::control_file; break;
            case 'f': arg >> opt::output_file; break;
            case '?': die = true; break;
            case 't': arg >> opt::num_threads; break;
            case OPT_TABLE_OUT: opt::table_out = 1; break;
            case OPT_ALIGN_OUT: opt::align_out = 1; break;
            case OPT_PARASAIL: opt::parasail = 1; break;
            case OPT_EDLIB: opt::edlib = 1; break;
            case OPT_HELP:
                std::cout << ALIGN_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cout << ALIGN_MESSAGE;
                exit(EXIT_SUCCESS);
        }
    }

    /*if(argc - optind > 0) {
        opt::region = argv[optind++];
    }

    if (argc - optind > 0) {
        std::cerr << SUBPROGRAM ": too many arguments\n";
        die = true;
    }*/

    if(opt::num_threads <= 0) {
        std::cerr << "Align : invalid number of threads: " << opt::num_threads << "\n";
        die = true;
    }

    if(opt::control_file.empty()) {
        std::cerr << "Align: a --control file must be provided\n";
        die = true;
    }

    if(opt::reads_file.empty()) {
        std::cerr << "Align: a --reads file must be provided\n";
        die = true;
    }

    /*if(opt::output_file.empty()) {
        std::cerr << "Align: a --output file must be provided\n";
        die = true;
    }*/

    if(opt::parasail == 0 && opt::edlib == 0) {
        std::cerr << "Align: processing option must be selected. Please selected either --parasail or --edlib\n";
        die = true;
    }

    if (die)
    {
        std::cout << "\n" << ALIGN_MESSAGE;
        exit(EXIT_FAILURE);
    }
}

inline void subString_control(const std::string& s, int sub_length,const std::string& probe_name, auto& control_substrings)  { 
    
    int n = s.length();
    //#pragma omp parallel for
    for (int i = 0; i < n; i++)
	//#pragma omp critical  
        if(i+sub_length < n) {
            control_substrings.insert({s.substr(i, sub_length), probe_name});
        }
	
    
} 

inline void subString_test(const std::string& s, int sub_length, auto& test_substrings)  {
    int n = s.length();
    
    for (int i = 0; i < n; i++)  
        if(i+sub_length < n) {
            test_substrings.insert({s.substr(i, sub_length), 1});
        }
} 

inline bool find_substring(const std::string& to_find, auto& test_substrings) {
    auto got = test_substrings.find (to_find);

    if ( got == test_substrings.end() ) {
        return false;
    } else  {
        return true;
    }
}

std::string removeDupWord(std::string str)  {
   std::string word = "";
   for (auto x : str)  {
       if (x == ' ')  {
           word = "";
       } else  {
           word = word + x;
       }
   }

   return word;
}

std::string dna_reverse_complement(std::string seq) {
    reverse(seq.begin(),seq.end());
    for (std::size_t i = 0; i < seq.length(); ++i){
        switch (seq[i]){
        case 'A':
            seq[i] = 'T';
            break;
        case 'C':
            seq[i] = 'G';
            break;
        case 'G':
            seq[i] = 'C';
            break;
        case 'T':
            seq[i] = 'A';
            break;
        }
    }
    return seq;
}

int percentage_identity(std::string comp) {
    std::size_t matches = std::count(comp.begin(),comp.end(),'|');

    return ((matches*100)/comp.length())
}

void printAlignment(const char* query, const char* target,
                    const unsigned char* alignment, const int alignmentLength,
                    const int position, const EdlibAlignMode modeCode) {
    int tIdx = -1;
    int qIdx = -1;
    if (modeCode == EDLIB_MODE_HW) {
        tIdx = position;
        for (int i = 0; i < alignmentLength; i++) {
            if (alignment[i] != EDLIB_EDOP_INSERT)
                tIdx--;
        }
    }
    for (int start = 0; start < alignmentLength; start += 50) {
        // target
        printf("T: ");
        int startTIdx = -1;
        for (int j = start; j < start + 50 && j < alignmentLength; j++) {
            if (alignment[j] == EDLIB_EDOP_INSERT)
                printf("-");
            else
                printf("%c", target[++tIdx]);
            if (j == start)
                startTIdx = tIdx;
        }
        printf(" (%d - %d)\n", std::max(startTIdx, 0), tIdx);

        // match / mismatch
        printf("   ");
        for (int j = start; j < start + 50 && j < alignmentLength; j++) {
            printf(alignment[j] == EDLIB_EDOP_MATCH ? "|" : " ");
        }
        printf("\n");

        // query
        printf("Q: ");
        int startQIdx = qIdx;
        for (int j = start; j < start + 50 && j < alignmentLength; j++) {
            if (alignment[j] == EDLIB_EDOP_DELETE)
                printf("-");
            else
                printf("%c", query[++qIdx]);
            if (j == start)
                startQIdx = qIdx;
        }
        printf(" (%d - %d)\n\n", std::max(startQIdx, 0), qIdx);
    }
}

int main(int argc, char *argv[])  {
    gzFile fp;
    kseq_t *seq;
    int l;

    parse_align_options(argc , argv);

    fp = gzopen(opt::reads_file.c_str(), "r");
    seq = kseq_init(fp);
    std::ifstream control_fd(opt::control_file);
	
    std::unordered_map<std::string,std::string> control_substrings;
	
    parasail_result_t *best_result = NULL;
    parasail_matrix_t *matrix = parasail_matrix_create("ACGT", 5, -1);	
    
    std::string best_query,best_ref,best_comp,best_oligo, best_probe_name;
    std::string second_best_query,second_best_ref,second_best_comp,second_best_oligo;

    int read_aligns;

    int k = 15;
    int max = 0,second_max = 0;	
    int control_count = 0;
    std::string probe_name;

    // discard header
    std::string one_line;
    getline(control_fd, one_line);

    // read the control sequences
    while(control_fd) {
        getline(control_fd,one_line);
        
	    std::string non_rep_section;
        if(one_line.length() > 75) {
            std::string sequence2 = one_line.substr(one_line.find('A')); //used to separate the sequence from the probe no and name
	        non_rep_section = sequence2.substr(61,68);   //indexed used to extract the non repeating sequence from the sequence obtained above
            std::string temp = one_line.substr(one_line.find('F')); //used to extract the probe name
            probe_name = temp.substr(0,temp.find('A'));
        }
        
        control_count++;
        subString_control(non_rep_section, k, one_line, control_substrings);
    }

    std::cerr << "Read " << control_count << " control sequences\n";
    if(!opt::print_alignment)
        std::cout << "read_name\toligo_name\tnum_alignments\tbest_sccore\tpercentage_identity\torientation" << std::endl;

    unsigned char  *best_alignment, *second_best_alignment;
    int best_alignmentLength, second_best_alignmentLength;
    int *best_endLocations, *second_best_endLocations;

    char orientation;

    int read_count = 0;
    while ((l = kseq_read(seq)) >= 0) {
	    max = 0;
	    second_max = 0;
        std::string sequence(seq->seq.s);
        std::unordered_map<std::string,int> test_substrings;
        subString_test(sequence, k, test_substrings);

        // Find the set of control oligos that share a k-mer with the read
        std::vector<std::string> control_set;
        for(auto itr = control_substrings.begin(); itr != control_substrings.end(); ++itr) {
            if( find_substring(itr->first,test_substrings) ) {
                if (std::find(control_set.begin(), control_set.end(), itr->second) == control_set.end())
		            control_set.push_back(itr->second);
            } 
        }

	    int num_alignments = 0;
	    int vector_size = control_set.size();
	    //#pragma omp parallel for
    	//for(const auto& elem: control_set) {
        #pragma omp parallel for
        for(auto i = 0; i<vector_size ; i++) {
            std::string itr_elem = control_set[i];    
            std::string sequence2 = itr_elem.substr(itr_elem.find('A')); //used to separate the sequence from the probe no and name
            std::string temp = itr_elem.substr(itr_elem.find('F')); //used to extract the probe name
            std::string probe_name_new = temp.substr(0,temp.find('A'));
            if(opt::parasail) {
                parasail_result_t* result = parasail_sg_trace_scan_16(seq->seq.s,sequence.length(),sequence2.c_str(),sequence2.length(),5,4,matrix);
                parasail_traceback_t* traceback = parasail_result_get_traceback(result,seq->seq.s, sequence.length(), sequence2.c_str(), sequence2.length(),matrix,'|','*','*');
		        #pragma omp critical
                if(result->score > max) {
                    second_max = max;
                    max = result->score;
                    best_probe_name = probe_name_new;
                    second_best_query = best_query;
                    second_best_ref = best_ref;
                    second_best_comp = best_comp;
                    best_ref = traceback ->ref;
                    best_comp = traceback->comp;
                    best_query = traceback->query;
                    second_best_oligo = best_oligo;
                    best_oligo = itr_elem;
                    orientation = '+'
                }
                parasail_traceback_free(traceback);
                parasail_result_free(result);
            }
            if(opt::edlib) {
                EdlibAlignResult result = edlibAlign(seq->seq.s, sequence.length(), sequence2.c_str(), sequence2.length(), edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));
	    	    #pragma omp critical
                if(result.editDistance > max) {
                    second_max = max;
                    max = result.editDistance;
                    best_probe_name = probe_name_new;
                    second_best_oligo = best_oligo;
                    best_oligo = sequence2;
                    second_best_alignment = best_alignment;
                    best_alignment = result.alignment;
                    second_best_alignmentLength = best_alignmentLength;
                    best_alignmentLength = result.alignmentLength;
                    second_best_endLocations = best_endLocations;
                    best_endLocations = result.endLocations;
                }
                edlibFreeAlignResult(result);
            }
    	    
	        num_alignments++;
        }   
        
        read_count++;

        if(opt::print_alignment) {
            std::cout<<"\n\n\nSequence "<<read_count<<" : "<<seq->seq.s;
            std::cout<<"\nAligned to: "<<read_aligns<<" control oligos";
            std::cout<<"\nHighest Score: "<<max<<"\nControl Oligo: "<<best_oligo<<"\n\n";
            if(opt::parasail)
                std::cout<<best_query<<"\n"<<best_comp<<"\n"<<best_ref;
            if(opt::edlib)
                printAlignment(seq->seq.s, best_oligo.c_str(), best_alignment, best_alignmentLength, *(best_endLocations), EDLIB_MODE_NW);
            std::cout<<"\n\nSecond Highest Score:"<<second_max<<"\nControl Oligo: "<<second_best_oligo<<"\n\n";
            if(opt::parasail)
                std::cout<<second_best_query<<"\n"<<second_best_comp<<"\n"<<second_best_ref<<"\n\n";
            //if(opt::edlib)
                //printAlignment(seq->seq.s, second_best_oligo.c_str(), second_best_alignment, second_best_alignmentLength, *(second_best_endLocations), EDLIB_MODE_NW);
        }
        else
            std::cout << std::left << seq->name.s << "\t" << best_probe_name << "\t" << num_alignments << "\t" << max << "\t" << percentage_identity(best_comp) << "\t" << orientation << std::endl;	

	    /*output_file<<"\n\n\nSequence "<<c<<" : "<<seq->seq.s;
        output_file<<"\nHighest Score: "<<max<<"\nProbe Name:"<<best_probe_no;
        output_file<<best_query<<"\n"<<best_comp<<"\n"<<best_ref;
        output_file<<"\n\nSecond Highest Score:"<<second_max<<"\nProbe name:"<<second_best_probe_no<<std::endl;
        output_file<<second_best_query<<"\n"<<second_best_comp<<"\n"<<second_best_ref<<"\n\n";*/
    }

    parasail_matrix_free(matrix);
    //fclose(out_file);
    control_fd.close();
    //fclose(control_fd);
    kseq_destroy(seq);
    gzclose(fp);
    return 0;
}

