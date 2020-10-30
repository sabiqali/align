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
"      --spacer_inserted                flag used to indicate that a prefix/suffix has been inserted along with a spacer sequence\n"
"  -s, --spacer=FILE                    the Spacer sequence in fastq/fasta FILE\n"
"  -r, --reads=FILE                     the ONT reads are in fastq FILE\n"
"  -c, --control=FILE                   the control oligos that are in a text FILE\n"
"  -f, --output_file=FILE               the file to which the score will be output to\n"
"      --table_out                      the output to file will be in a tabular format including only the score and names of the read and the control\n"
"      --align_out                      the output to file will be in an alignment format which will include the read, control oligo and the SAM cigar\n"
"      --SAM_out                        the output to file will be a SAM file\n"
"  -t, --threads=NUM                    use NUM threads (default: 1)\n"
"      --parasail                       use the parasail library to generate the alignment sequence\n"
"      --edlib                          use the edlib library to produce the alignment sequence\n\n";

namespace opt {
    static unsigned int print_alignment = 0;
    static std::string reads_file;
    static std::string control_file;
    static std::string output_file;
    static std::string spacer_file;
    static int table_out = 0;
    static int align_out = 0;
    static int parasail = 0;
    static int edlib = 0;
    static int num_threads = 1;
    static int sam_out = 0;
    static int spacer_inserted = 0;
}

struct AlignmentResult {
    int score;
    std::string probe_name;
    std::string query;
    std::string ref;
    std::string comp;
    std::string oligo;
    char orientation;
    char* cigar;
    unsigned char* alignment;
    int* endLocation;
    int alignmentLength;
};

static const char* shortopts = "r:c:f:t:p:s";

enum { OPT_HELP = 1, OPT_VERSION, OPT_TABLE_OUT, OPT_ALIGN_OUT, OPT_PARASAIL, OPT_EDLIB, OPT_SAM_OUT, OPT_SPACER_INSERTED };

static const struct option longopts[] = {
    { "print_alignment",      no_argument,       NULL, 'p' },
    { "reads",                required_argument, NULL, 'r' },
    { "control",              required_argument, NULL, 'c' },
    { "output_file",          required_argument, NULL, 'f' },
    { "threads",              required_argument, NULL, 't' },
    { "spacer",               required_argument, NULL, 's' },
    { "spacer_inserted",      no_argument,       NULL, OPT_SPACER_INSERTED },
    { "table_out",            no_argument,       NULL, OPT_TABLE_OUT },
    { "align_out",            no_argument,       NULL, OPT_ALIGN_OUT },
    { "sam_out",              no_argument,       NULL, OPT_SAM_OUT },
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
            case 's': arg >> opt::spacer_file; break;
            case OPT_SPACER_INSERTED: opt::spacer_inserted = 1; break;
            case OPT_TABLE_OUT: opt::table_out = 1; break;
            case OPT_ALIGN_OUT: opt::align_out = 1; break;
            case OPT_SAM_OUT: opt::sam_out = 1; break;
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

    if(opt::spacer_file.empty() && opt::spacer_inserted) {
        std::cerr << "Align: a --spacer file must be provided\n";
        die = true;
    }

    if(opt::output_file.empty() && (opt::align_out || opt::sam_out || opt::table_out)) {
        std::cerr << "Align: an --output file must be provided\n";
        die = true;
    }

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

//omp_set_num_threads(opt::num_threads);

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

std::string cigar_parser(char* cigar_to_parse) {
    std::string cigar(cigar_to_parse);
    std::string new_cigar = "";
    for (char& elem : cigar) {
	if( elem == '=' || elem == 'X')
	    elem = 'M';
	new_cigar += elem;
    }
    new_cigar[new_cigar.find('I')] = 'S';
    new_cigar[new_cigar.length() - 1] = new_cigar[new_cigar.length() - 1] == 'D'?'D':'S';
    std::cout<<new_cigar<<"\n";
    std::string final_cigar = "";
    int num_mem1 = 0, num_mem2 = 0;
    char char_mem1, char_mem2;

    for(int i = 0; i < new_cigar.length(); i++)  {
	if(isdigit(new_cigar[i])) {
	    num_mem1 = (num_mem1*10) + (new_cigar[i] - '0');
	    //std::cout<<num_mem1; 
	}
	if(new_cigar[i] == 'S') {
	    if(num_mem2 != 0)
	    final_cigar = final_cigar + std::to_string(num_mem2) + char_mem2;			
	    num_mem2 = num_mem1;
	    char_mem2 = new_cigar[i];
	    num_mem1 = 0;
	}
	if(new_cigar[i] == 'M' || new_cigar[i] == 'I' || new_cigar[i] == 'D') {
	    //std::cout<<num_mem1;
	    if(char_mem2 == new_cigar[i]) {
	        num_mem2 += num_mem1;
	    }
	    if(char_mem2 != new_cigar[i]) {
		final_cigar = final_cigar + std::to_string(num_mem2) + char_mem2;
		num_mem2 = num_mem1;
		char_mem2 = new_cigar[i];
	    }
	    //final_cigar = final_cigar + std::to_string(num_mem1) + "M";
	    num_mem1 = 0;
	}	
    }
    final_cigar = final_cigar + std::to_string(num_mem2) + char_mem2;

    return final_cigar;
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

void find_pos(std::string seq, std::string input, std::string spacer) {
    std::string temp = input;

    int cut_seq = 0;

    while(temp.find(spacer) != std::string::npos) {
        std::cout<< seq << "\t" << (temp.find(spacer) + cut_seq)<<"\t"<<spacer.length()<<"\t+" << std::endl;

        cut_seq = temp.find(spacer) + spacer.length();

        temp = temp.substr(cut_seq);
    }

    temp = dna_reverse_complement(input);

    while(temp.find(spacer) != std::string::npos) {
        std::cout<< seq << "\t" << (temp.find(spacer) + cut_seq)<<"\t"<<spacer.length()<<"\t-" << std::endl;

        cut_seq = temp.find(spacer) + spacer.length();

        temp = temp.substr(cut_seq);
    }

    return;
}

float percentage_identity(std::string comp) {
    std::size_t matches = std::count(comp.begin(),comp.end(),'|');

    return (((float)(matches*100)/(float)comp.length()));
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

inline AlignmentResult align_parasail(std::string reads, std::string control, parasail_matrix_t *matrix) {
    AlignmentResult complete_result;
    std::string reverse = dna_reverse_complement(reads);
    
    parasail_result_t* result = parasail_sg_trace_scan_16(reads.c_str(),reads.length(),control.c_str(),control.length(),5,4,matrix);
    parasail_traceback_t* traceback = parasail_result_get_traceback(result,reads.c_str(), reads.length(), control.c_str(), control.length(),matrix,'|','*','*');
    parasail_result_t* result_reverse = parasail_sg_trace_scan_16(reverse.c_str(),reverse.length(), control.c_str(), control.length(),5,4,matrix);
    parasail_traceback_t* traceback_reverse = parasail_result_get_traceback(result_reverse,reverse.c_str(), reverse.length(), control.c_str(), control.length(),matrix,'|','*','*');
    parasail_cigar_t* cigar = result->score > result_reverse->score ? parasail_result_get_cigar(result, reads.c_str(), reads.length(), control.c_str(), control.length(), matrix) : parasail_result_get_cigar(result_reverse, reverse.c_str(), reverse.length(), control.c_str(), control.length(), matrix);
    char* cigar_decoded = parasail_cigar_decode(cigar);
    
    if(result->score > result_reverse->score) {
	
        complete_result.score = result->score;
        complete_result.ref = traceback ->ref;
        complete_result.comp = traceback->comp;
        complete_result.query = traceback->query;
        complete_result.orientation = '+';
        complete_result.cigar = cigar_decoded;
    }
    else {
        complete_result.score = result_reverse->score;
        complete_result.ref = traceback_reverse ->ref;
        complete_result.comp = traceback_reverse->comp;
        complete_result.query = traceback_reverse->query;
        complete_result.orientation = '-';
        complete_result.cigar = cigar_decoded;
    }
    parasail_traceback_free(traceback);
    parasail_result_free(result);
    parasail_traceback_free(traceback_reverse);
    parasail_result_free(result_reverse);

    return complete_result;
}

inline AlignmentResult align_edlib(std::string reads, std::string control) {
    AlignmentResult complete_result;

    EdlibAlignResult result = edlibAlign(reads.c_str(), reads.length(), control.c_str(), control.length(), edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));
    
    complete_result.score = result.editDistance;
    complete_result.alignment = result.alignment;
    complete_result.alignmentLength = result.alignmentLength;
    complete_result.endLocation = result.endLocations;
        
    edlibFreeAlignResult(result);

    return complete_result;
}

inline AlignmentResult align(std::string reads, std::string control, parasail_matrix_t *matrix) {
    if(opt::parasail)
        return align_parasail(reads, control, matrix);
    if(opt::edlib)
        return align_edlib(reads, control);
}

int main(int argc, char *argv[])  {
    gzFile fp;
    kseq_t *seq;
    int l;
    std::ofstream out_fd;

    parse_align_options(argc , argv);

    fp = gzopen(opt::reads_file.c_str(), "r");
    seq = kseq_init(fp);
    std::ifstream control_fd(opt::control_file);

    if(opt::sam_out || opt::align_out || opt::table_out) {
        out_fd.open(opt::output_file);
    }
	
    std::unordered_map<std::string,std::string> control_substrings;
	
    parasail_matrix_t *matrix = parasail_matrix_create("ACGT", 5, -1);	

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
        std::cout << "read_name\toligo_name\tnum_alignments\tbest_score\tpercentage_identity\torientation" << std::endl;

    if(opt::table_out)
	out_fd << "read_name\toligo_name\tnum_alignments\tbest_sccore\tpercentage_identity\torientation" << std::endl;

    AlignmentResult best;
    AlignmentResult second_best;

    std::string spacer_sequence = "CGTATTGCGTCGAACAACCTCAACCCGTACGGAAACGGCTGAGCAGCAACTGGCACTTCAACTCCCCTGAAATAATACGCGTAAATCTGTATAGAGTTTCCTGGGGTTTTGTTCAGGACGCTCATCCCCAATCCCCTCTTTGTGTCCGACGACTCAGTTGCTTCTTAGGACATACTGACGCACTCCCGCACGAATGAGCCCTCCATTACATAACAAAAGGCATTGATCCATCACTGGTGGTAGTCGTTCATTGGTTAGTAACAGGAATAGCACCTATATGCGCCAGCGCTTATGCACATAAATGTGCATGAAGCCGTAAAAGCGACACACATTACGCAGGCTCTATGATAGAGCACATCTCCATCCTATCCACCTGAGTGAGCCAACCATCGCCTGTCTGATGCCTTTAGCGCTGCGACTCCTATAGGCCCCCCGTTGGCAAGTAAACCCCTACTCGCCCACCTTGACAGAGCTCTTCGGTATAGTCAGTATGACCAGGAACCATTGGCGCGGGCTAAGGCCAAAAAGTAACGCTACTGGAGGATTCGTTCATGCTAATGATCGGTATCCTGATTTGCTGTCGGTTGATATATAGGTTCCTGGACAGCGCTAGTCGATCCGTGAAGTGGAACTACGTTCAAACGGGTCCTCATGTAAGAGCAACGCGTAGTGTTAATTAGGCAGCAGTACAGCCCGATCGCTCCGCAGTTACCCATTTGGGATACGAATCTGGCCTTCGGATTTGCTGATCACGCAACCGGGGTGGTCGCCCATGAAACTCGCCGACAGATTAAGAGAGGCCTGCGACAGCAGCGTCGAGAAAGCATTTGAAGCCCTGGAGCTTCCATAGAGTATTATTGCTACCCCACATGGAACTTCATGCCCAACCGAACCCCTTAATGGGGAACCAATTCAAGCAAAGGGACTAGCCTTTACACTAGGCCGATCTATATCCTTCATCCTCAGTGACCCTGGTCCACACGTAGTGCTTAAAAGCATGGTTATGTCGATCAAGCATTCTTGACGGGTGAGGCAATTTCGGTACACCGAACCGTTCCCCTTGCGGTTGATCCGTATTGTCTCTTAGCAGGAGTTGTGCCAGCTGCTAACTCGCCCCGAAGGATTAACGGTATCAGGAGTGCAATCAAAGCAGCCTCTGTTAACACCCGGCTGGGGACGGGCTAATTCGCATACCTTATAATACCGAAGGGAACACAGCCGAGGGTTCTTATAAGCACTATACCCATTTACCGAACGGATATACTCGCGTGGTTTAGCTTCGGATGCCCGTTTTCACGAACAGATGTACAGATCCCTGACCTCATAGGACTTTATTCTGGTATGCCATCAACATGATGAAGCTCAGACAGTCGGACGTTTTGGACGTCTTCCATAATTGTGTTCTTCCGTAAGTCCAGCAAGGGAAGAAATAACCCACGTTTTAGGCATAAGACGTCTAGTCCATATCTCTAGCCGGCACGAATCAATCAGGATTACGCCCCGCCAGTGCCAGACCCAGGTGTGCGCCGCGGAGCCATAAGATGAAACCTCGCGTTTTCTACGTTATTTTCTCCTCGAAAGATTAAACGGACCCGGATTCATCCGTATCGGCTACCACGTTATACCGCCTTTTCGACTCAAGTGTGATTTTTGATTCCAATGAGCCCTACTGTGGCGTCTCACTGAGACAGAGTCAGTAAGTTTTATCTCACCGCGTCTCCCCCCAAAATTTACCACGTACTCCCCTGTTACATGTTCAATTCTCGTGCCATTTTTCGAAGGTAAGTTATTAGGCCCGTCGGCCGAGATAGGATGCCACAGTTAGAGGTGTCTTGCGTGGAACCTCTGGGAATAGTAAGACATCGAGACCCTCATCATGGGAGCATACGGGGTGTTGCCAACGTATCTTAGGAAGGCCCTTCTTATGGTGTAGTGCTTTCAAGGCTCAAACGACGCGTTTGTAGCTGCTTTTACTCCCGTTCAATTTCGGCCAATATAGA";

    int read_count = 0;
    while ((l = kseq_read(seq)) >= 0) {
	max = 0;
	second_max = 0;
        std::string sequence(seq->seq.s);
        std::unordered_map<std::string,int> test_substrings;

        std::string cut_string_temp = sequence;
        int cut_seq;

        if(opt::spacer_inserted) {
            while(cut_string_temp.find(spacer_sequence) != std::string::npos) {
                std::cout<< seq << "\t" << (cut_string_temp.find(spacer_sequence) + cut_seq)<<"\t"<<spacer_sequence.length()<<"\t+" << std::endl;

                cut_seq = cut_string_temp.find(spacer_sequence) + spacer_sequence.length();

                cut_string_temp = cut_string_temp.substr(cut_seq);

                if(cut_string_temp.find(spacer_sequence) != std::string::npos) {
                    std::string oligo = cut_string_temp.substr(0,cut_string_temp.find(spacer_sequence));

                    subString_test(sequence, k, test_substrings);
                }
                else {
                    std::string oligo = cut_string_temp;

                    subString_test(oligo, k, test_substrings);
                }
            }

            cut_string_temp = dna_reverse_complement(sequence);

            while(cut_string_temp.find(spacer_sequence) != std::string::npos) {
                std::cout<< seq << "\t" << (cut_string_temp.find(spacer_sequence) + cut_seq)<<"\t"<<spacer_sequence.length()<<"\t-" << std::endl;

                cut_seq = cut_string_temp.find(spacer_sequence) + spacer_sequence.length();

                cut_string_temp = cut_string_temp.substr(cut_seq);

                if(cut_string_temp.find(spacer_sequence) != std::string::npos) {
                    std::string oligo = cut_string_temp.substr(0,cut_string_temp.find(spacer_sequence));

                    subString_test(sequence, k, test_substrings);
                }
                else {
                    std::string oligo = cut_string_temp;

                    subString_test(oligo, k, test_substrings);
                }
            }
        }
        else
            subString_test(sequence, k, test_substrings);

        // Find the set of control oligos that share a k-mer with the read
        std::vector<std::string> control_set;
        for(auto itr = control_substrings.begin(); itr != control_substrings.end(); ++itr) {
            if( find_substring(itr->first,test_substrings) ) {
                if (std::find(control_set.begin(), control_set.end(), itr->second) == control_set.end())
		            control_set.push_back(itr->second);
            } 
        }

	best.score = 0;

	int num_alignments = 0;
	int vector_size = control_set.size();
	//#pragma omp parallel for
    	//for(const auto& elem: control_set) {
        #pragma omp parallel for
        for(auto i = 0; i<vector_size ; i++) {
	    std::string sequence2,temp,probe_name_new;
            std::string itr_elem = control_set[i];
	    int index = itr_elem.find('A');
	    if(index != std::string::npos)    
                sequence2 = itr_elem.substr(index); //used to separate the sequence from the probe no and name
	    else
		sequence2 = itr_elem;
	    index = itr_elem.find('F');
	    if(index != std::string::npos)
	        temp = itr_elem.substr(index); //used to extract the probe name
	    else
		temp = itr_elem;
	    index = temp.find('A');
	    if(index != std::string::npos)
	        probe_name_new = temp.substr(0,index);
	    else
		probe_name_new = temp;

	    AlignmentResult result = align(sequence,sequence2,matrix);
	    result.oligo = itr_elem;
	    result.probe_name = probe_name_new;
            #pragma omp critical
            if(result.score > best.score) {
                //result.oligo = itr_elem;
                //result.probe_name = probe_name_new;
                second_best = best;
                best = result;
            }
    	    
	    num_alignments++;
        }   
        
        read_count++;
	
        std::string percentage_identity_comp;
        int ind1 = best.comp.find('*');
        int ind2 = best.comp.find('|');
        int ind3 = best.comp.find_last_of('*');
        int ind4 = best.comp.find_last_of('|');
	//std::cout<<ind1<<" "<<ind2<<" "<<ind3<<" "<<ind4;
	if(ind2 != -1) {
	if(ind1 >= 0 && ind1 <= best.comp.length()) {
            int first_char = ind1<ind2?ind1:ind2;
            percentage_identity_comp = best.comp.substr( first_char , ind3>ind4?(ind3 - first_char):(ind4 - first_char));
	}
	else {
	    int first_char = ind2;
	    percentage_identity_comp = best.comp.substr(first_char, (ind4 - first_char));
 	}
	}
	else 
	    percentage_identity_comp = ' ';

        if(opt::sam_out && opt::parasail) {
            out_fd << seq->name.s << "\t" << (best.orientation == '+' ? "4" : "16") << "\t*\t0\t255\t" << cigar_parser(best.cigar) << "\t*\t0\t0\t" << seq->seq.s << "\t*\n"; 
        }

	if(opt::table_out && opt::parasail)
	    out_fd << seq->name.s << "\t" << best.probe_name << "\t" << num_alignments << "\t" << best.score << "\t" << percentage_identity(percentage_identity_comp) << "\t" << best.orientation << std::endl;

        if(opt::print_alignment) {
            std::cout<<"\n\n\nSequence "<<read_count<<" : "<<seq->seq.s;
            std::cout<<"\nAligned to: "<<read_aligns<<" control oligos";
            std::cout<<"\nHighest Score: "<<best.score<<"\nControl Oligo: "<<best.oligo<<"\n\n";
            if(opt::parasail) 
                std::cout<<best.query<<"\n"<<best.comp<<"\n"<<best.ref;
            if(opt::edlib)
                printAlignment(seq->seq.s, best.oligo.c_str(), best.alignment, best.alignmentLength, *(best.endLocation), EDLIB_MODE_NW);
            std::cout<<"\n\nSecond Highest Score:"<<second_best.score<<"\nControl Oligo: "<<second_best.oligo<<"\n\n";
            if(opt::parasail)
                std::cout<<second_best.query<<"\n"<<second_best.comp<<"\n"<<second_best.ref<<"\n\n";
	    std::cout<<best.cigar<<std::endl;
	    //if(opt::edlib)
                //printAlignment(seq->seq.s, second_best_oligo.c_str(), second_best_alignment, second_best_alignmentLength, *(second_best_endLocations), EDLIB_MODE_NW);
        }
        else {
            if(opt::parasail)
                std::cout << std::left << seq->name.s << "\t" << best.probe_name << "\t" << num_alignments << "\t" << best.score << "\t" << percentage_identity(percentage_identity_comp) << "\t" << best.orientation << std::endl;	
            if(opt::edlib)
                std::cout << std::left << seq->name.s << "\t" << best.probe_name << "\t" << num_alignments << "\t" << best.score << std::endl;
        }
	/*output_file<<"\n\n\nSequence "<<c<<" : "<<seq->seq.s;
        output_file<<"\nHighest Score: "<<max<<"\nProbe Name:"<<best_probe_no;
        output_file<<best_query<<"\n"<<best_comp<<"\n"<<best_ref;
        output_file<<"\n\nSecond Highest Score:"<<second_max<<"\nProbe name:"<<second_best_probe_no<<std::endl;
        output_file<<second_best_query<<"\n"<<second_best_comp<<"\n"<<second_best_ref<<"\n\n";*/
    }

    parasail_matrix_free(matrix);
    //fclose(out_file);
    if(opt::sam_out)
        out_fd.close();
    control_fd.close();
    //fclose(control_fd);
    kseq_destroy(seq);
    gzclose(fp);
    return 0;
}
