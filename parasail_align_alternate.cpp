

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
//#include "edlib.h"

KSEQ_INIT(gzFile, gzread)

std::unordered_map<std::string,int> control_substrings, test_substrings;

void subString_control(std::string s, int n, int sub_length)  
{ 
    //std::unordered_map<std::string,int> temp_array;
    //#pragma omp parallel for
    for (int i = 0; i < n; i++)  
        if(i+sub_length < n) {
                control_substrings.insert({s.substr(i, sub_length), 1});
        }
	
    //control_substrings = temp_array;
} 

void subString_test(std::string s, int n, int sub_length)  
{
    for (int i = 0; i < n; i++)  
        if(i+sub_length < n) {
                test_substrings.insert({s.substr(i, sub_length), 1});
        }
} 

bool find_substring(std::string to_find) {
    std::unordered_map<std::string,int>::const_iterator got = control_substrings.find (to_find);

    if ( got == control_substrings.end() ) {
        return false;
    }
    else {
        return true;
    }
}

std::string removeDupWord(std::string str)
{
   std::string word = "";
   for (auto x : str)
   {
       if (x == ' ')
       {
           word = "";
       }
       else
       {
           word = word + x;
       }
   }

   return word;
}

int main(int argc, char *argv[])
{
	gzFile fp;
	kseq_t *seq;
	int l;
	if (argc == 1) {
		fprintf(stderr, "Usage: %s <in.seq>\n", argv[0]);
		return 1;
	}
	fp = gzopen(argv[1], "r");
	seq = kseq_init(fp);
	std::ifstream control_fd(argv[2]);
	
	std::string one_line;
        std::string sequence2;
	
	parasail_result_t *best_result = NULL;
	parasail_matrix_t *matrix = NULL;
	//parasail_traceback_t  *best_traceback = NULL;
	//parasail_cigar_t *cigar;
	std::string best_query,best_ref,best_comp;

        matrix = parasail_matrix_create("ACGT", 5, -1);	
	//EdlibAlignResult result;
	//FILE *out_file;
	//out_file = fopen("traceback_output.txt", "w");

	int k = 25;
	

	while ((l = kseq_read(seq)) >= 0) {
		
		std::string sequence(seq->seq.s);
		char *fastq_seq = seq->seq.s;
		//int c = 0;
		int max = 0;
		int second_max = 0;
		
		subString_test(sequence, sequence.length(), k);
	
		control_fd.clear();
		control_fd.seekg(0,std::ios::beg);
		getline(control_fd,one_line);
		while(control_fd) {
			getline(control_fd,one_line);
			sequence2 = removeDupWord(one_line);
			if(one_line.length() > 75)
				sequence2 = one_line.substr(one_line.find('A'));
			
			//std::cout<<sequence2<<std::endl;
			subString_control(sequence2, sequence2.length(), k);

			char seq_arr[sequence2.length() + 1];
			strcpy(seq_arr,sequence2.c_str());
			if(sequence2.length() == 0)
				break;
			std::unordered_map<std::string,int>::iterator itr = test_substrings.begin();

			while (itr != test_substrings.end()) {
				//printf("testing4");
				//std::cout<<itr->first<<std::endl;
				if(find_substring(itr->first) == true) {
					//std::cout<<itr->first<<" Found\n";
					//printf("testing2");
					//c++;
				//	std::cout<<c<<" "<<sequence2.length()<<std::endl;
					//std::cout<<seq->seq.s<<std::endl;
					//std::cout<<seq_arr<<std::endl;
					parasail_result_t *result = parasail_sg_trace_scan_sat(fastq_seq,sequence.length(),seq_arr,sequence2.length(),5,4,matrix);
					parasail_traceback_t *traceback = parasail_result_get_traceback(result,seq->seq.s, sequence.length(), seq_arr, sequence2.length(),matrix,'|','*','*');
					//parasail_traceback_generic_extra(fastq_seq, sequence.length(), seq_arr, sequence2.length(), "Query: ","Target: ", matrix, result, '|','*','*',60,7,1,8,out_file);
					//cigar = parasail_result_get_cigar(result, fastq_seq,sequence.length(), seq_arr, sequence2.length(), matrix);
					//char *cigar_decoded = parasail_cigar_decode(cigar);
					//std::cout<<traceback->query<<std::endl;
					if(result->score > max) {
						second_max = max;
						max = result->score;
						best_result = result;
						best_ref = traceback ->ref;
						best_comp = traceback->comp;
						best_query = traceback->query;
					} 

					parasail_result_free(result);
					parasail_traceback_free(traceback);
					//std::cout<< traceback->query<<"\n"<<traceback->comp<<"\n"<<traceback->ref<<std::endl;
					/*result = edlibAlign(seq->seq.s, sequence.length(), seq_arr, sequence2.length(), edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));
					if (result.status == EDLIB_STATUS_OK) {
						if(result.editDistance>max > 0)
							max = result.editDistance;
					}*/
					break;
				}
					//Call Parasail here!!
				//else {
					//std::cout<<itr->first<<" Not found\n";
					//printf("testing3");
				//}
				itr++;
			}

			control_substrings.erase(control_substrings.begin(),control_substrings.end());
			//std::cout<<one_line;
			//std::cout<<probe_no<<std::endl<<probe_name<<std::endl<<seq<<std::endl<<one_line;
			//std::cout<<"seq "<<seq<<"\nProbe_no "<<probe_no<<"\n Probe name "<<probe_name<<std::endl;
			//std::cout<<"Found matching current test_sequence "<<c<<std::endl;
			//break;
			
		}
		test_substrings.erase(test_substrings.begin(), test_substrings.end());
		std::cout<<"Highest Score: "<<max<<"\nSecond Highest Score:"<<second_max<<std::endl;
		std::cout<< best_query<<"\n"<<best_comp<<"\n"<<best_ref<<std::endl;


	}
	//parasail_traceback_free(traceback);
	//parasail_result_free(result);
	parasail_matrix_free(matrix);
	//fclose(out_file);
	control_fd.close();
	//fclose(control_fd);
	kseq_destroy(seq);
	gzclose(fp);
	return 0;
}

