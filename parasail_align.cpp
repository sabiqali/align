
#include <zlib.h>
#include <stdio.h>
#include "kseq.h"
#include "parasail.h"
#include <cstring>
#include <iomanip>
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
    // Pick starting point in outer loop 
    // and lengths of different strings for 
    // a given starting point  
    std::unordered_map<std::string,int> temp_array;
    //#pragma omp parallel for
    for (int i = 0; i < n; i++)  
        if(i+sub_length < n) {
                //std::cout << s.substr(i, sub_length) << std::endl;
                control_substrings.insert({s.substr(i, sub_length), 1});
        }
	
    //control_substrings = temp_array;
} 

void subString_test(std::string s, int n, int sub_length)  
{ 
    // Pick starting point in outer loop 
    // and lengths of different strings for 
    // a given starting point 
    for (int i = 0; i < n; i++)  
        if(i+sub_length < n) {
                //std::cout << s.substr(i, sub_length) << std::endl;
                test_substrings.insert({s.substr(i, sub_length), 1});
        }
} 

bool find_substring(std::string to_find) {
    std::unordered_map<std::string,int>::const_iterator got = control_substrings.find (to_find);

    if ( got == control_substrings.end() ) {
        //cout << "not found";
        //cout << endl;
        return false;
    }
    else {
        //cout << got->first << " is " << got->second;
        //cout << endl;
        return true;
    }
}

int removeDupWord(std::string str)
{
   int index = 0;
   char seq_arr[str.length() + 1];
   strcpy(seq_arr,str.c_str());
   //std::cout<<str.find('A');
   for(int i = 0; i<str.length(); i++) {
	//std::cout<<seq_arr[i];
   	if(str[i] = 'A') {
		//std::cout<<"testing"<<str[i]<<i;
		index = i;
		break;
	}
   }
   return index;
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
	//std::FILE* control_fd(argv[2]);
	std::string one_line;
        std::string sequence2;
	std::ofstream output_file(argv[3]);
	parasail_result_t *best_result = NULL;
	parasail_matrix_t *matrix = NULL;
	//parasail_traceback_t *traceback,*best_traceback,*second_best_traceback;
	std::string best_query,best_ref,best_comp,second_best_query,second_best_ref,second_best_comp;

        matrix = parasail_matrix_create("ACGT", 5, -1);	
	//EdlibAlignResult result;
	output_file<<std::left<<std::setw(55)<<"sequence_name"<<std::setw(35)<<"best_oligo_name"<<std::setw(35)<<"best_score"<<std::setw(35)<<"second_best_oligo_name"<<std::setw(35)<<"second_best_score"<<std::endl;
	int k = 25;
	//std::cout<<k;
	int c = 1;

	while ((l = kseq_read(seq)) >= 0) {
		//printf("name: %s\n", seq->name.s);
		//if (seq->comment.l) printf("comment: %s\n", seq->comment.s);
		//printf("seq: %s\n", seq->seq.s);
		std::string sequence(seq->seq.s);
		int max = 0;
		int second_max = 0;
		//std::cout<<seq->seq.s<<std::endl;
		//std::ifstream control_fd(argv[2])
		//std::cout<<sequence<<std::endl;
		//strcpy(sequence,seq->seq.s);
		//subString_test(sequence, sequence.length(), k);
		//fseek(control_fd, 0, SEEK_SET);
		control_fd.clear();
		control_fd.seekg(0,std::ios::beg);
		getline(control_fd,one_line);
		std::string probe_no,probe_name,best_probe_no,best_probe_name,second_best_probe_no,second_best_probe_name;
		
		while(control_fd) {
			getline(control_fd,one_line);
			std::string temporary(one_line);
		
			if(one_line.length() > 75) {
                                sequence2 = temporary.substr(temporary.find('A'));
				std::string temp2 = temporary.substr(temporary.find('F'));
				probe_no = one_line.substr(0,one_line.find(' '));
				probe_name = temp2.substr(0,temp2.find('A'));
				//std::cout<<temp2<<std::endl;

			}
			//std::cout<<probe_name<<std::endl;
			//subString_control(sequence2, sequence2.length(), k);

			char seq_arr[sequence2.length() + 1];
			strcpy(seq_arr,sequence2.c_str());
			if(sequence2.length() == 0)
				break;
			//std::unordered_map<std::string,int>::iterator itr = test_substrings.begin();

			//while (itr != test_substrings.end()) {
				//printf("testing4");
				//if(find_substring(itr->first) == true) {
					//std::cout<<itr->first<<" Found\n";
					//printf("testing2");
					//c++;
					//std::cout<<seq->seq.s<<std::endl;
					//std::cout<<"testing1"<<seq_arr<<std::endl;
					parasail_result_t *result = parasail_sg_trace_scan_16(seq->seq.s,sequence.length(),seq_arr,sequence2.length(),5,4,matrix);
					//std::cout<<"testing2"<<seq_arr<<std::endl;
					parasail_traceback_t *traceback = parasail_result_get_traceback(result,seq->seq.s, sequence.length(), seq_arr, sequence2.length(),matrix,'|','*','*');
					//std::cout<<traceback->ref<<std::endl<<std::endl;
					if(result->score > max) {
						second_max = max;
						max = result->score;
						best_result = result;
						second_best_probe_no = best_probe_no;
						second_best_probe_name = best_probe_name;
						best_probe_no = probe_no;
						best_probe_name = probe_name;
						second_best_query = best_query;
						second_best_ref = best_ref;
						second_best_comp = best_comp;
						best_ref = traceback ->ref;
						best_comp = traceback->comp;
						best_query = traceback->query;
						//std::cout<<traceback->ref<<std::endl;
					} 
					//std::cout<<traceback->ref<<std::endl;
					parasail_traceback_free(traceback);
					parasail_result_free(result);
					/*result = edlibAlign(seq->seq.s, sequence.length(), seq_arr, sequence2.length(), edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));
					if (result.status == EDLIB_STATUS_OK) {
						if(result.editDistance>max > 0)
							max = result.editDistance;
					}*/
					//break;
				//}
					//Call Parasail here!!
				//else {
					//std::cout<<itr->first<<" Not found\n";
					//printf("testing3");
				//}
				//itr++;
			//}

			//control_substrings.erase(control_substrings.begin(),control_substrings.end());
			//std::cout<<one_line;
			//std::cout<<probe_no<<std::endl<<probe_name<<std::endl<<seq<<std::endl<<one_line;
			//std::cout<<"seq "<<seq<<"\nProbe_no "<<probe_no<<"\n Probe name "<<probe_name<<std::endl;
			//std::cout<<"Found matching current test_sequence "<<c<<std::endl;
			//break;
		}
		//subString_test(sequence, sequence.length(), 10);
		//test_substrings.erase(test_substrings.begin(),test_substrings.end());
		//std::cout<<"Testing"<<best_ref<<std::endl;
		//std::cout<<"testing2"<<second_best_ref<<std::endl;
		/*std::cout<<"\n\n\nSequence "<<c<<" : "<<seq->seq.s;
		std::cout<<"\nHighest Score: "<<max<<"\nProbe Name:"<<best_probe_no<<"\n\n";
		std::cout<<best_query<<"\n"<<best_comp<<"\n"<<best_ref;
		std::cout<<"\n\nSecond Highest Score:"<<second_max<<"\nProbe name:"<<second_best_probe_no<<"\n\n";
		std::cout<<second_best_query<<"\n"<<second_best_comp<<"\n"<<second_best_ref<<"\n\n";*/
		std::cout<<"Calculating "<<c<<"....."<<std::endl;
		//std::cout<<best_probe_name<<"\t"<<second_best_probe_name<<std::endl;
		
		output_file<<std::left<<std::setw(55)<<seq->name.s<<std::setw(35)<<best_probe_name<<std::setw(35)<<max<<std::setw(35)<<second_best_probe_name<<std::setw(35)<<second_max<<std::endl;
                /*output_file<<"\n\n\nSequence "<<c<<" : "<<seq->seq.s;
                output_file<<"\nHighest Score: "<<max<<"\nProbe Name:"<<best_probe_no;
		output_file<<best_query<<"\n"<<best_comp<<"\n"<<best_ref;
                output_file<<"\n\nSecond Highest Score:"<<second_max<<"\nProbe name:"<<second_best_probe_no<<std::endl;
		output_file<<second_best_query<<"\n"<<second_best_comp<<"\n"<<second_best_ref<<"\n\n";
		*/
		c++;
		//output_file<<"\nHighest Score: ";
		//output_file<<max;
		//output_file<<"\nSecond Highest Score: ";
		//output_file<<second_max;
		//break;
		//if (seq->qual.l) printf("qual: %s\n", seq->qual.s);
	}

	/*while ((l = kseq_read(seq)) >= 0) {
		//printf("name: %s\n", seq->name.s);
		//if (seq->comment.l) printf("comment: %s\n", seq->comment.s);
		//printf("seq: %s\n", seq->seq.s);
		std::string sequence(seq->seq.s);
		//std::cout<<sequence<<std::endl;
		//strcpy(sequence,seq->seq.s);
		subString_test(sequence, sequence.length(), 10);
		
		//subString_test(sequence, sequence.length(), 10);
		break;
		//if (seq->qual.l) printf("qual: %s\n", seq->qual.s);
	}

	while(control_fd) {
		getline(control_fd,one_line);
		sequence2 = removeDupWord(one_line);
		subString_control(sequence2, sequence2.length(), 10);
		//std::cout<<one_line;
		//std::cout<<probe_no<<std::endl<<probe_name<<std::endl<<seq<<std::endl<<one_line;
		//std::cout<<"seq "<<seq<<"\nProbe_no "<<probe_no<<"\n Probe name "<<probe_name<<std::endl;
		break;
	}

    control_fd.close();

	std::unordered_map<std::string,int>::iterator itr = test_substrings.begin();

	//printf("testing");

    while (itr != test_substrings.end()) {
		//printf("testing4");
        if(find_substring(itr->first) == true) {
            std::cout<<itr->first<<" Found\n";
			//printf("testing2");
		}
            //Call Parasail here!!
        else {
            std::cout<<itr->first<<" Not found\n";
			//printf("testing3");
		}
        itr++;
    }*/

	//std::string to_find("AATATCAAAT");
	/*
	if(find_substring(to_find) == true) 
		std::cout<<to_find<<" Found\n";
		//Call Parasail here!!
	else 
		std::cout<<to_find<<" Not found\n";
	*/
	parasail_matrix_free(matrix);
	output_file.close();
	control_fd.close();
	//fclose(control_fd);
	kseq_destroy(seq);
	gzclose(fp);
	return 0;
}

