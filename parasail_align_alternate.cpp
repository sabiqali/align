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
#include <vector>
#include <iomanip>
#include <algorithm>
#include <unordered_set>
#include <chrono>

KSEQ_INIT(gzFile, gzread)

//std::unordered_map<std::string,std::string> control_substrings;
//std::unordered_map<std::string,int> test_substrings;

inline void subString_control(const std::string& s, int sub_length,const std::string& probe_name, auto& control_substrings)  { 
    //std::unordered_map<std::string,int> temp_array;
    int n = s.length();
    //#pragma omp parallel for
    for (int i = 0; i < n; i++)  
        if(i+sub_length < n) {
            control_substrings.insert({s.substr(i, sub_length), probe_name});
        }
	
    //control_substrings = temp_array;
} 

inline void subString_test(const std::string& s, int sub_length, auto& test_substrings)  {
    int n = s.length();
    //#pragma omp parallel for
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

int main(int argc, char *argv[])  {
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

   std::unordered_map<std::string,std::string> control_substrings;
   std::unordered_map<std::string,int> test_substrings;
	
    parasail_result_t *best_result = NULL;
    parasail_matrix_t *matrix = NULL;
    
    std::string best_query,best_ref,best_comp,best_oligo;
    std::string second_best_query,second_best_ref,second_best_comp,second_best_oligo;


    matrix = parasail_matrix_create("ACGT", 5, -1);	

    int read_aligns;

    std::ofstream output_file(argv[3]);

    output_file<<std::left<<std::setw(55)<<"sequence_name"<<std::setw(35)<<"best_oligo_name"<<std::setw(35)<<"best_score"<<std::setw(35)<<"second_best_oligo_name"<<std::setw(35)<<"second_best_score"<<std::endl;

    int k = 15;
    int max = 0,second_max = 0;	
    int c = 1;
    std::string probe_name;
    getline(control_fd, one_line);
    while(control_fd) {
        getline(control_fd,one_line);
        
	std::string non_rep_section;
        if(one_line.length() > 75) {
            sequence2 = one_line.substr(one_line.find('A')); //used to separate the sequence from the probe no and name
	    non_rep_section = sequence2.substr(61,68);   //indexed used to extract the non repeating sequence from the sequence obtained above
	    
            std::string temp = one_line.substr(one_line.find('F')); //used to extract the probe name
            probe_name = temp.substr(0,temp.find('A'));
        }
        
        c++;
        subString_control(non_rep_section, k, sequence2, control_substrings);
    }
    std::cout<<"Calculated "<<c<<" control sequence.....\n";
    
    c = 0;
    while ((l = kseq_read(seq)) >= 0) {
	max = 0;
	second_max = 0;
        std::string sequence(seq->seq.s);
        
        subString_test(sequence, k, test_substrings);

        std::unordered_set<std::string> control_set;

        auto itr = control_substrings.begin();

        while (itr != control_substrings.end()) {                    
            if(find_substring(itr->first,test_substrings) == true) {
                
                /*if(std::find(control_set.begin(), control_set.end(), itr->second) == control_set.end()) {
                    control_set.push_back(itr->second);
                }*/
		control_set.insert(itr->second);
            } 
            itr++;
        }
	read_aligns = 0;
	
	for(const auto& elem: control_set){
	    //char seq_arr[elem.length() + 1];
            //strcpy(seq_arr,elem.c_str());
	    //auto t1 = std::chrono::high_resolution_clock::now();
            parasail_result_t *result = parasail_sg_trace_scan_16(seq->seq.s,sequence.length(),elem.c_str(),elem.length(),5,4,matrix);
            //auto t2 = std::chrono::high_resolution_clock::now();
	    parasail_traceback_t *traceback = parasail_result_get_traceback(result,seq->seq.s, sequence.length(), elem.c_str(), elem.length(),matrix,'|','*','*');
            //auto t3 = std::chrono::high_resolution_clock::now();
	    if(result->score > max) {
                second_max = max;
                max = result->score;
                //best_result = result;
                //second_best_probe_no = best_probe_no;
                //second_best_probe_name = best_probe_name;
                //best_probe_no = probe_no;
                //best_probe_name = probe_name;
                second_best_query = best_query;
                second_best_ref = best_ref;
                second_best_comp = best_comp;
                best_ref = traceback ->ref;
                best_comp = traceback->comp;
                best_query = traceback->query;
		second_best_oligo = best_oligo;
		best_oligo = elem;
                //std::cout<<traceback->ref<<std::endl;
            }
	    //auto duration1 = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();

	    //auto duration2 = std::chrono::duration_cast<std::chrono::microseconds>( t3 - t2 ).count();
	
	    //std::cout<<"Result Time: "<<duration1<<"\nTraceback Time: "<<duration2<<std::endl;
	    read_aligns++;
            parasail_traceback_free(traceback);
            parasail_result_free(result);
        }
        //std::cout<<"Calculating sequence "<<c<<".....\n";
	//std::cout<<"Highest Score:"<<max<<std::endl;
	//std::cout<<"Second Highest Score:"<<second_max<<std::endl;
	//output_file<<std::left<<std::setw(55)<<seq->name.s<<std::setw(35)<<best_probe_name<<std::setw(35)<<max<<std::setw(35)<<second_best_probe_name<<std::setw(35)<<second_max<<std::endl;
        c++;

	std::cout<<"\n\n\nSequence "<<c<<" : "<<seq->seq.s;
	std::cout<<"\nAligned to: "<<read_aligns<<" control oligos";
        std::cout<<"\nHighest Score: "<<max<<"\nControl Oligo: "<<best_oligo<<"\n\n";
        std::cout<<best_query<<"\n"<<best_comp<<"\n"<<best_ref;
        std::cout<<"\n\nSecond Highest Score:"<<second_max<<"\nControl Oligo: "<<second_best_oligo<<"\n\n";
        std::cout<<second_best_query<<"\n"<<second_best_comp<<"\n"<<second_best_ref<<"\n\n";

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

