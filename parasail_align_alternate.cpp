#include <zlib.h>
#include <stdio.h>
#include "kseq.h"
//#include "parasail.h"
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
	
    //parasail_result_t *best_result = NULL;
    //parasail_matrix_t *matrix = NULL;
    
    std::string best_query,best_ref,best_comp,best_oligo;
    std::string second_best_query,second_best_ref,second_best_comp,second_best_oligo;


    //matrix = parasail_matrix_create("ACGT", 5, -1);	

    int read_aligns;

    std::ofstream output_file(argv[3]);

    output_file<<std::left<<std::setw(55)<<"sequence_name"<<std::setw(35)<<"best_oligo_name"<<std::setw(35)<<"best_score"<<std::setw(35)<<"second_best_oligo_name"<<std::setw(35)<<"second_best_score"<<std::endl;

    int k = 15;
    std::string best_probe_name;
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
        subString_control(non_rep_section, k, one_line, control_substrings);
    }
    std::cout<<"Calculated "<<c<<" control sequence.....\n";
    
    unsigned char * best_alignment;
    int best_alignmentLength;
    int* best_endLocations;

    std::cout<<std::left<<std::setw(55)<<"Sequence Name"<<std::setw(35)<<"Oligo Name"<<std::setw(35)<<"Reads Aligned"<<std::setw(35)<<"Score"<<std::endl;

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
	   
	    sequence2 = elem.substr(elem.find('A')); //used to separate the sequence from the probe no and name
            //non_rep_section = sequence2.substr(61,68);   //indexed used to extract the non repeating sequence from the sequence obtained above

            std::string temp = elem.substr(elem.find('F')); //used to extract the probe name
            probe_name = temp.substr(0,temp.find('A')); 
	   
	    EdlibAlignResult result = edlibAlign(seq->seq.s, sequence.length(), sequence2.c_str(), sequence2.length(), edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));
	    
	    if(result.editDistance > max) {
			second_max = max;
			max = result.editDistance;
			//best_result = result;
			//second_best_probe_no = best_probe_no;
			//second_best_probe_name = best_probe_name;
			//best_probe_no = probe_no;
			best_probe_name = probe_name;
			//second_best_query = best_query;
			//second_best_ref = best_ref;
			//second_best_comp = best_comp;
			//best_ref = traceback ->ref;
			//best_comp = traceback->comp;
			//best_query = traceback->query;
			second_best_oligo = best_oligo;
			best_oligo = sequence2;
			best_alignment = result.alignment;
            		best_alignmentLength = result.alignmentLength;
            		best_endLocations = result.endLocations;
			//std::cout<<traceback->ref<<std::endl;
            }
	    //auto duration1 = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();

	    //auto duration2 = std::chrono::duration_cast<std::chrono::microseconds>( t3 - t2 ).count();
	
	    //std::cout<<"Result Time: "<<duration1<<"\nTraceback Time: "<<duration2<<std::endl;
	    read_aligns++;
	    edlibFreeAlignResult(result);
		//parasail_traceback_free(traceback);
		//parasail_result_free(result);
	}
        //std::cout<<"Calculating sequence "<<c<<".....\n";
	//std::cout<<"Highest Score:"<<max<<std::endl;
	//std::cout<<"Second Highest Score:"<<second_max<<std::endl;
	//output_file<<std::left<<std::setw(55)<<seq->name.s<<std::setw(35)<<best_probe_name<<std::setw(35)<<max<<std::setw(35)<<second_best_probe_name<<std::setw(35)<<second_max<<std::endl;
        c++;

	//std::cout<<"\n\n\nSequence "<<c<<" : "<<seq->seq.s;
	//std::cout<<"\nAligned to: "<<read_aligns<<" control oligos";
        //std::cout<<"\nHighest Score: "<<max<<"\nControl Oligo: "<<best_oligo<<"\n\n";
        //std::cout<<best_query<<"\n"<<best_comp<<"\n"<<best_ref;
        //std::cout<<"\n\nSecond Highest Score:"<<second_max<<"\nControl Oligo: "<<second_best_oligo<<"\n\n";
        //std::cout<<second_best_query<<"\n"<<second_best_comp<<"\n"<<second_best_ref<<"\n\n";
	//printAlignment(seq->seq.s, best_oligo.c_str(), best_alignment, best_alignmentLength, *(best_endLocations), EDLIB_MODE_NW);
        std::cout<<std::left<<std::setw(55)<<seq->name.s<<std::setw(35)<<best_probe_name<<std::setw(35)<<read_aligns<<std::setw(35)<<max<<std::endl;
	
	/*output_file<<"\n\n\nSequence "<<c<<" : "<<seq->seq.s;
        output_file<<"\nHighest Score: "<<max<<"\nProbe Name:"<<best_probe_no;
        output_file<<best_query<<"\n"<<best_comp<<"\n"<<best_ref;
        output_file<<"\n\nSecond Highest Score:"<<second_max<<"\nProbe name:"<<second_best_probe_no<<std::endl;
        output_file<<second_best_query<<"\n"<<second_best_comp<<"\n"<<second_best_ref<<"\n\n";*/
	test_substrings.erase(test_substrings.begin(),test_substrings.end());
    }

    //parasail_matrix_free(matrix);
    //fclose(out_file);
    control_fd.close();
    //fclose(control_fd);
    kseq_destroy(seq);
    gzclose(fp);
    return 0;
}
