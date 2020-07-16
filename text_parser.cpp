#include <iostream>
#include <fstream>
#include <cstring>
/*
void removeDupWord(std::string str) 
{ 
    // Used to split string around spaces. 
    std::istringstream ss(str); 
  
    // Traverse through all words 
    do { 
        // Read a word 
        std::string word; 
        ss >> word; 
  
        // Print the read word 
        std::cout << word << std::endl; 
  
        // While there is more to read 
    } while (ss); 
}*/ 

std::string removeDupWord(std::string str) 
{ 
   std::string word = ""; 
   for (auto x : str) 
   { 
       if (x == ' ') 
       { 
           std::cout << word << std::endl; 
           word = ""; 
       } 
       else
       { 
           word = word + x; 
       } 
   }  
   std::cout << word << std::endl; 
   return word;
}

int main() {
	std::ifstream control_fd("./data/Control_oligos(-)38lines.txt");

	std::string probe_no;
	std::string probe_name;
	std::string seq;
	getline(control_fd,probe_no);
	while(control_fd) {
		//getline(control_fd,probe_no);
		getline(control_fd,probe_no);
		//getline(control_fd, probe_name,' ');
		//getline(control_fd,seq,'\n');
		//seq = removeDupWord(probe_no);
		seq = probe_no.substr(75,61);
		//std::cout<<probe_no<<std::endl;
		std::cout<<seq<<std::endl;
		//std::cout<<"seq "<<seq<<"\nProbe_no "<<probe_no<<"\n Probe name "<<probe_name<<std::endl;
		//break;
	}

	control_fd.close();
	return 0;
}
