
// a simple driver for CCJ
#include "cmdline.hh"
#include "W_final.hh"
#include "part_func.hh"
#include "h_globals.hh"
//
#include <iostream>
#include <fstream>
#include <algorithm>
#include <stdio.h>
#include <sys/stat.h>
#include <string>
#include <getopt.h>

bool exists (const std::string path) {
  struct stat buffer;   
  return (stat (path.c_str(), &buffer) == 0); 
}

//check if sequence is valid with regular expression
//check length and if any characters other than GCAUT
void validateSequence(std::string sequence){

	if(sequence.length() == 0){
		std::cout << "sequence is missing" << std::endl;
		exit(EXIT_FAILURE);
	}
  // return false if any characters other than GCAUT -- future implement check based on type
  for(char c : sequence) {
    if (!(c == 'G' || c == 'C' || c == 'A' || c == 'U' || c == 'T')) {
		std::cout << "Sequence contains character " << c << " that is not G,C,A,U, or T." << std::endl;
		exit(EXIT_FAILURE);
    }
  }
}

void seqtoRNA(std::string &sequence){
    for (char &c : sequence) {
      	if (c == 'T') c = 'U';
    }
}

std::string ccj(std::string seq,double &energy, int dangle){
	W_final min_fold(seq, dangle);
	energy = min_fold.ccj();
    std::string structure = min_fold.structure;
    return structure;
}

std::string ccj_pf(std::string seq,pf_t &energy,std::string mfe_structure,double mfe_energy, int dangle, int num_samples,bool PSplot){
	W_final_pf min_fold(seq,mfe_structure,mfe_energy,dangle,num_samples,PSplot);
	energy = min_fold.ccj_pf();
    std::string structure = min_fold.structure;
    return structure;
}

int main (int argc, char *argv[])
{
    args_info args_info;

	// get options (call getopt command line parser)
	if (cmdline_parser (argc, argv, &args_info) != 0) {
	    exit(1);
	}

    std::string seq;
	if (args_info.inputs_num>0) {
	    seq=args_info.inputs[0];
	} else {
		if(!args_info.input_file_given) std::getline(std::cin,seq);
	}

    int dangles = args_info.dangles_given ? dangle_model : 2;

    std::transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
	if(!args_info.noConv_given) seqtoRNA(seq);

    validateSequence(seq);

    std::string file= args_info.paramFile_given ? parameter_file : "params/rna_DirksPierce09.par";
	if(exists(file)){
		vrna_params_load(file.c_str(), VRNA_PARAMETER_FORMAT_DEFAULT);
	}
	else if (seq.find('T') != std::string::npos){
		vrna_params_load_DNA_Mathews2004();
	}
	int num_samples = 1000; 
	bool PSplot = true;
	double energy;
	pf_t pf_energy;
    std::string structure = ccj(seq,energy,dangles);
	std::string pf_structure = ccj_pf(seq,pf_energy,structure,energy,dangles,num_samples,PSplot);

    std::cout << seq << std::endl;
    std::cout << structure << " (" << energy << ")" << std::endl;
	std::cout << pf_structure << " (" << pf_energy << ")" << std::endl;

    cmdline_parser_free(&args_info);


    return 0;
}



