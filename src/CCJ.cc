
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

bool exists(const std::string path) {
    struct stat buffer;
    return (stat(path.c_str(), &buffer) == 0);
}

// check if sequence is valid with regular expression
// check length and if any characters other than GCAUT
void validateSequence(std::string sequence) {

    if (sequence.length() == 0) {
        std::cout << "sequence is missing" << std::endl;
        exit(EXIT_FAILURE);
    }
    // return false if any characters other than GCAUT -- future implement check based on type
    for (char c : sequence) {
        if (!(c == 'G' || c == 'C' || c == 'A' || c == 'U' || c == 'T' || c == 'N')) {
            std::cout << "Sequence contains character " << c << " that is not N,G,C,A,U, or T." << std::endl;
            exit(EXIT_FAILURE);
        }
    }
}

inline void trim(std::string& s) {
    s.erase(0, s.find_first_not_of(" \t\n\r\f\v\""));
    s.erase(s.find_last_not_of(" \t\n\r\f\v\"") + 1);
}

/**
 * @brief Represents an RNA entry with name, sequence, and structure.
 *
 * This struct is designed to store information about an RNA molecule, 
 * including its name, nucleotide sequence, and secondary structure. 
 * It also provides utility functions for checking the size of the RNA entry.
 */

struct RNAEntry {
    std::string name;
    std::string sequence;

    RNAEntry() = default;

    RNAEntry(std::string rna_name, std::string rna_sequence)
        : name{rna_name},
          sequence{rna_sequence} {}

    RNAEntry(std::string rna_sequence)
        : RNAEntry("N/A", rna_sequence) {}

    size_t size() const {
        return sequence.size();
    }
};

std::vector<RNAEntry> get_all_file_entries(const std::string& file){
    if(!exists(file)){
        std::cerr << "Error: Input file not found: " << file << std::endl;
        std::exit(EXIT_FAILURE);
    }

    // state machine to parse the file
    #define UNINITIALIZED -1
    #define NAME 0
    #define SEQUENCE 1

    std::ifstream in(file.c_str());
    RNAEntry current;
    std::vector<RNAEntry> entries;
    int state = UNINITIALIZED;
    int line_number = 0;

    std::string str;
    while(getline(in, str)){
        ++line_number;
        trim(str);
        if (str.empty()) continue;

        // Check if the line is the name of the entry
        if ((state == UNINITIALIZED) && (str[0] != '>')) {
            std::cerr << "Error: Expected '>' at the beginning of the line: " << str << ". Line number: " << line_number <<std::endl;
            exit(EXIT_FAILURE);
        }

        if (str[0] == '>'){

            if (state != UNINITIALIZED) { // valid entry, save it
                if (current.sequence.empty()) {
                    std::cerr << "Warning: Sequence is empty for entry: " << current.name << ". Line number: " << line_number << ". Skipping..."<<  std::endl;
                }
                validateSequence(current.sequence);
                entries.push_back(current);
                current = RNAEntry();
            }

            current.name = str.substr(1);
            current.sequence = "";
            state = SEQUENCE;

        } else if (state == SEQUENCE){
            if(str.find('(') != std::string::npos || str.find(')') != std::string::npos || str.find('.') != std::string::npos){
                state = NAME;
            } else{
                current.sequence += str;
            }
            
        }
    }

    // Handle the last entry
    if (!current.name.empty() && !current.sequence.empty()) {
        entries.push_back(current);
    }

    return entries;
}

std::vector<RNAEntry> get_all_inputs(const std::string& fileI, const std::string& seq) {
    std::vector<RNAEntry> entries;
    if (!seq.empty()) {
        entries.emplace_back("Console Sequence", seq);
    }
    if (!fileI.empty()){
        std::vector<RNAEntry> file_entries = get_all_file_entries(fileI);
        entries.insert(entries.end(), file_entries.begin(), file_entries.end());
    }
    if (entries.empty()) throw std::runtime_error("Sequence is missing");
    return entries;
}

void print_results(std::string &fileO, std::string &sequence, std::string &structure, pf_t MFE, std::string &structure_pf, double EE, std::string MEA_structure, pf_t MEA, std::string centroid_structure, pf_t distance,std::vector<std::pair<std::string,double>> fatgraphs, pf_t frequency, pf_t diversity){
    if (fileO != "") {
        std::ofstream out(fileO,std::fstream::app);
        out << sequence << std::endl;
        out << structure << " (" << MFE << ")" << std::endl;
        out << structure_pf << " (" << EE << ")" << std::endl;
        out << MEA_structure << " (" << MEA << ")" << std::endl;
        out << centroid_structure << " (" << distance << ")" << std::endl;
        for(size_t j=0; j<fatgraphs.size();++j){
               out << fatgraphs[j].first << "\t(" << fatgraphs[j].second << ")\t";
        }
        out << std::endl;
        out << "frequency of MFE structure in ensemble: " << frequency << "; ensemble diversity " << diversity << std::endl;
    } else {
        std::cout << sequence << std::endl;
        std::cout << structure << " (" << MFE << ")" << std::endl;
        std::cout << structure_pf << " (" << EE << ")" << std::endl;
        std::cout << MEA_structure << " (" << MEA << ")" << std::endl;
        std::cout << centroid_structure << " (" << distance << ")" << std::endl;
        for(size_t j=0; j<fatgraphs.size();++j){
               std::cout << fatgraphs[j].first << "\t(" << fatgraphs[j].second << ")\t";
        }
        std::cout << std::endl;
        std::cout << "frequency of MFE structure in ensemble: " << frequency << "; ensemble diversity " << diversity << std::endl;
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

std::string ccj_pf(std::string seq,double &energy,std::string &MFE_structure, double MFE, std::string &MEA_structure, pf_t &MEA, std::string &centroid_structure,pf_t &distance, std::vector<std::pair<std::string,double>> &fatgraphs, pf_t &frequency, pf_t &diversity, int dangle, int num_samples, bool PSplot){
	W_final_pf min_fold(seq,MFE_structure,MFE,dangle,num_samples,PSplot);
    cand_pos_t n = seq.length();
	energy = min_fold.ccj_pf();
    std::string structure = min_fold.structure;
    MEA_structure = std::string(n,'.');
    MEA = 0;
    centroid_structure = std::string(n,'.');
    distance = 0;
    frequency = min_fold.frequency;
    diversity = 0;
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

    std::string fileI = args_info.input_file_given ? args_info.input_file_arg : "";

    std::string fileO = args_info.output_file_given ? args_info.output_file_arg : "";

    int num_samples = args_info.samples_arg;

    int num_fatgraph = args_info.fatgraph_arg;

    bool PSplot = !args_info.noPS_given;

    noGU = args_info.noGU_given;

	std::vector<RNAEntry> Inputs = get_all_inputs(fileI,seq);
    for(RNAEntry current : Inputs){
        std::transform(current.sequence.begin(), current.sequence.end(), current.sequence.begin(), ::toupper);
	    if(!args_info.noConv_flag) seqtoRNA(current.sequence);

        if(args_info.paramFile_given){
        std::string file = args_info.paramFile_arg;
        if (exists(file)) vrna_params_load(file.c_str(), VRNA_PARAMETER_FORMAT_DEFAULT);
            else{
                std::cerr << "Not a valid parameter file!" << std::endl;
                exit(EXIT_FAILURE);
            }
        } else {
            if (current.sequence.find('T') != std::string::npos) {
                noGU = 1;
                vrna_params_load_DNA_Mathews2004();
            } else{
                std::string file = std::string(PARAMS_DIR) + "/rna_DirksPierce09.par";
                if (exists(file)) vrna_params_load(file.c_str(), VRNA_PARAMETER_FORMAT_DEFAULT);
                else{
                    std::cerr << "Not a valid parameter file!" << std::endl;
                    exit(EXIT_FAILURE);
                }
            }
        }
        pf_t energy,pf_energy,MEA,distance,frequency,diversity;
        std::string MEA_structure,centroid_structure;
        std::vector<std::pair<std::string,double>> fatgraphs(num_fatgraph);
        std::string structure = ccj(current.sequence,energy,args_info.dangles_arg);
        std::string pf_structure = ccj_pf(current.sequence,pf_energy,structure,energy,MEA_structure,MEA,centroid_structure,distance,fatgraphs,frequency,diversity,args_info.dangles_arg,num_samples,PSplot); // I am seeing this give nan sometimes with GGGGGGAAGGGGGGGGAACCCCCCACCCCCCCC currently
        print_results(fileO,current.sequence,structure,energy,pf_structure,pf_energy,MEA_structure,MEA,centroid_structure,distance,fatgraphs,frequency,diversity);
    }

    cmdline_parser_free(&args_info);


    return 0;
}



