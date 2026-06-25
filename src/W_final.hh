#ifndef W_FINAL_H_
#define W_FINAL_H_

#include "pseudo_loop.hh"
#include "base_types.hh"
#include "s_energy_matrix.hh"
#include "constants.hh"
#include <cstring>
#include <string>
#include <vector>

#include "ViennaRNA/loops.hh"
#include "ViennaRNA/pair_mat.hh"
#include "ViennaRNA/params/io.hh"


class W_final{
	public:
		W_final(std::string seq, int dangle);
        // constructor for the restricted mfe case

        ~W_final ();
        // The destructor

        // double hfold (sparse_tree &tree);
        double ccj ();

        vrna_param_t *params_;
        std::string structure;        // MFE structure

        void Trace_W(cand_pos_t i, cand_pos_t j, energy_t e);
        void Trace_V(cand_pos_t i, cand_pos_t j, energy_t e);
        void Trace_WM(cand_pos_t i, cand_pos_t j, energy_t e);
        void Trace_WMv(cand_pos_t i, cand_pos_t j, energy_t e);
        void Trace_WMp(cand_pos_t i, cand_pos_t j, energy_t e);


    protected:
    	// Hosna: June 18th, 2007:
        // this pointer is the main part of the Hierarchical fold program
        // and corresponds to WMB recurrence
        pseudo_loop *P;

        s_energy_matrix *V;     // the V object
        std::vector<energy_t> W; // size n+1 so left as non-Trianglematrix
        // PARAMTYPE *W;                 // the W exterior loop array
        cand_pos_t n;     // sequence length (number of nucleotides)
        seq_interval *stack_interval;  // used for backtracking
        minimum_fold *f;        // the minimum folding, see structs.h
        std::string seq_;
        short *S_;
	    short *S1_;
        

        void insert_node (cand_pos_t i, cand_pos_t j, char type);

        void space_allocation();

        // allocate the necessary memory
        double fold_sequence_restricted ();

        void backtrack();

        energy_t E_ext_Stem(const energy_t& vij,const energy_t& vi1j,const energy_t& vij1,const energy_t& vi1j1,const short* S, vrna_param_t* params, const cand_pos_t i,const cand_pos_t j, cand_pos_t n);

        void fill_structure();
};

#endif /*W_FINAL_H_*/
