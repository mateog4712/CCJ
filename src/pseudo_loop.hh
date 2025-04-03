#ifndef PSEUDO_LOOP_H_
#define PSEUDO_LOOP_H_
#include "base_types.hh"
#include "h_struct.hh"
#include "constants.hh"
#include <stdio.h>
#include <string.h>
#include "s_energy_matrix.hh"

class VM_final;
class V_final;
class pseudo_loop{

public:
	// constructor
	pseudo_loop(std::string seq, s_energy_matrix *V, short *S, short *S1, vrna_param_t *params);

	// destructor
	~pseudo_loop();

    void compute_energies(cand_pos_t i, cand_pos_t j);

	void backtrack(minimum_fold *f, seq_interval *cur_interval);
    void set_stack_interval(seq_interval *stack_interval);
    seq_interval *get_stack_interval(){return stack_interval;}
    std::string get_structure(){return structure;}
    minimum_fold *get_minimum_fold(){return f;}

	energy_t get_WB(cand_pos_t i, cand_pos_t j);
	energy_t get_WBP(cand_pos_t i, cand_pos_t j);
	
	// nested substr in a pseudoloop
	energy_t get_WP(cand_pos_t i, cand_pos_t j);
	energy_t get_WPP(cand_pos_t i, cand_pos_t j);
	
	energy_t get_energy(cand_pos_t i, cand_pos_t j){return get_P(i,j);}
	energy_t get_P(cand_pos_t i, cand_pos_t j);
	energy_t get_PK(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	energy_t get_PL(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	energy_t get_PR(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	energy_t get_PM(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	energy_t get_PO(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	
	
	energy_t get_PfromL(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
    energy_t get_PfromR(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	energy_t get_PfromM(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	energy_t get_PfromO(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	
	
	energy_t get_PLiloop(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	energy_t get_PLiloop5(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l, cand_pos_t s);
	energy_t get_PLmloop(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	energy_t get_PLmloop0(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	energy_t get_PLmloop1(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	
	energy_t get_PRiloop(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	energy_t get_PRiloop5(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l, cand_pos_t s);
	energy_t get_PRmloop(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	energy_t get_PRmloop0(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	energy_t get_PRmloop1(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	
	energy_t get_PMiloop(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	energy_t get_PMiloop5(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l, cand_pos_t s);
	energy_t get_PMmloop(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	energy_t get_PMmloop0(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	energy_t get_PMmloop1(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	
	energy_t get_POiloop(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	energy_t get_POiloop5(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l, cand_pos_t s);
	energy_t get_POmloop(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	energy_t get_POmloop0(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	energy_t get_POmloop1(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);

	std::vector<energy_t> P;					// the main loop for pseudoloops and bands
private:

	cand_pos_t n;
	std::string res;
	std::string seq;

    s_energy_matrix *V;		        // the V object

	seq_interval *stack_interval;
	std::string structure;
	minimum_fold *f;
	vrna_param_t *params_;

    std::vector<energy_t> WI;				// the loop inside a pseudoknot (in general it looks like a W but is inside a pseudoknot)
	std::vector<cand_pos_t> index;				// the array to keep the index of two dimensional arrays like WI and weakly_closed

	short *S_;
	short *S1_;

	std::vector<energy_t> WPP;				// similar to WP but has at least one base pair
	std::vector<energy_t> WBP;				// similar to WB but has at least one base pair
	
	std::vector<std::vector<energy_t> > PK;				// MFE of a TGB structure over gapped region [i,j] U [k,l]
	std::vector<std::vector<energy_t> > PL;				// MFE of a TGB structure s.t. i.j is paired
	std::vector<std::vector<energy_t> > PR;				// MFE of a TGB structure s.t. k.l is paired
	std::vector<std::vector<energy_t> > PM;				// MFE of a TGB structure s.t. j.k is paired
	std::vector<std::vector<energy_t> > PO;				// MFE of a TGB structure s.t. i.l is paired
	
	// transition recurrences
	std::vector<std::vector<energy_t> > PfromL;
	std::vector<std::vector<energy_t> > PfromR;
	std::vector<std::vector<energy_t> > PfromM;
	std::vector<std::vector<energy_t> > PfromO;
	
	// internal loops and multi loops that span a band
	std::vector<std::vector<energy_t> > PLiloop;
	std::vector<std::vector< std::vector<energy_t> > > PLiloop5;
	std::vector<std::vector<energy_t> > PLmloop;
	std::vector<std::vector<energy_t> > PLmloop0;
	std::vector<std::vector<energy_t> > PLmloop1;
	
	
	std::vector<std::vector<energy_t> > PRiloop;
	std::vector<std::vector< std::vector<energy_t> > > PRiloop5;
	std::vector<std::vector<energy_t> > PRmloop;
	std::vector<std::vector<energy_t> > PRmloop0;
	std::vector<std::vector<energy_t> > PRmloop1;
	
	
	std::vector<std::vector<energy_t> > PMiloop;
	std::vector<std::vector< std::vector<energy_t> > > PMiloop5;
	std::vector<std::vector<energy_t> > PMmloop;
	std::vector<std::vector<energy_t> > PMmloop0;
	std::vector<std::vector<energy_t> > PMmloop1;
	
	
	std::vector<std::vector<energy_t> > POiloop;
	std::vector<std::vector< std::vector<energy_t> > > POiloop5;
	std::vector<std::vector<energy_t> > POmloop;
	std::vector<std::vector<energy_t> > POmloop0;
	std::vector<std::vector<energy_t> > POmloop1;

	void compute_WBP(cand_pos_t i, cand_pos_t l);
	void compute_WPP(cand_pos_t i, cand_pos_t l);
		
	void compute_P(cand_pos_t i, cand_pos_t l);
	void compute_PK(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	void compute_PL(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	void compute_PR(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	void compute_PM(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	void compute_PO(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	
	
	void compute_PfromL(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
    void compute_PfromR(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	void compute_PfromM(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	void compute_PfromO(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	
	
	void compute_PLiloop(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	void compute_PLiloop5(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l, cand_pos_t s);
	void compute_PLmloop(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	void compute_PLmloop0(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	void compute_PLmloop1(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	
	void compute_PRiloop(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	void compute_PRiloop5(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l, cand_pos_t s);
	void compute_PRmloop(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	void compute_PRmloop0(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	void compute_PRmloop1(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	
	void compute_PMiloop(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	void compute_PMiloop5(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l, cand_pos_t s);
	void compute_PMmloop(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	void compute_PMmloop0(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	void compute_PMmloop1(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	
	void compute_POiloop(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	void compute_POiloop5(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l, cand_pos_t s);
	void compute_POmloop(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	void compute_POmloop0(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	void compute_POmloop1(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);

    // function to allocate space for the arrays
    void allocate_space();


	energy_t get_e_stP(cand_pos_t i, cand_pos_t j);
	energy_t get_e_intP(cand_pos_t i,cand_pos_t ip, cand_pos_t jp, cand_pos_t j);
	energy_t compute_int(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l);

	energy_t alpha1P(cand_pos_t z);
	energy_t alpha2P( cand_pos_t i, cand_pos_t l);
	energy_t alpha3P(cand_pos_t z);
	energy_t beta2(cand_pos_t i, cand_pos_t l);
	energy_t beta2P(cand_pos_t i, cand_pos_t l);
	energy_t gamma2(cand_pos_t i, cand_pos_t l);

  	// Hosna: Feb 19th 2007
  	// used for backtracking
  	void insert_node (cand_pos_t i, cand_pos_t j, char type);//, seq_interval *stack_interval);

};
#endif /*PSEUDO_LOOP_H_*/
