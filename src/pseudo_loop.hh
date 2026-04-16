#ifndef PSEUDO_LOOP_H_
#define PSEUDO_LOOP_H_
#include "base_types.hh"
#include "h_struct.hh"
#include "constants.hh"
#include <stdio.h>
#include <string.h>
#include "s_energy_matrix.hh"
#include "matrices.hh"

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
	
	// nested substr in a pseudoloop
	energy_t get_WP(cand_pos_t i, cand_pos_t j);
	
	energy_t get_energy(cand_pos_t i, cand_pos_t j){return P.get(i,j);}

	energy_t get_PfromLdoubleprime(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	energy_t get_PfromRdoubleprime(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	energy_t get_PfromMdoubleprime(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	energy_t get_PfromOdoubleprime(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	
	
	energy_t get_PLiloop(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	energy_t get_PLmloop(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	
	energy_t get_PRiloop(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	energy_t get_PRmloop(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	
	energy_t get_PMiloop(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	energy_t get_PMmloop(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	
	energy_t get_POiloop(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	energy_t get_POmloop(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);

	TriangleMatrix P;					// the main loop for pseudoloops and bands
private:

	cand_pos_t n;
	std::string res;
	std::string seq;

    s_energy_matrix *V;		        // the V object

	seq_interval *stack_interval;
	std::string structure;
	minimum_fold *f;
	vrna_param_t *params_;

	std::vector<cand_pos_t> index;				// the array to keep the index of two dimensional arrays like WI and weakly_closed
	index_offset_t index3D;

	short *S_;
	short *S1_;

	TriangleMatrix WPP;	// similar to WP but has at least one base pair
	TriangleMatrix WBP;	// similar to WB but has at least one base pair
	
	Matrix4D PK;		// MFE of a TGB structure over gapped region [i,j] U [k,l]
	Matrix4D PL;		// MFE of a TGB structure s.t. i.j is paired
	Matrix4D PR;		// MFE of a TGB structure s.t. k.l is paired
	Matrix4D PM;		// MFE of a TGB structure s.t. j.k is paired
	Matrix4D PO;		// MFE of a TGB structure s.t. i.l is paired
	// Matrix4D POs;	// MFE of a TGB structure s.t. i.l is paired single
	// Matrix4D POm;	// MFE of a TGB structure s.t. i.l is paired multiple
	
	// transition recurrences
	Matrix4D PfromL;
	Matrix4D PfromLprime;
	Matrix4D PfromR;
	Matrix4D PfromRprime;
	Matrix4D PfromM;
	Matrix4D PfromMprime;
	Matrix4D PfromO;
	Matrix4D PfromOprime;
	
	// internal loops and multi loops that span a band
	Matrix4D PLmloop00;
	Matrix4D PLmloop01;
	Matrix4D PLmloop10;
	
	
	Matrix4D PRmloop00;
	Matrix4D PRmloop01;
	Matrix4D PRmloop10;
	
	
	Matrix4D PMmloop00;
	Matrix4D PMmloop01;
	Matrix4D PMmloop10;
	
	
	Matrix4D POmloop00;
	Matrix4D POmloop01;
	Matrix4D POmloop10;

	std::vector<char> can_pair_;
    void init_can_pair() {
        can_pair_.resize((n+1)*(n+1));
        for (cand_pos_t i=1; i<=n; i++) {
            for (cand_pos_t j=i; j<=i+TURN && j<=n; j++) {
                cand_pos_t ij = i*(n+1)+j;
                can_pair_[ij] = false;
            }
            for (cand_pos_t j=i+TURN+1; j<=n; j++) {
                cand_pos_t ij = i*(n+1)+j;
                can_pair_[ij] = (pair[S_[i]][S_[j]]>0);
            }
        }
    }

	//! @brief check whether two positions can pair
    //! checks for canonical base pairing *and* distance (TURN)
    inline int can_pair(cand_pos_t i, cand_pos_t j) const {
        assert(i<=j);
        return can_pair_[i*(n+1)+j];
    }

	void compute_WBP(cand_pos_t i, cand_pos_t l);
	void compute_WPP(cand_pos_t i, cand_pos_t l);
		
	void compute_P(cand_pos_t i, cand_pos_t l);
	void compute_PK(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	void compute_PX(Index4D &x, MType type);
	energy_t compute_PX_helper(Index4D &x, MType type);
	void compute_PL(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	void compute_PR(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	void compute_PM(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	void compute_PO(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	
	
	void compute_PfromL(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
    void compute_PfromR(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	void compute_PfromM(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	void compute_PfromMprime(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	void compute_PfromO(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	
	void compute_PLmloop00(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	void compute_PLmloop01(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	void compute_PLmloop10(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	
	void compute_PRmloop00(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	void compute_PRmloop01(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	void compute_PRmloop10(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	
	void compute_PMmloop00(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	void compute_PMmloop01(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	void compute_PMmloop10(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	
	void compute_POmloop00(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	void compute_POmloop01(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	void compute_POmloop10(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);

	energy_t calc_PXiloop(const Index4D &x, MType type);
	energy_t calc_PXmloop(const Index4D &x, MType type);
	energy_t calc_PfromX(const Index4D &x, MType type);

    // function to allocate space for the arrays
    void allocate_space();

	energy_t get_e_stP(cand_pos_t i, cand_pos_t j);
	energy_t get_e_intP(cand_pos_t i,cand_pos_t ip, cand_pos_t jp, cand_pos_t j);
	energy_t compute_int(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l);

	inline energy_t beta2(cand_pos_t i, cand_pos_t l);
	inline energy_t beta2P(cand_pos_t i, cand_pos_t l);
	// penalty for closing pair i.l or l.i of a pseudoloop
	constexpr energy_t gamma2(cand_pos_t i, cand_pos_t l){
		return 0;
	}
	inline bool impossible_case(Index4D &x) const;
	inline Matrix4D& PX_by_mtype(MType type);

  	// Hosna: Feb 19th 2007
  	// used for backtracking
  	void insert_node (cand_pos_t i, cand_pos_t j, char type);//, seq_interval *stack_interval);
	void insert_node(int i, int j, int k, int l, char type);

};
#endif /*PSEUDO_LOOP_H_*/
