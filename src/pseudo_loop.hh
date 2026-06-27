#ifndef PSEUDO_LOOP_H_
#define PSEUDO_LOOP_H_
#include "base_types.hh"
#include "h_struct.hh"
#include "constants.hh"
#include <stdio.h>
#include <string.h>
#include "s_energy_matrix.hh"
#include "matrices.hh"
// #define debug 0

#ifdef NDEBUG
	#define UNREACHABLE() \
		do { \
			std::cerr << "Reached unreachable at line " << __LINE__ << " in File: " << __FILE__ << std::endl; \
			abort(); \
		} while(0)
#else
    #define UNREACHABLE() __builtin_unreachable()
#endif

class W_final;
class pseudo_loop{

public:
	// constructor
	pseudo_loop(std::string seq, s_energy_matrix *V, W_final *W, short *S, short *S1, vrna_param_t *params);

	// destructor
	~pseudo_loop();

    void compute_energies(cand_pos_t i, cand_pos_t j);

	void set_fold(std::vector<int> &fres);
	
	energy_t get_energy(cand_pos_t i, cand_pos_t j){return P.get(i,j);}

	TriangleMatrix P;					// the main loop for pseudoloops and bands
	void Trace_P(cand_pos_t i, cand_pos_t l, energy_t e);
private:

	cand_pos_t n;
	std::string res;
	std::string seq;

    s_energy_matrix *V;		        // the V object
	W_final *W;		        // the W object

	std::string structure;
	std::vector<int> *fres;
	vrna_param_t *params_;

	std::vector<cand_pos_t> index;				// the array to keep the index of two dimensional arrays like WPP and WBP
	index_offset_t index3D;

	short *S_;
	short *S1_;

	TriangleMatrix WPP;	// similar to WP but has at least one base pair
	TriangleMatrix WBP;	// similar to WB but has at least one base pair
	
	Matrix4D PL;		// MFE of a TGB structure s.t. i.j is paired
	Matrix4D PLmloop00;
	Matrix4D PLmloop01;
	Matrix4D PLmloop10;
	Matrix4D PfromL;
	Matrix4D PfromLprime;

	Matrix4D PR;		// MFE of a TGB structure s.t. k.l is paired
	Matrix4D PRmloop00;
	Matrix4D PRmloop01;
	Matrix4D PRmloop10;
	Matrix4D PfromR;
	Matrix4D PfromRprime;

	Matrix4D PM;		// MFE of a TGB structure s.t. j.k is paired
	Matrix4D PfromM;
	Matrix4D PfromMprime;
	Matrix4D PMmloop00;
	Matrix4D PMmloop01;
	Matrix4D PMmloop10;

	Matrix4D POm;	// MFE of a TGB structure s.t. i.l is paired multiple
	Matrix4D POmmloop00;
	Matrix4D POmmloop01;
	Matrix4D POmmloop10;
	Matrix4D PfromO;
	Matrix4D PfromOprime;
	
	Matrix4D POs;	// MFE of a TGB structure s.t. i.l is paired single
	Matrix4D POsmloop00;
	Matrix4D POsmloop01;
	Matrix4D POsmloop10;

	Matrix4D PLreO;
	Matrix4D PLreOmloop00;
	Matrix4D PLreOmloop01;
	Matrix4D PLreOmloop10;
	Matrix4D PfromLreO;
	Matrix4D PfromLreOprime;

	Matrix4D PLreR;
	Matrix4D PLreRmloop00;
	Matrix4D PLreRmloop01;
	Matrix4D PLreRmloop10;
	Matrix4D PfromLreR;
	Matrix4D PfromLreRprime;

	Matrix4D PMreO;
	Matrix4D PMreOmloop00;
	Matrix4D PMreOmloop01;
	Matrix4D PMreOmloop10;
	Matrix4D PfromMreO;
	Matrix4D PfromMreOprime;

	Matrix4D PMreR;
	Matrix4D PMreRmloop00;
	Matrix4D PMreRmloop01;
	Matrix4D PMreRmloop10;
	Matrix4D PfromMreR;
	Matrix4D PfromMreRprime;

	Matrix4D PLMreR;
	Matrix4D PLMreRmloop00;
	Matrix4D PLMreRmloop01;
	Matrix4D PLMreRmloop10;
	Matrix4D PfromLMreR;
	Matrix4D PfromLMreRprime;

	Matrix4D PLMorO;
	Matrix4D PLMorOmloop00;
	Matrix4D PLMorOmloop01;
	Matrix4D PLMorOmloop10;
	Matrix4D PfromLMorO;
	Matrix4D PfromLMorOprime;

	Matrix4D PK1Om;		// MFE of a TGB structure over gapped region [i,j] U [k,l]
	Matrix4D PK2Om;

	Matrix4D PK1Os;
	Matrix4D PK2Os;

	Matrix4D PK1LreO;
	Matrix4D PK2LreO;

	Matrix4D PK1MreO;
	Matrix4D PK2MreO;

	Matrix4D PK1LreR;
	Matrix4D PK2LreR;

	Matrix4D PK1MreR;
	Matrix4D PK2MreR;

	Matrix4D PK1R;
	Matrix4D PK2R;

	Matrix4D PK1LMreR;
	Matrix4D PK1LMorO;

	void compute_WBP(cand_pos_t i, cand_pos_t l);
	void compute_WPP(cand_pos_t i, cand_pos_t l);
		
	void compute_P(cand_pos_t i, cand_pos_t l);
	void compute_PK1X(const Index4D &x, MType type);
	void compute_PK2X(const Index4D &x, MType type);
	void compute_PX(const Index4D &x, MType type);
	void compute_PfromX(const Index4D &x, MType type);
	void compute_PfromXprime(const Index4D &x, MType type);
	energy_t calc_PfromXdoubleprime(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l, MType type);
	void compute_PXmloop00(const Index4D &x, MType type);
	void compute_PXmloop01(const Index4D &x, MType type);
	void compute_PXmloop10(const Index4D &x, MType type);

	energy_t calc_PLiloop(const Index4D &x, MType type);
	energy_t calc_PRiloop(const Index4D &x, MType type);
	energy_t calc_PMiloop(const Index4D &x, MType type);
	energy_t calc_POiloop(const Index4D &x, MType type);

	void compute_PLmloop00(const Index4D &x, MType type);
	void compute_PLmloop01(const Index4D &x, MType type);
	void compute_PLmloop10(const Index4D &x, MType type);
	void compute_PMmloop00(const Index4D &x, MType type);
	void compute_PMmloop01(const Index4D &x, MType type);
	void compute_PMmloop10(const Index4D &x, MType type);
	void compute_PRmloop00(const Index4D &x, MType type);
	void compute_PRmloop01(const Index4D &x, MType type);
	void compute_PRmloop10(const Index4D &x, MType type);
	void compute_POmloop00(const Index4D &x, MType type);
	void compute_POmloop01(const Index4D &x, MType type);
	void compute_POmloop10(const Index4D &x, MType type);

	void compute_PfromL(const Index4D &x, MType type);
	void compute_PfromM(const Index4D &x,MType type);
	void compute_PfromR(const Index4D &x,MType type);
	void compute_PfromO(const Index4D &x,MType type);
	
	void compute_PfromLprime(const Index4D &x,MType type);
	void compute_PfromMprime(const Index4D &x,MType type);
	void compute_PfromRprime(const Index4D &x,MType type);
	void compute_PfromOprime(const Index4D &x,MType type);

	inline energy_t calc_PO(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l){return std::min(POm.get(i,j,k,l),POs.get(i,j,k,l));};
	energy_t calc_WB(cand_pos_t i, cand_pos_t l);
	energy_t calc_WP(cand_pos_t i, cand_pos_t l);
	energy_t calc_PfromLdoubleprime(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	energy_t calc_PfromRdoubleprime(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	energy_t calc_PfromMdoubleprime(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	energy_t calc_PfromOdoubleprime(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	energy_t calc_PfromLreOdoubleprime(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	energy_t calc_PfromLreRdoubleprime(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	energy_t calc_PfromMreOdoubleprime(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	energy_t calc_PfromMreRdoubleprime(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	energy_t calc_PfromLMreRdoubleprime(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	energy_t calc_PfromLMorOdoubleprime(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);


	// Traceback //
	void Trace_WB(cand_pos_t i, cand_pos_t l, energy_t e);
	void Trace_WBP(cand_pos_t i, cand_pos_t l, energy_t e);
	void Trace_WP(cand_pos_t i, cand_pos_t l, energy_t e);
	void Trace_WPP(cand_pos_t i, cand_pos_t l, energy_t e);

	void Trace_PX1(cand_pos_t i,cand_pos_t j,cand_pos_t k, cand_pos_t l, MType type, energy_t e);
	void Trace_PX2(cand_pos_t i,cand_pos_t j,cand_pos_t k, cand_pos_t l, MType type, energy_t e);
	void Trace_PX(cand_pos_t i,cand_pos_t j,cand_pos_t k, cand_pos_t l, MType type, energy_t e);
	void Trace_PXiloop(const Index4D &x, MType type, energy_t e);
	void Trace_PXmloop(const Index4D &x, MType type, energy_t e);
	void Trace_PXmloop00(const Index4D &x, MType type, energy_t e);
	void Trace_PXmloop01(const Index4D &x, MType type, energy_t e);
	void Trace_PXmloop10(const Index4D &x, MType type, energy_t e);
	void Trace_PfromX(const Index4D &x, MType type, energy_t e);
	void Trace_PfromXprime(cand_pos_t i,cand_pos_t j,cand_pos_t k, cand_pos_t l, MType type, energy_t e);
	void Trace_PfromXdoubleprime(cand_pos_t i,cand_pos_t j,cand_pos_t k, cand_pos_t l, MType type, energy_t e);

	void Trace_PLiloop(const Index4D &x, MType type, energy_t e);
	void Trace_PMiloop(const Index4D &x, MType type, energy_t e);
	void Trace_PRiloop(const Index4D &x, MType type, energy_t e);
	void Trace_POiloop(const Index4D &x, MType type, energy_t e);

	void Trace_PLmloop00(const Index4D &x, MType type, energy_t e);
	void Trace_PMmloop00(const Index4D &x, MType type, energy_t e);
	void Trace_PRmloop00(const Index4D &x, MType type, energy_t e);
	void Trace_POmloop00(const Index4D &x, MType type, energy_t e);

	void Trace_PLmloop01(const Index4D &x, MType type, energy_t e);
	void Trace_PMmloop01(const Index4D &x, MType type, energy_t e);
	void Trace_PRmloop01(const Index4D &x, MType type, energy_t e);
	void Trace_POmloop01(const Index4D &x, MType type, energy_t e);

	void Trace_PLmloop10(const Index4D &x, MType type, energy_t e);
	void Trace_PMmloop10(const Index4D &x, MType type, energy_t e);
	void Trace_PRmloop10(const Index4D &x, MType type, energy_t e);
	void Trace_POmloop10(const Index4D &x, MType type, energy_t e);

	void Trace_PfromL(cand_pos_t i,cand_pos_t j,cand_pos_t k, cand_pos_t l,MType type, energy_t e);
	void Trace_PfromM(cand_pos_t i,cand_pos_t j,cand_pos_t k, cand_pos_t l,MType type, energy_t e);
	void Trace_PfromR(cand_pos_t i,cand_pos_t j,cand_pos_t k, cand_pos_t l,MType type, energy_t e);
	void Trace_PfromO(cand_pos_t i,cand_pos_t j,cand_pos_t k, cand_pos_t l,MType type, energy_t e);

	void Trace_PfromLprime(cand_pos_t i,cand_pos_t j,cand_pos_t k, cand_pos_t l,MType type, energy_t e);
	void Trace_PfromMprime(cand_pos_t i,cand_pos_t j,cand_pos_t k, cand_pos_t l,MType type, energy_t e);
	void Trace_PfromRprime(cand_pos_t i,cand_pos_t j,cand_pos_t k, cand_pos_t l,MType type, energy_t e);
	void Trace_PfromOprime(cand_pos_t i,cand_pos_t j,cand_pos_t k, cand_pos_t l,MType type, energy_t e);

	void Trace_PfromLdoubleprime(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l, energy_t e);
	void Trace_PfromRdoubleprime(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l, energy_t e);
	void Trace_PfromMdoubleprime(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l, energy_t e);
	void Trace_PfromOdoubleprime(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l, energy_t e);
	void Trace_PfromLreOdoubleprime(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l, energy_t e);
	void Trace_PfromLreRdoubleprime(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l, energy_t e);
	void Trace_PfromMreOdoubleprime(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l, energy_t e);
	void Trace_PfromMreRdoubleprime(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l, energy_t e);
	void Trace_PfromLMreRdoubleprime(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l, energy_t e);
	void Trace_PfromLMorOdoubleprime(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l, energy_t e);


    // function to allocate space for the arrays
    void allocate_space();

	energy_t get_e_stP(cand_pos_t i, cand_pos_t j);
	energy_t get_e_intP(cand_pos_t i,cand_pos_t ip, cand_pos_t jp, cand_pos_t j);
	energy_t compute_int(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l);

	// penalty for closing pair i.l or l.i of a pseudoloop
	static constexpr energy_t gamma2(cand_pos_t i, cand_pos_t l){
		return 0;
	}
	template<class Penalty> energy_t penalty(const Index4D &x, Penalty p, MType type);
	energy_t calc_PXiloop(const Index4D &x, MType type);
	energy_t calc_PXmloop(const Index4D &x, MType type);
	inline bool impossible_case(const Index4D &x) const {
    	return !x.is_valid(n);
	}
	inline Matrix4D& PX_by_mtype(MType type) {
		static std::array<Matrix4D*,11> matrices{&PL, &PM, &PR, &POm, &POs, &PLreO, &PLreR, &PMreO, &PMreR, &PLMreR, &PLMorO};
		return *matrices[static_cast<int>(type)];
	}
	inline Matrix4D& PfromX_by_mtype(MType type) {
		static std::array<Matrix4D*,11> matrices{&PfromL, &PfromM, &PfromR, &PfromO, &PfromO, &PfromLreO, &PfromLreR, &PfromMreO, &PfromMreR, &PfromLMreR, &PfromLMorO};
		return *matrices[static_cast<int>(type)];
	}
	inline Matrix4D& PfromXprime_by_mtype(MType type) {
		static std::array<Matrix4D*,11> matrices{&PfromLprime, &PfromMprime, &PfromRprime, &PfromOprime, &PfromOprime, &PfromLreOprime, &PfromLreRprime, &PfromMreOprime, &PfromMreRprime, &PfromLMreRprime, &PfromLMorOprime};
		return *matrices[static_cast<int>(type)];
	}
	inline Matrix4D& PXmloop00_by_mtype(MType type) {
		static std::array<Matrix4D*,11> matrices{&PLmloop00, &PMmloop00, &PRmloop00, &POmmloop00, &POsmloop00, &PLreOmloop00, &PLreRmloop00, &PMreOmloop00, &PMreRmloop00, &PLMreRmloop00, &PLMorOmloop00};
		return *matrices[static_cast<int>(type)];
	}
	inline Matrix4D& PXmloop01_by_mtype(MType type) {
		static std::array<Matrix4D*,11> matrices{&PLmloop01, &PMmloop01, &PRmloop01, &POmmloop01, &POsmloop01, &PLreOmloop01, &PLreRmloop01, &PMreOmloop01, &PMreRmloop01, &PLMreRmloop01, &PLMorOmloop01};
		return *matrices[static_cast<int>(type)];
	}
	inline Matrix4D& PXmloop10_by_mtype(MType type) {
		static std::array<Matrix4D*,11> matrices{&PLmloop10, &PMmloop10, &PRmloop10, &POmmloop10, &POsmloop10, &PLreOmloop10, &PLreRmloop10, &PMreOmloop10, &PMreRmloop10, &PLMreRmloop10, &PLMorOmloop10};
		return *matrices[static_cast<int>(type)];
	}
	/**
	 * Because of the ordering of Mtype, L and M (the first two are the only two missing). This means that I can just pass Mtype-2 and get the matrices I want
	 */
	inline Matrix4D& PK1X_by_mtype(MType type) {
		static std::array<Matrix4D*,9> matrices{&PK1R, &PK1Om, &PK1Os, &PK1LreO, &PK1LreR, &PK1MreO, &PK1MreR, &PK1LMreR, &PK1LMorO};
		return *matrices[static_cast<int>(type)-2];
	}
	/**
	 * This follows the same logic as PK1 but it is missing the last two
	 */
	inline Matrix4D& PK2X_by_mtype(MType type) {
		static std::array<Matrix4D*,7> matrices{&PK2R, &PK2Om, &PK2Os, &PK2LreO, &PK2LreR, &PK2MreO, &PK2MreR};
		return *matrices[static_cast<int>(type)-2];
	}
};
#endif /*PSEUDO_LOOP_H_*/
