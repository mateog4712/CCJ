#ifndef PART_FUNC
#define PART_FUNC
#include "base_types.hh"
#include "matrices.hh"
#include <cstring>
#include <string>
#include <unordered_map>
#include <vector>

extern "C" {
#include "ViennaRNA/loops/all.h"
#include "ViennaRNA/pair_mat.h"
#include "ViennaRNA/params/default.h"
#include "ViennaRNA/params/io.h"
}

struct SzudzikHash {
    cand_pos_t operator()(const std::pair<cand_pos_t, cand_pos_t> pair) const {
        cand_pos_t a = pair.first;
        cand_pos_t b = pair.second;
        return (a >= b) ? (a * a + a + b) : (b * b + a);
    }
};

inline cand_pos_t boustrophedon_at(cand_pos_t start, cand_pos_t end, cand_pos_t pos);
std::vector<cand_pos_t> boustrophedon(cand_pos_t start, cand_pos_t end);

class W_final_pf {

  public:
    std::string structure;
    int num_samples;
    std::unordered_map<std::string, int> structures;

    W_final_pf(std::string &seq, std::string &MFE_structure, double MFE_energy, int dangle, int num_samples, bool PSplot);
    // constructor for the restricted mfe case

    ~W_final_pf();
    // The destructor

    pf_t ccj_pf();

    vrna_exp_param_t *exp_params_;

    pf_t get_WB(cand_pos_t i, cand_pos_t j);
	pf_t get_WBP(cand_pos_t i, cand_pos_t j);
	
	// nested substr in a pseudoloop
	pf_t get_WP(cand_pos_t i, cand_pos_t j);
	pf_t get_WPP(cand_pos_t i, cand_pos_t j);
	
	pf_t get_energy(cand_pos_t i, cand_pos_t j) {
        return V.get(i,j);
    }
	
	pf_t get_PLiloop(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	pf_t get_PLmloop(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	
	pf_t get_PRiloop(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	pf_t get_PRmloop(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	
	pf_t get_PMiloop(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	pf_t get_PMmloop(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	
	pf_t get_POiloop(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	pf_t get_POmloop(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);

	TriangleMatrix_PF P; // the main loop for pseudoloops and bands

  private:
    std::string seq;
    std::string MFE_structure;
    bool PSplot;
    cand_pos_t n;
    std::vector<cand_pos_t> index;
	index_offset_t index3D;

    short *S_;
    short *S1_;

    // PK Free
    TriangleMatrix_PF V;
    TriangleMatrix_PF VM;
    TriangleMatrix_PF WMv;
    TriangleMatrix_PF WMp;
    TriangleMatrix_PF WM;
    std::vector<pf_t> W;
    // PK
    TriangleMatrix_PF WPP;								// similar to WP but has at least one base pair
	TriangleMatrix_PF WBP;								// similar to WB but has at least one base pair
	
	Matrix4DPF PK;					// MFE of a TGB structure over gapped region [i,j] U [k,l]
	Matrix4DPF PL;					// MFE of a TGB structure s.t. i.j is paired
	Matrix4DPF PR;					// MFE of a TGB structure s.t. k.l is paired
	Matrix4DPF PM;					// MFE of a TGB structure s.t. j.k is paired
	Matrix4DPF PO;					// MFE of a TGB structure s.t. i.l is paired
	
	// transition recurrences
	Matrix4DPF PfromL;
	Matrix4DPF PfromR;
	Matrix4DPF PfromM;
	Matrix4DPF PfromO;
	
	// internal loops and multi loops that span a band
	Matrix4DPF PLmloop00;
	Matrix4DPF PLmloop01;
	Matrix4DPF PLmloop10;
	
	Matrix4DPF PRmloop00;
	Matrix4DPF PRmloop01;
	Matrix4DPF PRmloop10;
	
	Matrix4DPF PMmloop00;
	Matrix4DPF PMmloop01;
	Matrix4DPF PMmloop10;
	
	Matrix4DPF POmloop00;
	Matrix4DPF POmloop01;
	Matrix4DPF POmloop10;

    // Extra
    std::vector<pf_t> scale;
    std::vector<pf_t> expMLbase;
    std::vector<pf_t> expcp_pen;
    std::vector<pf_t> expPUP_pen;

    std::unordered_map<std::pair<cand_pos_t, cand_pos_t>, cand_pos_t, SzudzikHash> samples;

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

    double to_Energy(pf_t energy, cand_pos_t length);
    void rescale_pk_globals();

    void exp_params_rescale(double mfe);

    void compute_energy(cand_pos_t i, cand_pos_t j);

    void compute_WMv_WMp(cand_pos_t i, cand_pos_t j);

    void compute_energy_WM(cand_pos_t i, cand_pos_t j);

	pf_t compute_energy_VM(cand_pos_t i, cand_pos_t j);

    void compute_pk_energies(cand_pos_t i, cand_pos_t l);


    pf_t exp_Extloop(cand_pos_t i, cand_pos_t j);

    pf_t exp_MLstem(cand_pos_t i, cand_pos_t j);

    pf_t exp_Mbloop(cand_pos_t i, cand_pos_t j);

    pf_t HairpinE(cand_pos_t i, cand_pos_t j);

    pf_t compute_internal(cand_pos_t i, cand_pos_t j);

    pf_t compute_int(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l);

    pf_t get_e_stP(cand_pos_t i, cand_pos_t j);

    pf_t get_e_intP(cand_pos_t i, cand_pos_t ip, cand_pos_t jp, cand_pos_t j);

	pf_t beta2(cand_pos_t i, cand_pos_t l);
	pf_t beta2P(cand_pos_t i, cand_pos_t l);
	pf_t gamma2(cand_pos_t i, cand_pos_t l);
};


#endif