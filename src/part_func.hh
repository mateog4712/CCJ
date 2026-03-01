#ifndef PART_FUNC
#define PART_FUNC
#include "base_types.hh"
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
        if (i >= j) return 0;
        cand_pos_t ij = index[i] + j - i;
        return V[ij];
    }
    pf_t get_energy_VM(cand_pos_t i, cand_pos_t j) {
        if (i >= j) return 0;
        cand_pos_t ij = index[i] + j - i;
        return VM[ij];
    }
    pf_t get_energy_WM(cand_pos_t i, cand_pos_t j) {
        if (i >= j) return 0;
        cand_pos_t ij = index[i] + j - i;
        return WM[ij];
    }
    pf_t get_energy_WMv(cand_pos_t i, cand_pos_t j) {
        if (i >= j) return 0;
        cand_pos_t ij = index[i] + j - i;
        return WMv[ij];
    }
    pf_t get_energy_WMp(cand_pos_t i, cand_pos_t j) {
        if (i >= j) return 0;
        cand_pos_t ij = index[i] + j - i;
        return WMp[ij];
    }
	pf_t get_P(cand_pos_t i, cand_pos_t j);
	pf_t get_PK(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	pf_t get_PL(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	pf_t get_PR(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	pf_t get_PM(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	pf_t get_PO(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	
	
	pf_t get_PfromL(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
    pf_t get_PfromR(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	pf_t get_PfromM(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	pf_t get_PfromO(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	
	
	pf_t get_PLiloop(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	pf_t get_PLiloop5(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l, cand_pos_t s);
	pf_t get_PLmloop(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	pf_t get_PLmloop00(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	pf_t get_PLmloop01(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	pf_t get_PLmloop10(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	
	pf_t get_PRiloop(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	pf_t get_PRiloop5(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l, cand_pos_t s);
	pf_t get_PRmloop(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	pf_t get_PRmloop00(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	pf_t get_PRmloop01(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	pf_t get_PRmloop10(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	
	pf_t get_PMiloop(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	pf_t get_PMiloop5(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l, cand_pos_t s);
	pf_t get_PMmloop(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	pf_t get_PMmloop00(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	pf_t get_PMmloop01(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	pf_t get_PMmloop10(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	
	pf_t get_POiloop(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	pf_t get_POiloop5(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l, cand_pos_t s);
	pf_t get_POmloop(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	pf_t get_POmloop00(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	pf_t get_POmloop01(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	pf_t get_POmloop10(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);

	std::vector<pf_t> P; // the main loop for pseudoloops and bands

  private:
    std::string seq;
    std::string MFE_structure;
    bool PSplot;
    cand_pos_t n;
    std::vector<cand_pos_t> index;

    short *S_;
    short *S1_;

    // PK Free
    std::vector<pf_t> V;
    std::vector<pf_t> VM;
    std::vector<pf_t> WMv;
    std::vector<pf_t> WMp;
    std::vector<pf_t> WM;
    std::vector<pf_t> W;
    // PK
    std::vector<pf_t> WPP;								// similar to WP but has at least one base pair
	std::vector<pf_t> WBP;								// similar to WB but has at least one base pair
	
	std::vector<std::vector<pf_t> > PK;					// MFE of a TGB structure over gapped region [i,j] U [k,l]
	std::vector<std::vector< std::vector<pf_t> > > PL;	// MFE of a TGB structure s.t. i.j is paired
	std::vector<std::vector<pf_t> > PR;					// MFE of a TGB structure s.t. k.l is paired
	std::vector<std::vector<pf_t> > PM;					// MFE of a TGB structure s.t. j.k is paired
	std::vector<std::vector< std::vector<pf_t> > > PO;	// MFE of a TGB structure s.t. i.l is paired
	
	// transition recurrences
	std::vector<std::vector< std::vector<pf_t> > > PfromL;
	std::vector<std::vector<pf_t> > PfromR;
	std::vector<std::vector<pf_t> > PfromM;
	std::vector<std::vector< std::vector<pf_t> > > PfromO;
	
	// internal loops and multi loops that span a band
	std::vector<std::vector<pf_t> > PLiloop;
	std::vector<std::vector<pf_t> > PLmloop00;
	std::vector<std::vector< std::vector<pf_t> > > PLmloop01;
	std::vector<std::vector< std::vector<pf_t> > > PLmloop10;
	
	
	std::vector<std::vector<pf_t> > PRiloop;
	std::vector<std::vector<pf_t> > PRmloop00;
	std::vector<std::vector<pf_t> > PRmloop01;
	std::vector<std::vector<pf_t> > PRmloop10;
	
	
	std::vector<std::vector<pf_t> > PMiloop;
	std::vector<std::vector<pf_t> > PMmloop00;
	std::vector<std::vector<pf_t> > PMmloop01;
	std::vector<std::vector<pf_t> > PMmloop10;
	
	
	std::vector<std::vector<pf_t> > POiloop;
	std::vector<std::vector<pf_t> > POmloop00;
	std::vector<std::vector< std::vector<pf_t> > > POmloop01;
	std::vector<std::vector< std::vector<pf_t> > > POmloop10;

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
	
	
	void compute_PLiloop(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	void compute_PLmloop00(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	void compute_PLmloop01(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	void compute_PLmloop10(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	
	void compute_PRiloop(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	void compute_PRmloop00(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	void compute_PRmloop01(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	void compute_PRmloop10(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	
	void compute_PMiloop(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	void compute_PMmloop00(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	void compute_PMmloop01(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	void compute_PMmloop10(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	
	void compute_POiloop(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
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