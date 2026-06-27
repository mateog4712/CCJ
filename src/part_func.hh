#ifndef PART_FUNC
#define PART_FUNC
#include "base_types.hh"
#include "matrices.hh"
#include <cstring>
#include <string>
#include <unordered_map>
#include <vector>

#include "ViennaRNA/loops.hh"
#include "ViennaRNA/pair_mat.hh"
#include "ViennaRNA/params/io.hh"

#ifdef NDEBUG
	#define UNREACHABLE() \
		do { \
			std::cerr << "Reached unreachable at line " << __LINE__ << " in File: " << __FILE__ << std::endl; \
			abort(); \
		} while(0)
#else
    #define UNREACHABLE() __builtin_unreachable()
#endif

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
    std::string MEA_structure;
    std::string centroid_structure;
    int num_samples;
    pf_t frequency;
    pf_t ensemble_diversity;
    std::unordered_map<std::string, int> structures;
    double gamma;

    W_final_pf(std::string &seq, std::string &MFE_structure, double MFE_energy, int dangle, int num_samples, bool print_samples, bool PSplot);
    // constructor for the restricted mfe case

    ~W_final_pf();
    // The destructor

    pf_t ccj_pf();

	pf_t ccj_centroid();

	TriangleMatrixPF P; // the main loop for pseudoloops and bands

  private:
    std::string seq;
    std::string MFE_structure;
	double MFE;
	cand_pos_t n;
	vrna_exp_param_t *exp_params_;
	bool print_samples;
    bool PSplot;
    std::vector<cand_pos_t> index;				// the array to keep the index of two dimensional arrays like WPP and WBP
	index_offset_t index3D;

    short *S_;
    short *S1_;

	// PK Free
    TriangleMatrixPF V;
    TriangleMatrixPF VM;
    TriangleMatrixPF WMv;
    TriangleMatrixPF WMp;
    TriangleMatrixPF WM;
    std::vector<pf_t> W;

    TriangleMatrixPF WPP;	// similar to WP but has at least one base pair
	TriangleMatrixPF WBP;	// similar to WB but has at least one base pair
	
	Matrix4DPF PL;		// MFE of a TGB structure s.t. i.j is paired
	Matrix4DPF PLmloop00;
	Matrix4DPF PLmloop01;
	Matrix4DPF PLmloop10;
	Matrix4DPF PfromL;
	Matrix4DPF PfromLprime;

	Matrix4DPF PR;		// MFE of a TGB structure s.t. k.l is paired
	Matrix4DPF PRmloop00;
	Matrix4DPF PRmloop01;
	Matrix4DPF PRmloop10;
	Matrix4DPF PfromR;
	Matrix4DPF PfromRprime;

	Matrix4DPF PM;		// MFE of a TGB structure s.t. j.k is paired
	Matrix4DPF PfromM;
	Matrix4DPF PfromMprime;
	Matrix4DPF PMmloop00;
	Matrix4DPF PMmloop01;
	Matrix4DPF PMmloop10;

	Matrix4DPF POm;	// MFE of a TGB structure s.t. i.l is paired multiple
	Matrix4DPF POmmloop00;
	Matrix4DPF POmmloop01;
	Matrix4DPF POmmloop10;
	Matrix4DPF PfromO;
	Matrix4DPF PfromOprime;
	
	Matrix4DPF POs;	// MFE of a TGB structure s.t. i.l is paired single
	Matrix4DPF POsmloop00;
	Matrix4DPF POsmloop01;
	Matrix4DPF POsmloop10;

	Matrix4DPF PLreO;
	Matrix4DPF PLreOmloop00;
	Matrix4DPF PLreOmloop01;
	Matrix4DPF PLreOmloop10;
	Matrix4DPF PfromLreO;
	Matrix4DPF PfromLreOprime;

	Matrix4DPF PLreR;
	Matrix4DPF PLreRmloop00;
	Matrix4DPF PLreRmloop01;
	Matrix4DPF PLreRmloop10;
	Matrix4DPF PfromLreR;
	Matrix4DPF PfromLreRprime;

	Matrix4DPF PMreO;
	Matrix4DPF PMreOmloop00;
	Matrix4DPF PMreOmloop01;
	Matrix4DPF PMreOmloop10;
	Matrix4DPF PfromMreO;
	Matrix4DPF PfromMreOprime;

	Matrix4DPF PMreR;
	Matrix4DPF PMreRmloop00;
	Matrix4DPF PMreRmloop01;
	Matrix4DPF PMreRmloop10;
	Matrix4DPF PfromMreR;
	Matrix4DPF PfromMreRprime;

	Matrix4DPF PLMreR;
	Matrix4DPF PLMreRmloop00;
	Matrix4DPF PLMreRmloop01;
	Matrix4DPF PLMreRmloop10;
	Matrix4DPF PfromLMreR;
	Matrix4DPF PfromLMreRprime;

	Matrix4DPF PLMorO;
	Matrix4DPF PLMorOmloop00;
	Matrix4DPF PLMorOmloop01;
	Matrix4DPF PLMorOmloop10;
	Matrix4DPF PfromLMorO;
	Matrix4DPF PfromLMorOprime;

	Matrix4DPF PK1Om;		// MFE of a TGB structure over gapped region [i,j] U [k,l]
	Matrix4DPF PK2Om;

	Matrix4DPF PK1Os;
	Matrix4DPF PK2Os;

	Matrix4DPF PK1LreO;
	Matrix4DPF PK2LreO;

	Matrix4DPF PK1MreO;
	Matrix4DPF PK2MreO;

	Matrix4DPF PK1LreR;
	Matrix4DPF PK2LreR;

	Matrix4DPF PK1MreR;
	Matrix4DPF PK2MreR;

	Matrix4DPF PK1R;
	Matrix4DPF PK2R;

	Matrix4DPF PK1LMreR;
	Matrix4DPF PK1LMorO;

    // Extra
    std::vector<pf_t> scale;
    std::vector<pf_t> expMLbase;
    std::vector<pf_t> expcp_pen;
    std::vector<pf_t> expPUP_pen;
    std::unordered_map<std::pair<cand_pos_t, cand_pos_t>, cand_pos_t, SzudzikHash> samples;

	///////// Functions //////

    inline pf_t to_Energy(pf_t energy, cand_pos_t length) {
    return ((-log(energy) - length * log(exp_params_->pf_scale)) * exp_params_->kT / 1000.0);
	}
	inline pf_t to_PF(pf_t energy, cand_pos_t length) {
		return exp(-(energy*1000/exp_params_->kT)-length*log(exp_params_->pf_scale));
	}
    void rescale_pk_globals();

    void exp_params_rescale();

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

	void compute_WBP(cand_pos_t i, cand_pos_t l);
	void compute_WPP(cand_pos_t i, cand_pos_t l);
		
	void compute_P(cand_pos_t i, cand_pos_t l);
	void compute_PK1X(const Index4D &x, MType type);
	void compute_PK2X(const Index4D &x, MType type);
	void compute_PX(const Index4D &x, MType type);
	void compute_PfromX(const Index4D &x, MType type);
	void compute_PfromXprime(const Index4D &x, MType type);
	pf_t calc_PfromXdoubleprime(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l, MType type);
	void compute_PXmloop00(const Index4D &x, MType type);
	void compute_PXmloop01(const Index4D &x, MType type);
	void compute_PXmloop10(const Index4D &x, MType type);

	pf_t calc_PLiloop(const Index4D &x, MType type);
	pf_t calc_PRiloop(const Index4D &x, MType type);
	pf_t calc_PMiloop(const Index4D &x, MType type);
	pf_t calc_POiloop(const Index4D &x, MType type);

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

	// inline pf_t calc_PO(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l){return POm.get(i,j,k,l)*POs.get(i,j,k,l);};
	pf_t calc_WB(cand_pos_t i, cand_pos_t l);
	pf_t calc_WP(cand_pos_t i, cand_pos_t l);
	pf_t calc_PfromLdoubleprime(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	pf_t calc_PfromRdoubleprime(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	pf_t calc_PfromMdoubleprime(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	pf_t calc_PfromOdoubleprime(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	pf_t calc_PfromLreOdoubleprime(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	pf_t calc_PfromLreRdoubleprime(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	pf_t calc_PfromMreOdoubleprime(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	pf_t calc_PfromMreRdoubleprime(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	pf_t calc_PfromLMreRdoubleprime(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);
	pf_t calc_PfromLMorOdoubleprime(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l);

	// penalty for closing pair i.l or l.i of a pseudoloop
	static constexpr pf_t gamma2(cand_pos_t i, cand_pos_t l){
		return 1.0;
	}

	pf_t calc_PXiloop(const Index4D &x, MType type);
	pf_t calc_PXmloop(const Index4D &x, MType type);
	inline bool impossible_case(const Index4D &x) const {
    	return !x.is_valid(n);
	}
	template<class Penalty> pf_t penalty(const Index4D &x, Penalty p, MType type) {
		switch(type) {
		case MType::L: return p(x.j(),x.i());
		case MType::M: return p(x.j(),x.k());
		case MType::R: return p(x.l(),x.k());
		case MType::Om: return p(x.l(),x.i());
		case MType::Os: return p(x.l(),x.i());
		case MType::LreO: return p(x.j(),x.i());
		case MType::LreR: return p(x.j(),x.i());
		case MType::MreO: return p(x.j(),x.k());
		case MType::MreR: return p(x.j(),x.k());
		case MType::LMreR: return p(x.j(),x.i());
		case MType::LMorO: return p(x.j(),x.i());
		}
		__builtin_unreachable();
	}
	inline Matrix4DPF& PX_by_mtype(MType type) {
		static std::array<Matrix4DPF*,11> matrices{&PL, &PM, &PR, &POm, &POs, &PLreO, &PLreR, &PMreO, &PMreR, &PLMreR, &PLMorO};
		return *matrices[static_cast<int>(type)];
	}
	inline Matrix4DPF& PfromX_by_mtype(MType type) {
		static std::array<Matrix4DPF*,11> matrices{&PfromL, &PfromM, &PfromR, &PfromO, &PfromO, &PfromLreO, &PfromLreR, &PfromMreO, &PfromMreR, &PfromLMreR, &PfromLMorO};
		return *matrices[static_cast<int>(type)];
	}
	inline Matrix4DPF& PfromXprime_by_mtype(MType type) {
		static std::array<Matrix4DPF*,11> matrices{&PfromLprime, &PfromMprime, &PfromRprime, &PfromOprime, &PfromOprime, &PfromLreOprime, &PfromLreRprime, &PfromMreOprime, &PfromMreRprime, &PfromLMreRprime, &PfromLMorOprime};
		return *matrices[static_cast<int>(type)];
	}
	inline Matrix4DPF& PXmloop00_by_mtype(MType type) {
		static std::array<Matrix4DPF*,11> matrices{&PLmloop00, &PMmloop00, &PRmloop00, &POmmloop00, &POsmloop00, &PLreOmloop00, &PLreRmloop00, &PMreOmloop00, &PMreRmloop00, &PLMreRmloop00, &PLMorOmloop00};
		return *matrices[static_cast<int>(type)];
	}
	inline Matrix4DPF& PXmloop01_by_mtype(MType type) {
		static std::array<Matrix4DPF*,11> matrices{&PLmloop01, &PMmloop01, &PRmloop01, &POmmloop01, &POsmloop01, &PLreOmloop01, &PLreRmloop01, &PMreOmloop01, &PMreRmloop01, &PLMreRmloop01, &PLMorOmloop01};
		return *matrices[static_cast<int>(type)];
	}
	inline Matrix4DPF& PXmloop10_by_mtype(MType type) {
		static std::array<Matrix4DPF*,11> matrices{&PLmloop10, &PMmloop10, &PRmloop10, &POmmloop10, &POsmloop10, &PLreOmloop10, &PLreRmloop10, &PMreOmloop10, &PMreRmloop10, &PLMreRmloop10, &PLMorOmloop10};
		return *matrices[static_cast<int>(type)];
	}
	/**
	 * Because of the ordering of Mtype, L and M (the first two are the only two missing). This means that I can just pass Mtype-2 and get the matrices I want
	 */
	inline Matrix4DPF& PK1X_by_mtype(MType type) {
		static std::array<Matrix4DPF*,9> matrices{&PK1R, &PK1Om, &PK1Os, &PK1LreO, &PK1LreR, &PK1MreO, &PK1MreR, &PK1LMreR, &PK1LMorO};
		return *matrices[static_cast<int>(type)-2];
	}
	/**
	 * This follows the same logic as PK1 but it is missing the last two
	 */
	inline Matrix4DPF& PK2X_by_mtype(MType type) {
		static std::array<Matrix4DPF*,7> matrices{&PK2R, &PK2Om, &PK2Os, &PK2LreO, &PK2LreR, &PK2MreO, &PK2MreR};
		return *matrices[static_cast<int>(type)-2];
	}


	// Stochastic Backtracking/ Sampling
	void Sample_W(cand_pos_t start, cand_pos_t end, std::vector<int> &fres);
	void Sample_V(cand_pos_t i, cand_pos_t j, std::vector<int> &fres);
	void Sample_VM(cand_pos_t i, cand_pos_t j, std::vector<int> &fres);
	void Sample_WM(cand_pos_t i, cand_pos_t j, std::vector<int> &fres);
	void Sample_WMv(cand_pos_t i, cand_pos_t j, std::vector<int> &fres);
	void Sample_WMp(cand_pos_t i, cand_pos_t j, std::vector<int> &fres);
	// Util
	void Sample_WP(cand_pos_t i, cand_pos_t l, std::vector<int> &fres);
	void Sample_WPP(cand_pos_t i, cand_pos_t l, std::vector<int> &fres);
	void Sample_WB(cand_pos_t i, cand_pos_t l, std::vector<int> &fres);
	void Sample_WBP(cand_pos_t i, cand_pos_t l, std::vector<int> &fres);
	// Main
	void Sample_P(cand_pos_t i, cand_pos_t l, std::vector<int> &fres);
	void Sample_PK1X(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l, MType type, std::vector<int> &fres);
	void Sample_PK2X(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l, MType type, std::vector<int> &fres);
	void Sample_PX(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l, MType type, std::vector<int> &fres);
	void Sample_PXmloop(const Index4D &x, MType type,std::vector<int> &fres);
	void Sample_PXiloop(const Index4D &x, MType type,std::vector<int> &fres);
	void Sample_PfromX(const Index4D &x, MType type,std::vector<int> &fres);
	void Sample_PfromXprime(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l, MType type,std::vector<int> &fres);
	void Sample_PfromXdoubleprime(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l, MType type,std::vector<int> &fres);
	void Sample_PXmloop00(const Index4D &x, MType type,std::vector<int> &fres);
	void Sample_PXmloop10(const Index4D &x, MType type,std::vector<int> &fres);
	void Sample_PXmloop01(const Index4D &x, MType type,std::vector<int> &fres);

	void Sample_PLiloop(const Index4D &x, MType type,std::vector<int> &fres);
	void Sample_PMiloop(const Index4D &x, MType type,std::vector<int> &fres);
	void Sample_PRiloop(const Index4D &x, MType type,std::vector<int> &fres);
	void Sample_POiloop(const Index4D &x, MType type,std::vector<int> &fres);

	void Sample_PfromL(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l, MType type,std::vector<int> &fres);
	void Sample_PfromM(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l, MType type,std::vector<int> &fres);
	void Sample_PfromR(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l, MType type,std::vector<int> &fres);
	void Sample_PfromO(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l, MType type,std::vector<int> &fres);

	void Sample_PfromLprime(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l, MType type,std::vector<int> &fres);
	void Sample_PfromMprime(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l, MType type,std::vector<int> &fres);
	void Sample_PfromRprime(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l, MType type,std::vector<int> &fres);
	void Sample_PfromOprime(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l, MType type,std::vector<int> &fres);

	void Sample_PfromLdoubleprime(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l,std::vector<int> &fres);
	void Sample_PfromMdoubleprime(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l,std::vector<int> &fres);
	void Sample_PfromRdoubleprime(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l,std::vector<int> &fres);
	void Sample_PfromOdoubleprime(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l,std::vector<int> &fres);
	void Sample_PfromLreOdoubleprime(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l,std::vector<int> &fres);
	void Sample_PfromLreRdoubleprime(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l,std::vector<int> &fres);
	void Sample_PfromMreOdoubleprime(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l,std::vector<int> &fres);
	void Sample_PfromMreRdoubleprime(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l,std::vector<int> &fres);
	void Sample_PfromLMreRdoubleprime(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l,std::vector<int> &fres);
	void Sample_PfromLMorOdoubleprime(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l,std::vector<int> &fres);

	void Sample_PLmloop00(const Index4D &x, MType type,std::vector<int> &fres);
	void Sample_PMmloop00(const Index4D &x, MType type,std::vector<int> &fres);
	void Sample_PRmloop00(const Index4D &x, MType type,std::vector<int> &fres);
	void Sample_POmloop00(const Index4D &x, MType type,std::vector<int> &fres);

	void Sample_PLmloop10(const Index4D &x, MType type,std::vector<int> &fres);
	void Sample_PMmloop10(const Index4D &x, MType type,std::vector<int> &fres);
	void Sample_PRmloop10(const Index4D &x, MType type,std::vector<int> &fres);
	void Sample_POmloop10(const Index4D &x, MType type,std::vector<int> &fres);

	void Sample_PLmloop01(const Index4D &x, MType type,std::vector<int> &fres);
	void Sample_PMmloop01(const Index4D &x, MType type,std::vector<int> &fres);
	void Sample_PRmloop01(const Index4D &x, MType type,std::vector<int> &fres);
	void Sample_POmloop01(const Index4D &x, MType type,std::vector<int> &fres);


	// Structure Stuff
	std::string compute_centroid(pf_t &dist, pf_t &diversity);
	char bpp_symbol(pf_t *P);
	void pairing_tendency();
};


#endif