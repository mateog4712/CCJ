#include "part_func.hh"
#include "dot_plot.hh"
#include "h_globals.hh"
#include "pf_globals.hh"

#include <algorithm>
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <string>
#include <cassert>
#include <list>

#define debug 0
/*
 * If the global use_mfelike_energies flag is set, truncate doubles to int
 * values and cast back to double. This makes the energy parameters of the
 * partition (folding get_scaled_exp_params()) compatible with the mfe folding
 * parameters (get_scaled_exp_params()), e.g. for explicit partition function
 * computations.
 */
#define TRUNC_MAYBE(X) ((!pf_smooth) ? (double)((int)(X)) : (X))
/* Rescale Free energy contribution according to deviation of temperature from measurement conditions */
#define RESCALE_dG(dG, dH, dT) ((dH) - ((dH) - (dG)) * dT)

/*
 * Rescale Free energy contribution according to deviation of temperature from measurement conditions
 * and convert it to Boltzmann Factor for specific kT
 */
#define RESCALE_BF(dG, dH, dT, kT) (exp(-TRUNC_MAYBE((double)RESCALE_dG((dG), (dH), (dT))) * 10. / kT))

W_final_pf::W_final_pf(std::string &seq, std::string &MFE_structure, double MFE_energy, int dangle, int num_samples, bool print_samples, bool PSplot)
    : seq(seq), MFE_structure(MFE_structure), MFE(MFE_energy), exp_params_(vrna_exp_params(NULL)), print_samples(print_samples), PSplot(PSplot) {
    this->n = seq.length();
    this->num_samples = num_samples;

    make_pair_matrix();
    exp_params_->model_details.dangles = dangle;
    S_ = encode_sequence(seq.c_str(), 0);
    S1_ = encode_sequence(seq.c_str(), 1);

    index.resize(n + 1);
    scale.resize(n + 1);
    expMLbase.resize(n + 1);
    expcp_pen.resize(n + 1);
    expPUP_pen.resize(n + 1);
	TriangleMatrixPF::new_index(index,n+1);
	Matrix4DPF::construct_index(index3D,n);
	init_can_pair();
    // Allocate space
    V.init(n+1,index);
    VM.init(n+1,index);
    WM.init(n+1,index);
    WMv.init(n+1,index);
    WMp.init(n+1,index);

    // PK
    WBP.init(n+1,index);
	WPP.init(n+1,index);
	P.init(n+1,index);

	// 4D matrix initialization
	PL.init(n,index3D);
	PfromL.init(n,index3D);
	PfromLprime.init(n,index3D);
	PLmloop00.init(n,index3D);
	PLmloop01.init(n,index3D);
	PLmloop10.init(n,index3D);

	PR.init(n,index3D);
	PfromR.init(n,index3D);
	PfromRprime.init(n,index3D);
	PRmloop00.init(n,index3D);
	PRmloop01.init(n,index3D);
	PRmloop10.init(n,index3D);

	PM.init(n,index3D);
	PfromM.init(n,index3D);
	PfromMprime.init(n,index3D);
	PMmloop00.init(n,index3D);
	PMmloop01.init(n,index3D);
	PMmloop10.init(n,index3D);

	POm.init(n,index3D);
	PfromO.init(n,index3D);
	PfromOprime.init(n,index3D);
	POmmloop00.init(n,index3D);
	POmmloop01.init(n,index3D);
	POmmloop10.init(n,index3D);

	POs.init(n,index3D);
	POsmloop00.init(n,index3D);
	POsmloop01.init(n,index3D);
	POsmloop10.init(n,index3D);

	PLreO.init(n,index3D);
	PfromLreO.init(n,index3D);
	PfromLreOprime.init(n,index3D);
	PLreOmloop00.init(n,index3D);
	PLreOmloop01.init(n,index3D);
	PLreOmloop10.init(n,index3D);

	PLreR.init(n,index3D);
	PfromLreR.init(n,index3D);
	PfromLreRprime.init(n,index3D);
	PLreRmloop00.init(n,index3D);
	PLreRmloop01.init(n,index3D);
	PLreRmloop10.init(n,index3D);

	PMreO.init(n,index3D);
	PfromMreO.init(n,index3D);
	PfromMreOprime.init(n,index3D);
	PMreOmloop00.init(n,index3D);
	PMreOmloop01.init(n,index3D);
	PMreOmloop10.init(n,index3D);

	PMreR.init(n,index3D);
	PfromMreR.init(n,index3D);
	PfromMreRprime.init(n,index3D);
	PMreRmloop00.init(n,index3D);
	PMreRmloop01.init(n,index3D);
	PMreRmloop10.init(n,index3D);

	PLMreR.init(n,index3D);
	PfromLMreR.init(n,index3D);
	PfromLMreRprime.init(n,index3D);
	PLMreRmloop00.init(n,index3D);
	PLMreRmloop01.init(n,index3D);
	PLMreRmloop10.init(n,index3D);

	PLMorO.init(n,index3D);
	PfromLMorO.init(n,index3D);
	PfromLMorOprime.init(n,index3D);
	PLMorOmloop00.init(n,index3D);
	PLMorOmloop01.init(n,index3D);
	PLMorOmloop10.init(n,index3D);

	PK1Om.init(n,index3D);
	PK2Om.init(n,index3D);

	PK1Os.init(n,index3D);
	PK2Os.init(n,index3D);

	PK1LreO.init(n,index3D);
	PK2LreO.init(n,index3D);

	PK1MreO.init(n,index3D);
	PK2MreO.init(n,index3D);

	PK1LreR.init(n,index3D);
	PK2LreR.init(n,index3D);

	PK1MreR.init(n,index3D);
	PK2MreR.init(n,index3D);

	PK1R.init(n,index3D);
	PK2R.init(n,index3D);

	PK1LMreR.init(n,index3D);
	PK1LMorO.init(n,index3D);

    rescale_pk_globals();
    exp_params_rescale();
    W.resize(n+1,scale[1]);
}

W_final_pf::~W_final_pf() {}

void W_final_pf::exp_params_rescale() {
    double e_per_nt, kT;
    kT = exp_params_->kT;

    e_per_nt = MFE * 1000. / this->n;

    exp_params_->pf_scale = exp(-(exp_params_->model_details.sfact * e_per_nt) / kT);

    if (exp_params_->pf_scale < 1.) exp_params_->pf_scale = 1.;

    exp_params_->pf_scale = 1.; // I don't think it even needs scaling for the lengths it can do

    this->scale[0] = 1.;
    this->scale[1] = (pf_t)(1. / exp_params_->pf_scale);
    this->expMLbase[0] = 1;
    this->expMLbase[1] = (pf_t)(exp_params_->expMLbase / exp_params_->pf_scale);

    this->expcp_pen[0] = 1;
    this->expcp_pen[1] = (pf_t)(expcp_penalty / exp_params_->pf_scale);
    this->expPUP_pen[0] = 1;
    this->expPUP_pen[1] = (pf_t)(expPUP_penalty / exp_params_->pf_scale);

    for (cand_pos_t i = 2; i <= this->n; i++) {
        this->scale[i] = this->scale[i / 2] * this->scale[i - (i / 2)];
        this->expMLbase[i] = (pf_t)pow(exp_params_->expMLbase, (double)i) * this->scale[i];
        this->expcp_pen[i] = (pf_t)pow(expcp_penalty, (double)i) * this->scale[i];
        this->expPUP_pen[i] = (pf_t)pow(expPUP_penalty, (double)i) * this->scale[i];
    }
}

void W_final_pf::rescale_pk_globals() {
    double kT = exp_params_->model_details.betaScale * (exp_params_->model_details.temperature + K0) * GASCONST; /* kT in cal/mol  */
    double TT = (exp_params_->model_details.temperature + K0) / (Tmeasure);
    int pf_smooth = exp_params_->model_details.pf_smooth;

    expPS_penalty = RESCALE_BF(PS_penalty, PS_penalty, TT, kT);
    expPSM_penalty = RESCALE_BF(PSM_penalty, PSM_penalty, TT, kT);
    expPSP_penalty = RESCALE_BF(PSP_penalty, PSP_penalty, TT, kT);
    expPB_penalty = RESCALE_BF(PB_penalty, PB_penalty, TT, kT);
    expPUP_penalty = RESCALE_BF(PUP_penalty, PUP_penalty, TT, kT);
    expPPS_penalty = RESCALE_BF(PPS_penalty, PPS_penalty, TT, kT);

    expa_penalty = RESCALE_BF(a_penalty, ML_closingdH, TT, kT); //Should I do a,b,c here?
    expb_penalty = RESCALE_BF(b_penalty, ML_interndH, TT, kT);
    expc_penalty = RESCALE_BF(c_penalty, ML_BASEdH, TT, kT);

    expap_penalty = RESCALE_BF(ap_penalty, ap_penalty, TT, kT);
    expbp_penalty = RESCALE_BF(bp_penalty, bp_penalty, TT, kT);
    expcp_penalty = RESCALE_BF(cp_penalty, cp_penalty, TT, kT);
}
pf_t W_final_pf::ccj_pf(){

	for (cand_pos_t i = n; i>=1; --i){	
		for (cand_pos_t j =i; j<=n; ++j){
			compute_energy(i, j);
			compute_pk_energies(i,j);
			compute_WMv_WMp(i,j);
			compute_energy_WM(i,j);
		}
	}

	for (cand_pos_t j= TURN+1; j <= n; j++){
		pf_t contributions = 0;
		contributions += W[j-1]*scale[1];
		for (cand_pos_t k=1; k<=j-TURN-1; ++k){
			pf_t acc = (k>1) ? W[k-1]: 1;
			contributions += acc*V.get(k,j)*exp_Extloop(k, j);
			contributions += acc*P.get(k,j)*expPS_penalty;
		}
		W[j] = contributions;
	}
    pf_t energy = to_Energy(W[n], n);

    structure = std::string(n, '.');
	for (cand_pos_t i = 0; i < num_samples; ++i) {
		std::string sample_structure(n,'.');
		std::vector<int> fres(n,-2);
        Sample_W(1, n, fres);
		fill_structure(fres,sample_structure);
        structures[sample_structure]++;
    }

	if (print_samples) {
        std::vector<std::pair<std::string,int>> str_list;
        for (const auto &s : structures) {
            str_list.push_back(std::make_pair(s.first,s.second));
        }
        sort(str_list.begin(), str_list.end(),[](auto &x,auto &y) {return x.second>y.second;} );
        for (const auto &s : str_list) {
            std::cout << s.first << " " << s.second << std::endl;
        }
    }
	pairing_tendency();
	this->frequency = (pf_t)structures[MFE_structure] / num_samples;

	if (PSplot) {
        create_dot_plot(seq, MFE_structure, samples, num_samples);
    }

    return energy;
}

pf_t W_final_pf::ccj_centroid(){
    pf_t dist = 0;
    pf_t diversity = 0;
    std::string centroid = compute_centroid(dist,diversity);
    this->centroid_structure = centroid;
    this->ensemble_diversity = diversity;
    return dist;
}

pf_t W_final_pf::exp_Extloop(cand_pos_t i, cand_pos_t j) {
    pair_type tt = pair[S_[i]][S_[j]];

    if (exp_params_->model_details.dangles == 1 || exp_params_->model_details.dangles == 2) {
        base_type si1 = i > 1 ? S_[i-1] : -1;
        base_type sj1 = j < n ? S_[j+1] : -1;
        return exp_E_ExtLoop(tt, si1, sj1, exp_params_);
    } else {
        return exp_E_ExtLoop(tt, -1, -1, exp_params_);
    }
}

pf_t W_final_pf::exp_MLstem(cand_pos_t i, cand_pos_t j) {
    pair_type tt = pair[S_[i]][S_[j]];
    if (exp_params_->model_details.dangles == 1 || exp_params_->model_details.dangles == 2) {
        base_type si1 = i > 1 ? S_[i-1] : -1;
        base_type sj1 = j < n ? S_[j+1] : -1;
        return exp_E_MLstem(tt, si1, sj1, exp_params_);
    } else {
        return exp_E_MLstem(tt, -1, -1, exp_params_);
    }
}

pf_t W_final_pf::exp_Mbloop(cand_pos_t i, cand_pos_t j) {
    pair_type tt = pair[S_[j]][S_[i]];
    if (exp_params_->model_details.dangles == 1 || exp_params_->model_details.dangles == 2) {
        base_type si1 = i > 1 ? S_[i+1] : -1;
        base_type sj1 = j < n ? S_[j-1] : -1;
        return exp_E_MLstem(tt, sj1, si1, exp_params_);
    } else {
        return exp_E_MLstem(tt, -1, -1, exp_params_);
    }
}

pf_t W_final_pf::HairpinE(cand_pos_t i, cand_pos_t j) {
    const int ptype_closing = pair[S_[i]][S_[j]];
    if (ptype_closing == 0) return 0;
    pf_t e_h = static_cast<pf_t>(exp_E_Hairpin(j - i - 1, ptype_closing, S1_[i + 1], S1_[j - 1], &seq.c_str()[i - 1], exp_params_));
    e_h *= scale[j-i+1];
    return e_h;
}

pf_t W_final_pf::compute_internal(cand_pos_t i, cand_pos_t j) {
    pf_t v_iloop = 0;
    cand_pos_t max_k = std::min(j - TURN - 2, i + MAXLOOP + 1);
    const pair_type ptype_closing = pair[S_[i]][S_[j]];
    for (cand_pos_t k = i + 1; k <= max_k; ++k) {
        cand_pos_t min_l = std::max(k + TURN + 1 + MAXLOOP + 2, k + j - i) - MAXLOOP - 2;
        for (cand_pos_t l = j - 1; l >= min_l; --l) {
            pf_t v_iloop_kl = V.get(k, l)
                                * exp_E_IntLoop(k - i - 1, j - l - 1, ptype_closing, rtype[pair[S_[k]][S_[l]]], S1_[i + 1], S1_[j - 1],
                                                S1_[k - 1], S1_[l + 1], exp_params_);
            cand_pos_t u1 = k - i - 1;
            cand_pos_t u2 = j - l - 1;
            v_iloop_kl *= scale[u1 + u2 + 2];
            v_iloop += v_iloop_kl;
        }
    }

    return v_iloop;
}

void W_final_pf::compute_WMv_WMp(cand_pos_t i, cand_pos_t j) {
    if (j - i - 1 < TURN) return;
    cand_pos_t ij = index[(i)] + (j) - (i);

    pf_t WMv_contributions = 0;
    pf_t WMp_contributions = 0;

    WMv_contributions += (V.get(i,j) * exp_MLstem(i, j));
    WMp_contributions += (P.get(i,j) * expPSM_penalty * expb_penalty);
    WMv_contributions += (WMv.get(i,j-1) * expMLbase[1]);
    WMp_contributions += (WMp.get(i,j-1) * expMLbase[1]);

    WMv[ij] = WMv_contributions;
    WMp[ij] = WMp_contributions;
}

void W_final_pf::compute_energy_WM(cand_pos_t i, cand_pos_t j) {
    if (j - i + 1 < 4) return;
    pf_t contributions = 0;
    cand_pos_t ij = index[(i)] + (j) - (i);
    cand_pos_t ijminus1 = index[(i)] + (j)-1 - (i);

    for (cand_pos_t k = i; k < j - TURN; ++k) {
        pf_t qbt1 = V.get(k,j) * exp_MLstem(k, j);
        pf_t qbt2 = P.get(k,j) * expPSM_penalty * expb_penalty;
        contributions += (static_cast<pf_t>(expMLbase[k - i]) * qbt1);
        contributions += (static_cast<pf_t>(expMLbase[k - i]) * qbt2);
        contributions += (WM.get(i,k-1) * qbt1);
        contributions += (WM.get(i,k-1) * qbt2);
    }
    contributions += WM[ijminus1] * expMLbase[1];
    WM[ij] = contributions;
}

pf_t W_final_pf::compute_energy_VM(cand_pos_t i, cand_pos_t j) {
    pf_t contributions = 0;
    cand_pos_t ij = index[(i)] + (j) - (i);
    for (cand_pos_t k = i + 1; k <= j - TURN - 1; ++k) {
        contributions += (WM.get(i+1,k-1) * WMv.get(k,j-1) * exp_Mbloop(i, j) * exp_params_->expMLclosing);
        contributions += (WM.get(i+1,k-1) * WMp.get(k,j-1) * exp_Mbloop(i, j) * exp_params_->expMLclosing);
        contributions += (expMLbase[k - i - 1] * WMp.get(k,j-1) * exp_Mbloop(i, j) * exp_params_->expMLclosing);
    }

    contributions *= scale[2];
    VM[ij] = contributions;
    return contributions;
}

void W_final_pf::compute_energy(cand_pos_t i, cand_pos_t j) {

    cand_pos_t ij = index[i] + j - i;

    pf_t contributions = 0;

    contributions += HairpinE(i, j);
    contributions += compute_internal(i, j);
    contributions += compute_energy_VM(i, j);
    V[ij] = contributions;
}

void W_final_pf::compute_pk_energies(cand_pos_t i, cand_pos_t l) {
	
    // 1) compute all energies over region [i,l]
	compute_P(i,l);
	compute_WBP(i,l);
	compute_WPP(i,l);

	//2) compute all energies over gapped region [i,j]U[k,l]
	for(cand_pos_t j = i; j<l; ++j){
		// Hosna, July 8, 2014
		// in original recurrences we have j< k-1, so I am changing k=j+1 to k=j+2
		for(cand_pos_t k = l; k>=j+2; --k){
			Index4D x(i,j,k,l);

			compute_PXmloop00(x,MType::L);
			compute_PXmloop00(x,MType::M);
			compute_PXmloop00(x,MType::R);
			compute_PXmloop00(x,MType::Om);
			compute_PXmloop00(x,MType::Os);
			compute_PXmloop00(x,MType::LreO);
			compute_PXmloop00(x,MType::LreR);
			compute_PXmloop00(x,MType::MreO);
			compute_PXmloop00(x,MType::MreR);
			compute_PXmloop00(x,MType::LMreR);
			compute_PXmloop00(x,MType::LMorO);

			compute_PXmloop01(x,MType::L);
			compute_PXmloop01(x,MType::M);
			compute_PXmloop01(x,MType::R);
			compute_PXmloop01(x,MType::Om);
			compute_PXmloop01(x,MType::Os);
			compute_PXmloop01(x,MType::LreO);
			compute_PXmloop01(x,MType::LreR);
			compute_PXmloop01(x,MType::MreO);
			compute_PXmloop01(x,MType::MreR);
			compute_PXmloop01(x,MType::LMreR);
			compute_PXmloop01(x,MType::LMorO);

			compute_PXmloop10(x,MType::L);
			compute_PXmloop10(x,MType::M);
			compute_PXmloop10(x,MType::R);
			compute_PXmloop10(x,MType::Om);
			compute_PXmloop10(x,MType::Os);
			compute_PXmloop10(x,MType::LreO);
			compute_PXmloop10(x,MType::LreR);
			compute_PXmloop10(x,MType::MreO);
			compute_PXmloop10(x,MType::MreR);
			compute_PXmloop10(x,MType::LMreR);
			compute_PXmloop10(x,MType::LMorO);

			compute_PX(x,MType::L);
			compute_PX(x,MType::M);
			compute_PX(x,MType::R);
			compute_PX(x,MType::Om);
			compute_PX(x,MType::Os);
			compute_PX(x,MType::LreO);
			compute_PX(x,MType::LreR);
			compute_PX(x,MType::MreO);
			compute_PX(x,MType::MreR);
			compute_PX(x,MType::LMreR);
			compute_PX(x,MType::LMorO);

			compute_PfromX(x,MType::L);
			compute_PfromX(x,MType::M);
			compute_PfromX(x,MType::R);
			compute_PfromX(x,MType::Om);
			compute_PfromX(x,MType::LreO);
			compute_PfromX(x,MType::LreR);
			compute_PfromX(x,MType::MreO);
			compute_PfromX(x,MType::MreR);
			compute_PfromX(x,MType::LMreR);
			compute_PfromX(x,MType::LMorO);

			compute_PfromXprime(x,MType::L);
			compute_PfromXprime(x,MType::R);
			compute_PfromXprime(x,MType::M);
			compute_PfromXprime(x,MType::Om);
			compute_PfromXprime(x,MType::LreO);
			compute_PfromXprime(x,MType::LreR);
			compute_PfromXprime(x,MType::MreO);
			compute_PfromXprime(x,MType::MreR);
			compute_PfromXprime(x,MType::LMreR);
			compute_PfromXprime(x,MType::LMorO);

			compute_PK1X(x,MType::Om);
			compute_PK1X(x,MType::Os);
			compute_PK1X(x,MType::LreO);
			compute_PK1X(x,MType::LreR);
			compute_PK1X(x,MType::MreO);
			compute_PK1X(x,MType::MreR);
			compute_PK1X(x,MType::R);
			compute_PK1X(x,MType::LMreR);
			compute_PK1X(x,MType::LMorO);

			compute_PK2X(x,MType::Om);
			compute_PK2X(x,MType::Os);
			compute_PK2X(x,MType::LreO);
			compute_PK2X(x,MType::LreR);
			compute_PK2X(x,MType::MreO);
			compute_PK2X(x,MType::MreR);
			compute_PK2X(x,MType::R);
		}
	}
}

void W_final_pf::compute_WBP(cand_pos_t i, cand_pos_t l){
	pf_t contributions = 0;

	for(cand_pos_t d=i; d< l; ++d){
		contributions += calc_WB(i,d-1)*V.get(d,l)*expbp_penalty*expPPS_penalty;
		contributions += calc_WB(i,d-1)*P.get(d,l)*expPSM_penalty*expPPS_penalty;
	}
	contributions+=WBP.get(i,l-1)*expcp_pen[1];
	WBP.set(i,l) = contributions;
}

void W_final_pf::compute_WPP(cand_pos_t i, cand_pos_t l){
	pf_t contributions = 0;

	for(cand_pos_t d=i; d<l; ++d){
		contributions += calc_WP(i,d-1)*V.get(d,l)*gamma2(l,d)*expPPS_penalty;
		contributions += calc_WP(i,d-1)*P.get(d,l)*expPSP_penalty*expPPS_penalty;
	}
	contributions+=WBP.get(i,l-1)*expPUP_pen[1];
	WPP.set(i,l) = contributions;
}

void W_final_pf::compute_P(cand_pos_t i, cand_pos_t l){
    pf_t contributions = 0;
	for(cand_pos_t j=i; j<l; ++j){
		for (cand_pos_t d=j+1; d<l; ++d){
			for (cand_pos_t k=d+1; k<l; ++k){
				contributions += PK1Om.get(i,j,d+1,k)*PK2Om.get(j+1,d,k+1,l);
				contributions += PK1Om.get(i,j,d+1,k)*PK2Os.get(j+1,d,k+1,l);
				contributions += PK1Om.get(i,j,d+1,k)*PK2LreO.get(j+1,d,k+1,l);
				contributions += PK1Om.get(i,j,d+1,k)*PK2MreO.get(j+1,d,k+1,l);
				contributions += PK1Om.get(i,j,d+1,k)*PK2LreR.get(j+1,d,k+1,l);
				contributions += PK1Om.get(i,j,d+1,k)*PK2MreR.get(j+1,d,k+1,l);
				contributions += PK1Om.get(i,j,d+1,k)*PK2R.get(j+1,d,k+1,l);

				contributions += PK1Os.get(i,j,d+1,k)*PK2Om.get(j+1,d,k+1,l);
				contributions += PK1Os.get(i,j,d+1,k)*PK2Os.get(j+1,d,k+1,l);
				contributions += PK1Os.get(i,j,d+1,k)*PK2MreO.get(j+1,d,k+1,l);
				contributions += PK1Os.get(i,j,d+1,k)*PK2MreR.get(j+1,d,k+1,l);
				contributions += PK1Os.get(i,j,d+1,k)*PK2R.get(j+1,d,k+1,l);

				contributions += PK1LreO.get(i,j,d+1,k)*PK2Om.get(j+1,d,k+1,l);
				contributions += PK1LreO.get(i,j,d+1,k)*PK2LreO.get(j+1,d,k+1,l);
				contributions += PK1LreO.get(i,j,d+1,k)*PK2MreO.get(j+1,d,k+1,l);

				contributions += PK1LMreR.get(i,j,d+1,k)*PK2Om.get(j+1,d,k+1,l);
				contributions += PK1LMreR.get(i,j,d+1,k)*PK2LreO.get(j+1,d,k+1,l);
				contributions += PK1LMreR.get(i,j,d+1,k)*PK2MreO.get(j+1,d,k+1,l);

				contributions += PK1LMorO.get(i,j,d+1,k)*PK2LreR.get(j+1,d,k+1,l);
				contributions += PK1LMorO.get(i,j,d+1,k)*PK2MreR.get(j+1,d,k+1,l);
				contributions += PK1LMorO.get(i,j,d+1,k)*PK2R.get(j+1,d,k+1,l);
			}
		}
	}
	P.set(i,l) = contributions;
}
/**
 * This always chops from the interior with every single recurrence which means the code is the same and can be made generic
 */
void W_final_pf::compute_PK1X(const Index4D &x, MType type){
	pf_t contributions = 0;
	const cand_pos_t i = x.i(), j = x.j(), k = x.k(), l = x.l();
	Matrix4DPF &PX = PX_by_mtype(type);

	for(cand_pos_t d = i+1;d<=j;++d){
		contributions += PX.get(i,d,k,l)*calc_WP(d+1,j)*expPB_penalty;
	}
	Matrix4DPF &PK1X = PK1X_by_mtype(type);
	PK1X.set(i,j,k,l,contributions);
}
/**
 * Same as PK1
 */
void W_final_pf::compute_PK2X(const Index4D &x, MType type){
	pf_t contributions = 0;
	const cand_pos_t i = x.i(), j = x.j(), k = x.k(), l = x.l();
	Matrix4DPF &PK1X = PK1X_by_mtype(type);

	for(cand_pos_t d = k;d<l;++d){
		contributions += PK1X.get(i,j,d,l)*calc_WP(k,d-1);
	}
	Matrix4DPF &PK2X = PK2X_by_mtype(type);
	PK2X.set(i,j,k,l,contributions);
}
/**
 * This was a case where the code was almost the same for all recurrences.
 * Each time, you would take the respective recurrence, shrink the related indices
 * by one and then add the penalties. This could be solve by the Index4D shrink function
 * which shrinks based on the Mtype. And so 10 functions become one
 */
pf_t W_final_pf::calc_PXmloop(const Index4D &x, MType type){
	if(impossible_case(x)) return 0;

	Index4D xp(x);
	xp.shrink(type);
	Matrix4DPF &PXmloop00 = PXmloop00_by_mtype(type);
	return PXmloop00.get(xp.i(),xp.j(),xp.k(),xp.l())*expap_penalty*expbp_penalty;
}
pf_t W_final_pf::calc_PXiloop(const Index4D &x, MType type){
	switch(type) {
    case MType::L: return calc_PLiloop(x,type);
    case MType::M: return calc_PMiloop(x,type);
    case MType::R: return calc_PRiloop(x,type);
    case MType::Om: return calc_POiloop(x,type);
	case MType::Os: return calc_POiloop(x,type);
	case MType::LreO: return calc_PLiloop(x,type);
	case MType::LreR: return calc_PLiloop(x,type);
	case MType::MreO: return calc_PMiloop(x,type);
	case MType::MreR: return calc_PMiloop(x,type);
	case MType::LMreR: return calc_PLiloop(x,type);
	case MType::LMorO: return calc_PLiloop(x,type);
    }
    UNREACHABLE();
}
void W_final_pf::compute_PfromX(const Index4D &x, MType type){
    switch(type) {
    case MType::L: return compute_PfromL(x,type);
    case MType::M: return compute_PfromM(x,type);
    case MType::R: return compute_PfromR(x,type);
    case MType::Om: return compute_PfromO(x,type);
	case MType::Os: return;
	case MType::LreO: return compute_PfromL(x,type);
	case MType::LreR: return compute_PfromL(x,type);
	case MType::MreO: return compute_PfromM(x,type);
	case MType::MreR: return compute_PfromM(x,type);
	case MType::LMreR: return compute_PfromL(x,type);
	case MType::LMorO: return compute_PfromL(x,type);
    }
    UNREACHABLE();
}
void W_final_pf::compute_PfromXprime(const Index4D &x, MType type){
    switch(type) {
    case MType::L: return compute_PfromLprime(x,type);
    case MType::M: return compute_PfromMprime(x,type);
    case MType::R: return compute_PfromRprime(x,type);
    case MType::Om: return compute_PfromOprime(x,type);
	case MType::Os: return;
	case MType::LreO: return compute_PfromLprime(x,type);
	case MType::LreR: return compute_PfromLprime(x,type);
	case MType::MreO: return compute_PfromMprime(x,type);
	case MType::MreR: return compute_PfromMprime(x,type);
	case MType::LMreR: return compute_PfromLprime(x,type);
	case MType::LMorO: return compute_PfromLprime(x,type);
    }
    UNREACHABLE();
}
pf_t W_final_pf::calc_PfromXdoubleprime(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l, MType type){
    switch(type) {
    case MType::L: return calc_PfromLdoubleprime(i,j,k,l);
    case MType::M: return calc_PfromMdoubleprime(i,j,k,l);
    case MType::R: return calc_PfromRdoubleprime(i,j,k,l);
    case MType::Om: return calc_PfromOdoubleprime(i,j,k,l);
	case MType::Os: return 0;
	case MType::LreO: return calc_PfromLreOdoubleprime(i,j,k,l);
	case MType::LreR: return calc_PfromLreRdoubleprime(i,j,k,l);
	case MType::MreO: return calc_PfromMreOdoubleprime(i,j,k,l);
	case MType::MreR: return calc_PfromMreRdoubleprime(i,j,k,l);
	case MType::LMreR: return calc_PfromLMreRdoubleprime(i,j,k,l);
	case MType::LMorO: return calc_PfromLMorOdoubleprime(i,j,k,l);
    }
    UNREACHABLE();
}
void W_final_pf::compute_PXmloop00(const Index4D &x, MType type){
    switch(type) {
    case MType::L: return compute_PLmloop00(x,type);
    case MType::M: return compute_PMmloop00(x,type);
    case MType::R: return compute_PRmloop00(x,type);
    case MType::Om: return compute_POmloop00(x,type);
	case MType::Os: return compute_POmloop00(x,type);
	case MType::LreO: return compute_PLmloop00(x,type);
	case MType::LreR: return compute_PLmloop00(x,type);
	case MType::MreO: return compute_PMmloop00(x,type);
	case MType::MreR: return compute_PMmloop00(x,type);
	case MType::LMreR: return compute_PLmloop00(x,type);
	case MType::LMorO: return compute_PLmloop00(x,type);
    }
    UNREACHABLE();
}
void W_final_pf::compute_PXmloop01(const Index4D &x, MType type){
    switch(type) {
    case MType::L: return compute_PLmloop01(x,type);
    case MType::M: return compute_PMmloop01(x,type);
    case MType::R: return compute_PRmloop01(x,type);
    case MType::Om: return compute_POmloop01(x,type);
	case MType::Os: return compute_POmloop01(x,type);
	case MType::LreO: return compute_PLmloop01(x,type);
	case MType::LreR: return compute_PLmloop01(x,type);
	case MType::MreO: return compute_PMmloop01(x,type);
	case MType::MreR: return compute_PMmloop01(x,type);
	case MType::LMreR: return compute_PLmloop01(x,type);
	case MType::LMorO: return compute_PLmloop01(x,type);
    }
    UNREACHABLE();
}
void W_final_pf::compute_PXmloop10(const Index4D &x, MType type){
    switch(type) {
    case MType::L: return compute_PLmloop10(x,type);
    case MType::M: return compute_PMmloop10(x,type);
    case MType::R: return compute_PRmloop10(x,type);
    case MType::Om: return compute_POmloop10(x,type);
	case MType::Os: return compute_POmloop10(x,type);
	case MType::LreO: return compute_PLmloop10(x,type);
	case MType::LreR: return compute_PLmloop10(x,type);
	case MType::MreO: return compute_PMmloop10(x,type);
	case MType::MreR: return compute_PMmloop10(x,type);
	case MType::LMreR: return compute_PLmloop10(x,type);
	case MType::LMorO: return compute_PLmloop10(x,type);
    }
    UNREACHABLE();
}
/**
 * As this just calls the other functions, we can reduce it to just a PX.
 * ptype closing can be done because lend and rend give the correct indices
 * for each recurrence. The functions coming after just use their respective
 * PX functions
 */
void W_final_pf::compute_PX(const Index4D &x, MType type){
    pf_t contributions = 0;
	Matrix4DPF &PX = PX_by_mtype(type);

	const int ptype_closing = pair[S_[x.lend(type)]][S_[x.rend(type)]];

	if (ptype_closing>0){
		contributions += calc_PXiloop(x, type);
		contributions += calc_PXmloop(x,type);

		if(type == MType::Os){
			if(x.i()==x.j() && x.k()==x.l()){
				contributions+=gamma2(x.l(),x.i());
			}
		} else if (x.difference(type)>TURN){
			Index4D xp(x);
			xp.shrink(type);
			Matrix4DPF &PfromX = PfromX_by_mtype(type);
			contributions += PfromX.get(xp)*penalty(xp, gamma2, type);
		}
	}
	PX.setI(x, contributions);
}
/**
 * These should probably just be in a matrix at this point; there are so many arrays, why have four less and recompute every time;
 * It would mean not having four more though, but 11. Remember to probably switch this!!
 */
pf_t W_final_pf::calc_PLiloop(const Index4D &x, MType type){
	if(impossible_case(x)) return 0;
	if(!can_pair(x.lend(type),x.rend(type))) return 0;
	const cand_pos_t i = x.i(), j = x.j(), k = x.k(), l = x.l();

	Matrix4DPF &PX = PX_by_mtype(type);
	pf_t contributions = 0;
	if (i+TURN+2<j) { 
		contributions += PX.get(i+1,j-1,k,l)*get_e_stP(i,j);
	}
	cand_pos_t max_d = std::min(j,i+MAXLOOP);
	for(cand_pos_t d= i+1; d<max_d; ++d){
		cand_pos_t min_dp = std::max(d+TURN,j-MAXLOOP);
		for(cand_pos_t dp = j-1; dp > min_dp; --dp){
			if (!can_pair(d,dp)) continue;
			contributions += get_e_intP(i,d,dp,j)*PX.get(d,dp,k,l);
		}
	}
	return contributions;
}
pf_t W_final_pf::calc_PMiloop(const Index4D &x, MType type){
	if(impossible_case(x)) return 0;
	if(!can_pair(x.lend(type),x.rend(type))) return 0;
	const cand_pos_t i = x.i(), j = x.j(), k = x.k(), l = x.l();

	Matrix4DPF &PX = PX_by_mtype(type);
	pf_t contributions = 0;
	if (i<j && k<l) {
		contributions += PX.get(i,j-1,k+1,l)*get_e_stP(j-1,k+1);
	}
	cand_pos_t max_d = std::max(i,j-MAXLOOP);
	for(cand_pos_t d= j-1; d>max_d; --d){
		cand_pos_t min_dp = std::min(l,k+MAXLOOP); // could switch these here so that we are increasing in the first for like all the others
		for (cand_pos_t dp=k+1; dp <min_dp; ++dp) {
			if (!can_pair(d,dp)) continue;
			contributions += get_e_intP(d,j,k,dp)*PX.get(i,d,dp,l);
		}
	}
	return contributions;
}
pf_t W_final_pf::calc_PRiloop(const Index4D &x, MType type){
	if(impossible_case(x)) return 0;
	if(!can_pair(x.lend(type),x.rend(type))) return 0;
	const cand_pos_t i = x.i(), j = x.j(), k = x.k(), l = x.l();

	Matrix4DPF &PX = PX_by_mtype(type);
	pf_t contributions = 0;
	if (k+TURN+2<l) { 
		contributions += PX.get(i,j,k+1,l-1)*get_e_stP(k,l);
	}
	cand_pos_t max_d = std::min(l,k+MAXLOOP);
	for(cand_pos_t d= k+1; d<max_d; ++d){
		cand_pos_t min_dp = std::max(d+TURN,l-MAXLOOP);
		for(cand_pos_t dp=l-1; dp > min_dp; --dp){
			if (!can_pair(d,dp)) continue;
			contributions += get_e_intP(k,d,dp,l)*PX.get(i,j,d,dp);
		}
	}
	return contributions;
}
pf_t W_final_pf::calc_POiloop(const Index4D &x, MType type){
	if(impossible_case(x)) return 0;
	if(!can_pair(x.lend(type),x.rend(type))) return 0;
	const cand_pos_t i = x.i(), j = x.j(), k = x.k(), l = x.l();

	Matrix4DPF &PX = PX_by_mtype(type);
	pf_t contributions = 0;
	if (i<j && k<l ) { 
		contributions += PX.get(i+1,j,k,l-1)*get_e_stP(i,l);
	}
	cand_pos_t max_d = std::min(j,i+MAXLOOP);
	for(cand_pos_t d= i+1; d<max_d; ++d){
		cand_pos_t min_dp = std::max(l-MAXLOOP,k);
		for (cand_pos_t dp=l-1; dp >min_dp; --dp) {
			if (!can_pair(d,dp)) continue;
			contributions += get_e_intP(i,d,dp,l)*PX.get(d,j,k,dp);
		}
	}
	return contributions;
}
/**
 * 
 * 
 * 
 */
void W_final_pf::compute_PfromL(const Index4D &x, MType type){
	const cand_pos_t i = x.i(), j = x.j(), k = x.k(), l = x.l();
	pf_t contributions = 0;

	Matrix4DPF &PfromXprime = PfromXprime_by_mtype(type);
	for(cand_pos_t d=i+1; d<=j; ++d){
		contributions += calc_WP(i,d-1)*PfromXprime.get(d,j,k,l);
	}

	Matrix4DPF &PfromX = PfromX_by_mtype(type);
	PfromX.set(i,j,k,l,contributions);
}

void W_final_pf::compute_PfromM(const Index4D &x, MType type){
	const cand_pos_t i = x.i(), j = x.j(), k = x.k(), l = x.l();
	pf_t contributions = 0;

	Matrix4DPF &PfromXprime = PfromXprime_by_mtype(type);
	for(cand_pos_t d=i; d<j; ++d){
		contributions += PfromXprime.get(i,d,k,l)*calc_WP(d+1,j);
	}

	Matrix4DPF &PfromX = PfromX_by_mtype(type);
	PfromX.set(i,j,k,l,contributions);
}

void W_final_pf::compute_PfromR(const Index4D &x, MType type){
	const cand_pos_t i = x.i(), j = x.j(), k = x.k(), l = x.l();
	pf_t contributions = 0;

	Matrix4DPF &PfromXprime = PfromXprime_by_mtype(type);
	for(cand_pos_t d=k+1; d<=l; ++d){
		contributions += calc_WP(k,d-1)*PfromXprime.get(i,j,d,l);
	}

	Matrix4DPF &PfromX = PfromX_by_mtype(type);
	PfromX.set(i,j,k,l,contributions);
}

void W_final_pf::compute_PfromO(const Index4D &x, MType type){
	const cand_pos_t i = x.i(), j = x.j(), k = x.k(), l = x.l();
	pf_t contributions = 0;

	Matrix4DPF &PfromXprime = PfromXprime_by_mtype(type);
	for(cand_pos_t d=i+1; d<=j; ++d){
		contributions += calc_WP(i,d-1)*PfromXprime.get(d,j,k,l);
	}

	Matrix4DPF &PfromX = PfromX_by_mtype(type);
	PfromX.set(i,j,k,l,contributions);
}
/**
 * 
 * 
 * 
 */
void W_final_pf::compute_PfromLprime(const Index4D &x, MType type){
	const cand_pos_t i = x.i(), j = x.j(), k = x.k(), l = x.l();
	pf_t contributions = 0;

	for(cand_pos_t d=i; d<j; ++d){
		contributions += calc_PfromXdoubleprime(i,d,k,l,type)*calc_WP(d+1,j);
	}
	
	Matrix4DPF &PfromXprime = PfromXprime_by_mtype(type);
	PfromXprime.set(i,j,k,l,contributions);
}

void W_final_pf::compute_PfromMprime(const Index4D &x, MType type){
	const cand_pos_t i = x.i(), j = x.j(), k = x.k(), l = x.l();
	pf_t contributions = 0;

	for(cand_pos_t d=k+1; d<=l; ++d){
		contributions += calc_WP(k,d-1)*calc_PfromXdoubleprime(i,j,d,l,type);
	}
	
	Matrix4DPF &PfromXprime = PfromXprime_by_mtype(type);
	PfromXprime.set(i,j,k,l,contributions);
}

void W_final_pf::compute_PfromRprime(const Index4D &x, MType type){
	const cand_pos_t i = x.i(), j = x.j(), k = x.k(), l = x.l();
	pf_t contributions = 0;

	for(cand_pos_t d=k; d<l; ++d){
		contributions += calc_PfromXdoubleprime(i,j,k,d,type)*calc_WP(d+1,l);
	}
	
	Matrix4DPF &PfromXprime = PfromXprime_by_mtype(type);
	PfromXprime.set(i,j,k,l,contributions);
}

void W_final_pf::compute_PfromOprime(const Index4D &x, MType type){
	const cand_pos_t i = x.i(), j = x.j(), k = x.k(), l = x.l();
	pf_t contributions = 0;

	for(cand_pos_t d=k; d<l; ++d){
		contributions += calc_PfromXdoubleprime(i,j,k,d,type)*calc_WP(d+1,l);
	}
	
	Matrix4DPF &PfromXprime = PfromXprime_by_mtype(type);
	PfromXprime.set(i,j,k,l,contributions);
}
/**
 * 
 * 
 * 
 */
void W_final_pf::compute_PLmloop00(const Index4D &x, MType type){
	pf_t contributions = 0;
	const cand_pos_t i = x.i(), j = x.j(), k = x.k(), l = x.l();

	Matrix4DPF &PXmloop10 = PXmloop10_by_mtype(type);
	Matrix4DPF &PXmloop01 = PXmloop01_by_mtype(type);
	for(cand_pos_t d = i+1; d<=j; ++d){
		contributions += calc_WB(i,d-1)*PXmloop10.get(d,j,k,l);
	}
    for(cand_pos_t d = i; d<j; ++d){
		contributions += PXmloop01.get(i,d,k,l)*calc_WB(d+1,j);
    }
	Matrix4DPF &PXmloop00 = PXmloop00_by_mtype(type);
	PXmloop00.set(i,j,k,l,contributions);
}
void W_final_pf::compute_PMmloop00(const Index4D &x, MType type){
	pf_t contributions = 0;
	const cand_pos_t i = x.i(), j = x.j(), k = x.k(), l = x.l();

	Matrix4DPF &PXmloop10 = PXmloop10_by_mtype(type);
	Matrix4DPF &PXmloop01 = PXmloop01_by_mtype(type);
    for(cand_pos_t d=i; d<j; ++d){
        contributions += PXmloop10.get(i,d,k,l)*calc_WB(d+1,j);
    }
    for(cand_pos_t d=k+1; d<=l; ++d){
        contributions += PXmloop01.get(i,j,d,l)*calc_WB(k,d-1);
    }
	Matrix4DPF &PXmloop00 = PXmloop00_by_mtype(type);
	PXmloop00.set(i,j,k,l,contributions);
}
void W_final_pf::compute_PRmloop00(const Index4D &x, MType type){
	pf_t contributions = 0;
	const cand_pos_t i = x.i(), j = x.j(), k = x.k(), l = x.l();

	Matrix4DPF &PXmloop10 = PXmloop10_by_mtype(type);
	Matrix4DPF &PXmloop01 = PXmloop01_by_mtype(type);
	for(cand_pos_t d=k+1; d<=l; ++d){
		contributions += calc_WB(k,d-1)*PXmloop10.get(i,j,d,l);
	}
	for(cand_pos_t d=k; d<l; ++d){
		contributions += PXmloop01.get(i,j,k,d)*calc_WB(d+1,l);
    }
	Matrix4DPF &PXmloop00 = PXmloop00_by_mtype(type);
	PXmloop00.set(i,j,k,l,contributions);
}
void W_final_pf::compute_POmloop00(const Index4D &x, MType type){
	pf_t contributions = 0;
	const cand_pos_t i = x.i(), j = x.j(), k = x.k(), l = x.l();

	Matrix4DPF &PXmloop10 = PXmloop10_by_mtype(type);
	Matrix4DPF &PXmloop01 = PXmloop01_by_mtype(type);
    for(cand_pos_t d=i+1; d<=j; ++d){
        contributions += calc_WB(i,d-1)*PXmloop10.get(d,j,k,l);
    }
    for(cand_pos_t d=k; d<l; ++d){
        contributions += PXmloop01.get(i,j,k,d)*calc_WB(d+1,l);
    }
	Matrix4DPF &PXmloop00 = PXmloop00_by_mtype(type);
	PXmloop00.set(i,j,k,l,contributions);
}
/**
 * 
 * 
 * 
 */
void W_final_pf::compute_PLmloop10(const Index4D &x, MType type){
	pf_t contributions = 0;
	const cand_pos_t i = x.i(), j = x.j(), k = x.k(), l = x.l();

	Matrix4DPF &PX = PX_by_mtype(type);
	for(cand_pos_t d = i+1; d <= j; ++d){
        contributions += PX.get(i,d,k,l)*calc_WB(d+1,j);
    }
	Matrix4DPF &PXmloop10 = PXmloop10_by_mtype(type);
	PXmloop10.set(i,j,k,l,contributions);
}
void W_final_pf::compute_PMmloop10(const Index4D &x, MType type){
	pf_t contributions = 0;
	const cand_pos_t i = x.i(), j = x.j(), k = x.k(), l = x.l();

	Matrix4DPF &PX = PX_by_mtype(type);
    for(cand_pos_t d = i+1; d <=j; ++d){
        contributions += PX.get(i,d,k,l)*calc_WB(d+1,j);
    }
	Matrix4DPF &PXmloop10 = PXmloop10_by_mtype(type);
	PXmloop10.set(i,j,k,l,contributions);
}
void W_final_pf::compute_PRmloop10(const Index4D &x, MType type){
	pf_t contributions = 0;
	const cand_pos_t i = x.i(), j = x.j(), k = x.k(), l = x.l();

	Matrix4DPF &PX = PX_by_mtype(type);
    for(cand_pos_t d = k+1; d <= l; ++d){
        contributions += PX.get(i,j,k,d)*calc_WB(d+1,l);
    }
	Matrix4DPF &PXmloop10 = PXmloop10_by_mtype(type);
	PXmloop10.set(i,j,k,l,contributions);
}
void W_final_pf::compute_POmloop10(const Index4D &x, MType type){
	pf_t contributions = 0;
	const cand_pos_t i = x.i(), j = x.j(), k = x.k(), l = x.l();

	Matrix4DPF &PX = PX_by_mtype(type);
	for(cand_pos_t d=i+1; d<=j;++d){
        contributions += PX.get(i,j,k,d)*calc_WB(d+1,l);
    }
	Matrix4DPF &PXmloop10 = PXmloop10_by_mtype(type);
	PXmloop10.set(i,j,k,l,contributions);
}
 /**
  * 
  * 
  * 
  */
void W_final_pf::compute_PLmloop01(const Index4D &x, MType type){
	pf_t contributions = 0;
	const cand_pos_t i = x.i(), j = x.j(), k = x.k(), l = x.l();

	Matrix4DPF &PX = PX_by_mtype(type);
	for(cand_pos_t d = i; d < j; ++d){
        contributions += expcp_pen[d-i]*PX.get(d,j,k,l);
    }

	Matrix4DPF &PXmloop01 = PXmloop01_by_mtype(type);
	PXmloop01.set(i,j,k,l,contributions);
}
void W_final_pf::compute_PMmloop01(const Index4D &x, MType type){
	pf_t contributions = 0;
	const cand_pos_t i = x.i(), j = x.j(), k = x.k(), l = x.l();

	Matrix4DPF &PX = PX_by_mtype(type);
	for(cand_pos_t d = k; d < l; ++d){
        contributions += expcp_pen[d-k]*PX.get(i,j,d,l);
    }

	Matrix4DPF &PXmloop01 = PXmloop01_by_mtype(type);
	PXmloop01.set(i,j,k,l,contributions);
}
void W_final_pf::compute_PRmloop01(const Index4D &x, MType type){
	pf_t contributions = 0;
	const cand_pos_t i = x.i(), j = x.j(), k = x.k(), l = x.l();

	Matrix4DPF &PX = PX_by_mtype(type);
    for(cand_pos_t d = k; d < l; d++){
        contributions += expcp_pen[d-k]*PX.get(i,j,d,l);
    }

	Matrix4DPF &PXmloop01 = PXmloop01_by_mtype(type);
	PXmloop01.set(i,j,k,l,contributions);
}
void W_final_pf::compute_POmloop01(const Index4D &x, MType type){
	pf_t contributions = 0;
	const cand_pos_t i = x.i(), j = x.j(), k = x.k(), l = x.l();

	Matrix4DPF &PX = PX_by_mtype(type);
	for(cand_pos_t d = i; d < j; ++d){
        contributions += expcp_pen[d-i]*PX.get(d,j,k,l);
    }

	Matrix4DPF &PXmloop01 = PXmloop01_by_mtype(type);
	PXmloop01.set(i,j,k,l,contributions);
}
/**
 * 
 * 
 * 
 */
pf_t W_final_pf::calc_PfromLdoubleprime(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l){
	if (!(i <= j && j < k-1 && k <= l)){
		return 0;
	}
	assert(!(i<=0 || l> n));
	pf_t contributions = 0;
	contributions += PR.get(i,j,k,l)*gamma2(l,k)*expPB_penalty;
	contributions += PM.get(i,j,k,l)*gamma2(j,k)*expPB_penalty;
	contributions += POm.get(i,j,k,l)*gamma2(l,i)*expPB_penalty;
	contributions += POs.get(i,j,k,l)*gamma2(l,i)*expPB_penalty;

	return contributions;
}
pf_t W_final_pf::calc_PfromMdoubleprime(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l){
	if (!(i <= j && j < k-1 && k <= l)){
		return 0;
	}
	assert(!(i<=0 || l> n));
	pf_t contributions = 0;
	contributions += PL.get(i,j,k,l)*gamma2(j,i)*expPB_penalty;
	contributions += PR.get(i,j,k,l)*gamma2(l,k)*expPB_penalty;
	contributions += POm.get(i,j,k,l)*gamma2(l,i)*expPB_penalty;

	return contributions;
}
pf_t W_final_pf::calc_PfromRdoubleprime(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l){
	if (!(i <= j && j < k-1 && k <= l)){
		return 0;
	}
	assert(!(i<=0 || l> n));
	pf_t contributions = 0;
	contributions += PM.get(i,j,k,l)*gamma2(j,k)*expPB_penalty;
	contributions += POm.get(i,j,k,l)*gamma2(l,i)*expPB_penalty;
	contributions += POs.get(i,j,k,l)*gamma2(l,i)*expPB_penalty;

	return contributions;
}
pf_t W_final_pf::calc_PfromOdoubleprime(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l){
	if (!(i <= j && j < k-1 && k <= l)){
		return 0;
	}
	assert(!(i<=0 || l> n));
	pf_t contributions = 0;
	contributions += PL.get(i,j,k,l)*gamma2(j,i)*expPB_penalty;
	contributions += PR.get(i,j,k,l)*gamma2(l,k)*expPB_penalty;

	return contributions;
}
pf_t W_final_pf::calc_PfromLreOdoubleprime(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l){
	if (!(i <= j && j < k-1 && k <= l)){
		return 0;
	}
	assert(!(i<=0 || l> n));
	pf_t contributions = 0;
	contributions += PMreO.get(i,j,k,l)*gamma2(j,k)*expPB_penalty;
	contributions += POm.get(i,j,k,l)*gamma2(l,i)*expPB_penalty;
	contributions += POs.get(i,j,k,l)*gamma2(l,i)*expPB_penalty;

	return contributions;
}
pf_t W_final_pf::calc_PfromLreRdoubleprime(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l){
	if (!(i <= j && j < k-1 && k <= l)){
		return 0;
	}
	assert(!(i<=0 || l> n));
	pf_t contributions = 0;
	contributions += PMreR.get(i,j,k,l)*gamma2(j,k)*expPB_penalty;
	contributions += PR.get(i,j,k,l)*gamma2(l,k)*expPB_penalty;

	return contributions;
}
pf_t W_final_pf::calc_PfromMreOdoubleprime(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l){
	if (!(i <= j && j < k-1 && k <= l)){
		return 0;
	}
	assert(!(i<=0 || l> n));
	pf_t contributions = 0;
	contributions += PLreO.get(i,j,k,l)*gamma2(j,i)*expPB_penalty;
	contributions += POm.get(i,j,k,l)*gamma2(l,i)*expPB_penalty;
	contributions += POs.get(i,j,k,l)*gamma2(l,i)*expPB_penalty;

	return contributions;
}
pf_t W_final_pf::calc_PfromMreRdoubleprime(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l){
	if (!(i <= j && j < k-1 && k <= l)){
		return 0;
	}
	assert(!(i<=0 || l> n));
	pf_t contributions = 0;
	contributions += PLreR.get(i,j,k,l)*gamma2(j,i)*expPB_penalty;
	contributions += PR.get(i,j,k,l)*gamma2(l,k)*expPB_penalty;

	return contributions;
}
pf_t W_final_pf::calc_PfromLMreRdoubleprime(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l){
	if (!(i <= j && j < k-1 && k <= l)){
		return 0;
	}
	assert(!(i<=0 || l> n));
	pf_t contributions = 0;
	contributions += PMreR.get(i,j,k,l)*gamma2(j,k)*expPB_penalty;

	return contributions;
}
pf_t W_final_pf::calc_PfromLMorOdoubleprime(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l){
	if (!(i <= j && j < k-1 && k <= l)){
		return 0;
	}
	assert(!(i<=0 || l> n));
	pf_t contributions = 0;
	contributions += PM.get(i,j,k,l)*gamma2(j,k)*expPB_penalty;
	contributions += POm.get(i,j,k,l)*gamma2(l,i)*expPB_penalty;
	contributions += POs.get(i,j,k,l)*gamma2(l,i)*expPB_penalty;

	return contributions;
}
/**
 * 
 * 
 * 
 */
pf_t W_final_pf::calc_WB(cand_pos_t i, cand_pos_t l){
	if (i<=0 || l<=0 || i>n || l>n){
		return 0;
	}
	if (i>l) return 1;
	return expcp_pen[l-i+1] + WBP.get(i,l);
}
pf_t W_final_pf::calc_WP(cand_pos_t i, cand_pos_t l){
	if (i<=0 || l<=0 || i>n || l>n){
		return 0;
	}
	if (i>l) return 1; // needed as this will happen 
	return expPUP_pen[l-i+1] + WPP.get(i,l);
}
/**
 * 
 * 
 * 
 */
pf_t W_final_pf::compute_int(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l) {
    const pair_type ptype_closing = pair[S_[i]][S_[j]];
    return exp_E_IntLoop(k - i - 1, j - l - 1, ptype_closing, rtype[pair[S_[k]][S_[l]]], S1_[i + 1], S1_[j - 1], S1_[k - 1], S1_[l + 1], exp_params_);
}

pf_t W_final_pf::get_e_stP(cand_pos_t i, cand_pos_t j) {
    if (i + 1 == j - 1) { // TODO: do I need something like that or stack is taking care of this?
        return 0;
    }
    pf_t e_st = compute_int(i, j, i + 1, j - 1);

    return pow(e_st, e_stP_penalty);
}

pf_t W_final_pf::get_e_intP(cand_pos_t i, cand_pos_t ip, cand_pos_t jp, cand_pos_t j) {
    if (ip == i + 1 && jp == j - 1) return 0;
    pf_t e_int = compute_int(i, j, ip, jp);

    return pow(e_int, e_intP_penalty);
}

void W_final_pf::fill_structure(std::vector<int> &fres, std::string &structure){
	std::stack < brack_type > st;

    st.push(brack_type('<','>'));
    st.push(brack_type('{','}'));
    st.push(brack_type('[',']'));
    st.push(brack_type('(',')'));

    cand_pos_t isInABand=0;
    // cand_pos_t num_crossing_bands=0;

    std::list <band_elem > bands;
    bands.push_back(band_elem('|','|',0,0,0,0));

    for (cand_pos_t i = 0; i < n; i++){
        cand_pos_t j = fres[i];
        if (j == -1){ // i is unpaired
            structure[i]='.';
        }else if (i < j){
            isInABand=0;
            //			for (band_elem *p = head; p != NULL;  p = p->next){
            for (std::list<band_elem > ::iterator it = bands.begin(); it != bands.end(); it++){
                if(i> (*it).inner_start && j < (*it).inner_end){ // i.e. i is paired and i.pair[i] is nested in the band
                    (*it).inner_start = i;
                    (*it).inner_end =j;
                    structure[i] = (*it).open;
                    structure[j] = (*it).close;
                    isInABand=1;
                    break;
                }
            }

            if (!isInABand){
                brack_type e = st.top();
                st.pop();
                // num_crossing_bands++;

                bands.push_back(band_elem(e.open,e.close,i,j,i,j));
                structure[i] = e.open;
                structure[j] = e.close;
            }

        }else{ //having the closing base pair i>pair[i]
            for(std::list<band_elem > ::iterator current = bands.begin(); current != bands.end(); current++){
                if (i == (*current).outer_end){
                    st.push(brack_type((*current).open,(*current).close));
                    break;
                }
            }
        }
    }
}

char W_final_pf::bpp_symbol(pf_t *P) {
    if (P[0] > 0.667) return '.';
    if (P[0] > (P[1] + P[2] + P[3] + P[4])) return ',';

    if (P[1] > 0.667) return '(';
    if (P[2] > 0.667) return ')';
    if (P[3] > 0.667) return '[';
    if (P[4] > 0.667) return ']';
    if(P[1]+P[2] > P[3]+P[4]){
        if ((P[1] + P[2]) > P[0]) {
            if ((P[1] / (P[1] + P[2])) > 0.667) return '{';

            if ((P[2] / (P[1] + P[2])) > 0.667)
                return '}';
            else return '|';
        }
    } else{ 
        if ((P[3] + P[4]) > P[0]) {
            if ((P[3] / (P[3] + P[4])) > 0.667) return '/';

            if ((P[4] / (P[3] + P[4])) > 0.667)
                return '\\';
            else return '|';
        }
    }
    return ':';
}
void W_final_pf::pairing_tendency() {

    for (cand_pos_t j = 1; j <= n; j++) {
        pf_t P[5] = {1, 0, 0, 0, 0}; // unpaired, PK-free left, PK-free right, PK left, PK right
        for (cand_pos_t i = 1; i < j; i++) {
            // bool weakly_closed_ij = tree.weakly_closed(i, j);
            std::pair<cand_pos_tu, cand_pos_tu> base_pair(i, j);
            pf_t probability_ij = (pf_t)samples[base_pair] / num_samples;
            // if(weakly_closed_ij) P[2] += probability_ij; else P[4] += probability_ij;
			P[2] += probability_ij;
            P[0] -= probability_ij;
        }
        for (cand_pos_t i = j + 1; i <= n; i++) {
            // bool weakly_closed_ji = tree.weakly_closed(j, i);
            std::pair<cand_pos_tu, cand_pos_tu> base_pair(j, i);
            pf_t probability_ji = (pf_t)samples[base_pair] / num_samples;
            // if(weakly_closed_ji) P[1] += probability_ji; else P[3] += probability_ji;
			P[1] += probability_ji;
            P[0] -= probability_ji;
        }
        structure[j - 1] = bpp_symbol(P);
    }
}

// bool pair_is_pseudoknotted(cand_pos_t i, cand_pos_t j, const std::vector<std::tuple<cand_pos_t, cand_pos_t, pf_t>> &pairs)
// {
//     pf_t crossing_mass = 0;
//     pf_t nesting_mass  = 0;

//     for (auto &[c, d, prob_cd] : pairs) {
//         if (c == i && d == j) continue;
//         cand_pos_t k = std::min(c,d), l = std::max(c,d);

//         bool crosses = (i < k && k < j && j < l) ||
//                        (k < i && i < l && l < j);
//         bool nested  = (i <= k && l <= j) || (k <= i && j <= l);

//         if (crosses) crossing_mass += prob_cd;
//         if (nested)  nesting_mass  += prob_cd;
//     }
//     // PK if more probability mass crosses it than nests with it
//     return crossing_mass > nesting_mass;
// }
// void W_final_pf::pairing_tendency() {
//     // Step 1: collect all base pairs with non-trivial probability
//     const pf_t threshold = 0.01;
    
//     std::vector<std::tuple<cand_pos_t, cand_pos_t, pf_t>> pairs;
//     for (auto &kv : samples) {
//         pf_t prob = (pf_t)kv.second / num_samples;
//         if (prob >= threshold) {
//             pairs.emplace_back(kv.first.first, kv.first.second, prob);
//         }
//     }

//     // Step 2: for each position j, determine PK status of each pair involving j
//     for (cand_pos_t j = 1; j <= n; j++) {
//         pf_t P[5] = {1, 0, 0, 0, 0};

//         for (cand_pos_t i = 1; i < j; i++) {
//             std::pair<cand_pos_tu, cand_pos_tu> bp(i, j);
//             auto it = samples.find(bp);
//             if (it == samples.end()) continue;

//             pf_t prob_ij = (pf_t)it->second / num_samples;
//             if (prob_ij == 0) continue;

//             bool is_pk = pair_is_pseudoknotted(i, j, pairs);
//             if (is_pk) P[4] += prob_ij; else P[2] += prob_ij;
//             P[0] -= prob_ij;
//         }
//         for (cand_pos_t i = j + 1; i <= n; i++) {
//             std::pair<cand_pos_tu, cand_pos_tu> bp(j, i);
//             auto it = samples.find(bp);
//             if (it == samples.end()) continue;

//             pf_t prob_ji = (pf_t)it->second / num_samples;
//             if (prob_ji == 0) continue;

//             bool is_pk = pair_is_pseudoknotted(j, i, pairs);
//             if (is_pk) P[3] += prob_ji; else P[1] += prob_ji;
//             P[0] -= prob_ji;
//         }
//         structure[j - 1] = bpp_symbol(P);
//     }
// }