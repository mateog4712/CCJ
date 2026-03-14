#include "part_func.hh"
// #include "dot_plot.hh"
#include "h_externs.hh"
#include "pf_globals.hh"

#include <algorithm>
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <string>
#include <cassert>

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

W_final_pf::W_final_pf(std::string &seq, std::string &MFE_structure, double MFE_energy, int dangle, int num_samples, bool PSplot)
    : exp_params_(scale_pf_parameters()) {
    this->seq = seq;
    this->MFE_structure = MFE_structure;
    this->n = seq.length();
    this->PSplot = PSplot;
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
	TriangleMatrix_PF::new_index(index,n+1);
	Matrix4DPF::construct_index(index3D,n);
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
	PK.init(n,index3D);
	PL.init(n,index3D);
	PR.init(n,index3D);
	PM.init(n,index3D);
	PO.init(n,index3D);
	PfromL.init(n,index3D);
	PfromR.init(n,index3D);
	PfromM.init(n,index3D);
	PfromO.init(n,index3D);

	PLmloop00.init(n,index3D);
	PLmloop01.init(n,index3D);
	PLmloop10.init(n,index3D);

	PRmloop00.init(n,index3D);
	PRmloop01.init(n,index3D);
	PRmloop10.init(n,index3D);

	PMmloop00.init(n,index3D);
	PMmloop01.init(n,index3D);
	PMmloop10.init(n,index3D);

	POmloop00.init(n,index3D);
	POmloop01.init(n,index3D);
	POmloop10.init(n,index3D);

    rescale_pk_globals();
    exp_params_rescale(MFE_energy);
    W.resize(n+1,scale[1]);
}

W_final_pf::~W_final_pf() {}

void W_final_pf::exp_params_rescale(double mfe) {
    double e_per_nt, kT;
    kT = exp_params_->kT;

    e_per_nt = mfe * 1000. / this->n;

    exp_params_->pf_scale = exp(-(exp_params_->model_details.sfact * e_per_nt) / kT);

    if (exp_params_->pf_scale < 1.) exp_params_->pf_scale = 1.;

    exp_params_->pf_scale = 1.; // My scaling is bad so I need this to be 1 right now.

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

    expPS_penalty = RESCALE_BF(PS_penalty, PS_penalty * 3, TT, kT);
    expPSM_penalty = RESCALE_BF(PSM_penalty, PSM_penalty * 3, TT, kT);
    expPSP_penalty = RESCALE_BF(PSP_penalty, PSP_penalty * 3, TT, kT);
    expPB_penalty = RESCALE_BF(PB_penalty, PB_penalty * 3, TT, kT);
    expPUP_penalty = RESCALE_BF(PUP_penalty, PUP_penalty * 3, TT, kT);
    expPPS_penalty = RESCALE_BF(PPS_penalty, PPS_penalty * 3, TT, kT);

    expa_penalty = RESCALE_BF(a_penalty, ML_closingdH, TT, kT);
    expb_penalty = RESCALE_BF(b_penalty, ML_interndH, TT, kT);
    expc_penalty = RESCALE_BF(c_penalty, ML_BASEdH, TT, kT);

    expap_penalty = RESCALE_BF(ap_penalty, ap_penalty * 3, TT, kT);
    expbp_penalty = RESCALE_BF(bp_penalty, bp_penalty * 3, TT, kT);
    expcp_penalty = RESCALE_BF(cp_penalty, cp_penalty * 3, TT, kT);
}

inline pf_t W_final_pf::to_Energy(pf_t energy, cand_pos_t length) {
    return ((-log(energy) - length * log(exp_params_->pf_scale)) * exp_params_->kT / 1000.0);
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
			contributions += acc*get_energy(k, j)*exp_Extloop(k, j);
			contributions += acc*P.get(k,j)*expPS_penalty;
		}
		W[j] = contributions;
	}
    pf_t energy = to_Energy(W[n], n);

    structure = std::string(n, '.');

    return energy;
}

pf_t W_final_pf::exp_Extloop(cand_pos_t i, cand_pos_t j) {
    pair_type tt = pair[S_[i]][S_[j]];

    if (exp_params_->model_details.dangles == 1 || exp_params_->model_details.dangles == 2) {
        base_type si1 = i > 1 ? S_[i - 1] : -1;
        base_type sj1 = j < n ? S_[j + 1] : -1;
        return exp_E_ExtLoop(tt, si1, sj1, exp_params_);
    } else {
        return exp_E_ExtLoop(tt, -1, -1, exp_params_);
    }
}

pf_t W_final_pf::exp_MLstem(cand_pos_t i, cand_pos_t j) {
    pair_type tt = pair[S_[i]][S_[j]];
    if (exp_params_->model_details.dangles == 1 || exp_params_->model_details.dangles == 2) {
        base_type si1 = i > 1 ? S_[i - 1] : -1;
        base_type sj1 = j < n ? S_[j + 1] : -1;
        return exp_E_MLstem(tt, si1, sj1, exp_params_);
    } else {
        return exp_E_MLstem(tt, -1, -1, exp_params_);
    }
}

pf_t W_final_pf::exp_Mbloop(cand_pos_t i, cand_pos_t j) {
    pair_type tt = pair[S_[j]][S_[i]];
    if (exp_params_->model_details.dangles == 1 || exp_params_->model_details.dangles == 2) {
        base_type si1 = i > 1 ? S_[i + 1] : -1;
        base_type sj1 = j < n ? S_[j - 1] : -1;
        return exp_E_MLstem(tt, sj1, si1, exp_params_);
    } else {
        return exp_E_MLstem(tt, -1, -1, exp_params_);
    }
}

pf_t W_final_pf::HairpinE(cand_pos_t i, cand_pos_t j) {
    const int ptype_closing = pair[S_[i]][S_[j]];
    if (ptype_closing == 0) return 0;
    pf_t e_h = static_cast<pf_t>(exp_E_Hairpin(j - i - 1, ptype_closing, S1_[i + 1], S1_[j - 1], &seq.c_str()[i - 1], exp_params_));
    e_h *= scale[j - i + 1];
    return e_h;
}

pf_t W_final_pf::compute_internal(cand_pos_t i, cand_pos_t j) {
    pf_t v_iloop = 0;
    cand_pos_t max_k = std::min(j - TURN - 2, i + MAXLOOP + 1);
    const pair_type ptype_closing = pair[S_[i]][S_[j]];
    for (cand_pos_t k = i + 1; k <= max_k; ++k) {
        cand_pos_t min_l = std::max(k + TURN + 1 + MAXLOOP + 2, k + j - i) - MAXLOOP - 2;
        for (cand_pos_t l = j - 1; l >= min_l; --l) {
            pf_t v_iloop_kl = get_energy(k, l)
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

    WMv_contributions += (get_energy(i, j) * exp_MLstem(i, j));
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
        pf_t qbt1 = get_energy(k, j) * exp_MLstem(k, j);
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
		for(cand_pos_t k = l; k>=j+2; --k){

			compute_PLmloop00(i,j,k,l);

			compute_PLmloop01(i,j,k,l);

			compute_PLmloop10(i,j,k,l);

			compute_PRmloop00(i,j,k,l);

			compute_PRmloop01(i,j,k,l);

			compute_PRmloop10(i,j,k,l);

			compute_PMmloop00(i,j,k,l);

			compute_PMmloop01(i,j,k,l);

			compute_PMmloop10(i,j,k,l);

			compute_POmloop00(i,j,k,l);

			compute_POmloop01(i,j,k,l);

			compute_POmloop10(i,j,k,l);

			compute_PL(i,j,k,l);

			compute_PR(i,j,k,l);

			compute_PM(i,j,k,l);

			compute_PO(i,j,k,l);

			compute_PfromL(i,j,k,l);

			compute_PfromR(i,j,k,l);

			compute_PfromM(i,j,k,l);

			compute_PfromO(i,j,k,l);

			compute_PK(i,j,k,l);

		}
	}
}

void W_final_pf::compute_WBP(cand_pos_t i, cand_pos_t l){
	pf_t contributions = 0;

	for(cand_pos_t d=i; d< l; ++d){
		contributions += get_energy(d,l)*beta2P(l,d)*expPPS_penalty;
		contributions += P.get(d,l)*expPSM_penalty*expPPS_penalty;
	}
	contributions+=WBP.get(i,l-1)*expcp_pen[1];
	WBP.set(i,l) = contributions;
}

void W_final_pf::compute_WPP(cand_pos_t i, cand_pos_t l){
	pf_t contributions = 0;

	for(cand_pos_t d=i; d<l; ++d){
		contributions += get_WP(i,d-1)*get_energy(d,l)*gamma2(l,d)*expPPS_penalty;
		contributions += get_WP(i,d-1)*P.get(d,l)*expPSP_penalty*expPPS_penalty;
	}
	contributions+=WBP.get(i,l-1)*expPUP_pen[1];
	WPP.set(i,l) = contributions;
}

void W_final_pf::compute_P(cand_pos_t i, cand_pos_t l){
    pf_t contributions = 0;
	for(cand_pos_t j=i; j<l; ++j){
		for (cand_pos_t d=j+1; d<l; ++d){
			for (cand_pos_t k=d+1; k<l; ++k){
				contributions += PK.get(i,j,d+1,k)*PK.get(j+1,d,k+1,l);
			}
		}
	}
	P.set(i,l) = contributions;
}

void W_final_pf::compute_PK(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	pf_t contributions = 0;

	for(cand_pos_t d=i+1; d< j; ++d){
		contributions += PK.get(i,d,k,l)*get_WP(d+1,j);  //12G1
	}

	for(cand_pos_t d=k+1; d< l; ++d){
		contributions += PK.get(i,j,d,l)*get_WP(k,d-1); //1G21
	}

	contributions += PL.get(i,j,k,l)*gamma2(j,i)*expPB_penalty;
	contributions += PM.get(i,j,k,l)*gamma2(j,k)*expPB_penalty;
	contributions += PR.get(i,j,k,l)*gamma2(l,k)*expPB_penalty;
	contributions += PO.get(i,j,k,l)*gamma2(l,i)*expPB_penalty;
	
    PK.set(i,j,k,l,contributions);
}

void W_final_pf::compute_PL(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	pf_t contributions = 0;

	const int ptype_closing = pair[S_[i]][S_[j]];

	if (ptype_closing>0){

		contributions += get_PLiloop(i,j,k,l);

		contributions += get_PLmloop(i,j,k,l)*expbp_penalty;

		if (j>=(i+TURN+1)){
			contributions += PfromL.get(i+1,j-1,k,l)*gamma2(j,i);
		}
	}
	PL.set(i,j,k,l,contributions);
}

void W_final_pf::compute_PR(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	pf_t contributions = 0;

	const int ptype_closing = pair[S_[k]][S_[l]];

	if (ptype_closing>0){
		contributions += get_PRiloop(i,j,k,l);

		contributions += get_PRmloop(i,j,k,l)*expbp_penalty;

		if (l>=(k+TURN+1)){
			contributions += PfromR.get(i,j,k+1,l-1)*gamma2(l,k);
		}
	}
	PR.set(i,j,k,l,contributions);
}

void W_final_pf::compute_PM(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	pf_t contributions = 0;

	const int ptype_closing = pair[S_[j]][S_[k]];
	if (ptype_closing>0){
		contributions += get_PMiloop(i,j,k,l);

		contributions += get_PMmloop(i,j,k,l)*expbp_penalty;

		if (k>=(j+TURN-1)){
			contributions += PfromM.get(i,j-1,k+1,l)*gamma2(j,k);
		}

		if(i==j && k==l){
			contributions += gamma2(i,l);
		}
	}
	PM.set(i,j,k,l,contributions);
}

void W_final_pf::compute_PO(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	pf_t contributions = 0;

	const int ptype_closing = pair[S_[i]][S_[l]];

	if (ptype_closing>0){
		contributions += get_POiloop(i,j,k,l);

		contributions += get_POmloop(i,j,k,l)*expbp_penalty;

		// Hosna, July 11, 2014
		// To avoid addition of close base pairs we check for the following here
		if (l>=(i+TURN+1)){
			contributions += PfromO.get(i+1,j,k,l-1)*gamma2(l,i);
		}
	}
	PO.set(i,j,k,l,contributions);
}

void W_final_pf::compute_PfromL(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	pf_t contributions = 0;

	for(cand_pos_t d=i+1; d< j; ++d){
		contributions += PfromL.get(d,j,k,l)*get_WP(i,d-1);
		contributions += PfromL.get(i,d,k,l)*get_WP(d+1,j);
	}

	contributions += PR.get(i,j,k,l)*gamma2(l,k)*expPB_penalty;

	contributions += PM.get(i,j,k,l)*gamma2(j,k)*expPB_penalty;

	contributions += PO.get(i,j,k,l)*gamma2(l,i)*expPB_penalty;

	PfromL.set(i,j,k,l,contributions);
}

void W_final_pf::compute_PfromR(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	pf_t contributions = 0;

	for(cand_pos_t d=k+1; d< l; ++d){
		contributions += PfromR.get(i,j,d,l)*get_WP(k,d-1);
		contributions += PfromR.get(i,j,k,d)*get_WP(d+1,l);
	}

	contributions += PM.get(i,j,k,l)*gamma2(j,k)*expPB_penalty;

	contributions += PO.get(i,j,k,l)*gamma2(l,i)*expPB_penalty;

	PfromR.set(i,j,k,l,contributions);
}

void W_final_pf::compute_PfromM(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	pf_t contributions = 0;

	for(cand_pos_t d=i+1; d<j; ++d){
		contributions += PfromM.get(i,d,k,l)*get_WP(d+1,j);
	}
	for(cand_pos_t d=k+1; d<l; ++d){
		contributions += PfromM.get(i,j,d,l)*get_WP(k,d-1);
	}

	contributions += PL.get(i,j,k,l)*gamma2(j,i)*expPB_penalty;

	contributions += PR.get(i,j,k,l)*gamma2(l,k)*expPB_penalty;

	PfromM.set(i,j,k,l,contributions);
}

void W_final_pf::compute_PfromO(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	pf_t contributions = 0;

	for(cand_pos_t d=i+1; d< j; ++d){
		contributions += PfromO.get(d,j,k,l)*get_WP(i,d-1);
	}
	for(cand_pos_t d=k+1; d<l; ++d){
		contributions += PfromO.get(i,j,k,d)*get_WP(d+1,l);
	}

	contributions += PL.get(i,j,k,l)*gamma2(j,i)*expPB_penalty;

	contributions += PR.get(i,j,k,l)*gamma2(l,k)*expPB_penalty;

    PfromO.set(i,j,k,l,contributions);
}

void W_final_pf::compute_PLmloop00(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	pf_t contributions = 0;

	contributions += PL.get(i,j,k,l)*beta2P(j,i);
    for(cand_pos_t d = i; d<=j; ++d){
        if (d>i){
            contributions += get_WB(i,d-1)*PLmloop00.get(d,j,k,l);
        }
        if(d<j){
            contributions += PLmloop00.get(i,d,k,l)*get_WB(d+1,j);
        }

    }
	PLmloop00.set(i,j,k,l,contributions);	
}

void W_final_pf::compute_PLmloop01(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	pf_t contributions = 0;
	
	for(cand_pos_t d = i; d < j; ++d){
        contributions += PLmloop00.get(i,d,k,l)*WBP.get(d+1,j);
    }
	PLmloop01.set(i,j,k,l,contributions);

}

void W_final_pf::compute_PLmloop10(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	pf_t contributions = 0;

	for(cand_pos_t d = i+1; d <= j; ++d){
        contributions += WBP.get(i,d-1)*PLmloop00.get(d,j,k,l);
        if(d<j){
            contributions += PLmloop10.get(i,d,k,l)*get_WB(d+1,j);
        }
    }
	PLmloop10.set(i,j,k,l,contributions);
}

void W_final_pf::compute_PRmloop00(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	pf_t contributions = 0;

	contributions += PR.get(i,j,k,l)*beta2P(l,k);
	for(cand_pos_t d=k; d<=l; ++d){
        if(d>k){
            contributions += get_WB(k,d-1)*PRmloop00.get(i,j,d,l);
        }
        if (d<l){
            contributions += PRmloop00.get(i,j,k,d)*get_WB(d+1,l);
        }
    }
	PRmloop00.set(i,j,k,l,contributions);
}


void W_final_pf::compute_PRmloop01(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	pf_t contributions = 0;

	contributions += PRmloop01.get(i,j,k,l-1)*expcp_pen[1];
    for(cand_pos_t d = k; d < l; d++){
        contributions += PRmloop00.get(i,j,k,d)*WBP.get(d+1,l);
    }
	PRmloop01.set(i,j,k,l,contributions);
}

void W_final_pf::compute_PRmloop10(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	pf_t contributions = 0;

	contributions += PRmloop10.get(i,j,k+1,l)*expcp_pen[1];
    for(cand_pos_t d = k+1; d <= l; ++d){
        contributions += WBP.get(k,d-1)*PRmloop00.get(i,j,d,l);
    }
	PRmloop10.set(i,j,k,l,contributions);
}

void W_final_pf::compute_PMmloop00(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	pf_t contributions = 0;

	contributions += PM.get(i,j,k,l)*beta2P(j,k);
    for(cand_pos_t d=i; d<j; ++d){
        contributions += PMmloop00.get(i,d,k,l)*get_WB(d+1,j);
    }
    for(cand_pos_t d=k+1; d<=l; ++d){
        contributions += PMmloop00.get(i,j,d,l)*get_WB(k,d-1);
    }
	PMmloop00.set(i,j,k,l,contributions);
}


void W_final_pf::compute_PMmloop01(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	pf_t contributions = 0;

	contributions += PMmloop01.get(i,j,k+1,l)+expcp_pen[1];
	for(cand_pos_t d = k; d < l; ++d){
        contributions += PMmloop00.get(i,j,k,d)*WBP.get(d+1,l);
    }
	PMmloop01.set(i,j,k,l,contributions);
}

void W_final_pf::compute_PMmloop10(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	pf_t contributions = 0;

	contributions += PMmloop10.get(i,j-1,k,l)*expcp_pen[1];
	for(cand_pos_t d=i+1; d<=j;++d){
        contributions += WBP.get(i,d-1)*PMmloop00.get(d,j,k,l);
    }
    for(cand_pos_t d = k+1; d < l; ++d){
        contributions += POmloop10.get(i,j,k,d)*get_WB(d+1,l);
    }
	PMmloop10.set(i,j,k,l,contributions);
}

void W_final_pf::compute_POmloop00(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	pf_t contributions = 0;

	contributions += PO.get(i,j,k,l)*beta2P(l,i);
    for(cand_pos_t d=i+1; d<=j; ++d){
        contributions += get_WB(i,d-1)*POmloop00.get(d,j,k,l);
    }
    for(cand_pos_t d=k; d<l; ++d){
        contributions = POmloop00.get(i,j,k,d)*get_WB(d+1,l);
    }
	POmloop00.set(i,j,k,l,contributions);
}


void W_final_pf::compute_POmloop01(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	pf_t contributions = 0;

	for(cand_pos_t d = k; d < l; ++d){
        contributions += POmloop00.get(i,j,k,d)*WBP.get(d+1,l);
    }
	POmloop01.set(i,j,k,l,contributions);
}

void W_final_pf::compute_POmloop10(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	pf_t contributions = 0;

	for(cand_pos_t d=i+1; d<=j;++d){
        contributions += WBP.get(i,d-1)*POmloop00.get(d,j,k,l);
    }
    for(cand_pos_t d = k+1; d < l; ++d){
        contributions += POmloop10.get(i,j,k,d) + get_WB(d+1,l);
    }

	POmloop10.set(i,j,k,l,contributions);
}

pf_t W_final_pf::get_WB(cand_pos_t i, cand_pos_t j){
	if (i<=0 || j<=0 || i>n || j>n){
		return 0;
	}
	if (i>j) return 1;
	return expcp_pen[j-i+1]+WBP.get(i,j);
}

pf_t W_final_pf::get_WP(cand_pos_t i, cand_pos_t j){
	if (i<=0 || j<=0 || i>n || j>n){
		return 0;
	}
	if (i>j) return 1; // needed as this will happen 
	return expPUP_pen[j-i+1]+WPP.get(i,j);
}

// pf_t W_final_pf::get_PM(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
// 	if (!(i <= j && j < k-1 && k <= l)){
// 		return 0;
// 	}
// 	assert(!(i<=0 || l> n));

// 	int ptype_closing = pair[S_[j]][S_[k]];
// 	if (ptype_closing == 0){
// 		return 0;
// 	}
// 	cand_pos_t ij = index[i]+j-i;
// 	cand_pos_t kl = index[k]+l-k;
// 	if (i ==j && k ==l){
// 		return gamma2(i,l);   ////// Does this suggest I could be missing a case or would it be set to this already?
// 	}

// 	return PM[ij][kl];
// }

pf_t W_final_pf::get_PLiloop(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	if (!(i <= j && j < k-1 && k <= l)){
		return 0;
	}
	assert(!(i<=0 || l> n));

	const int ptype_closing = pair[S_[i]][S_[j]];
	if(ptype_closing == 0) return 0;

	pf_t contributions = 0;
	contributions += PL.get(i+1,j-1,k,l)*get_e_stP(i,j);

	for(cand_pos_t d= i+1; d<std::min(j,i+MAXLOOP); ++d){
		for(cand_pos_t dp = j-1; dp > std::max(d+TURN,j-MAXLOOP); --dp){
			cand_pos_t u1 = d-(i)-1; // check these
			cand_pos_t u2 = (j)-dp-1;
			contributions += get_e_intP(i,d,dp,j)*PL.get(d,dp,k,l)*scale[u1 + u2 + 2];
		}
	}
	return contributions;
}

pf_t W_final_pf::get_PLmloop(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	if (!(i <= j && j < k-1 && k <= l)){
		return 0;
	}
	assert(!(i<=0 || l> n));
    pf_t contributions = 0;
	contributions += PLmloop10.get(i+1,j-1,k,l)*expap_penalty*beta2P(j,i);
    contributions += PLmloop01.get(i+1,j-1,k,l)*ap_penalty*beta2P(j,i);

	return contributions;
}

pf_t W_final_pf::get_PRiloop(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	if (!(i <= j && j < k-1 && k <= l)){
		return 0;
	}
	assert(!(i<=0 || l> n));

	const int ptype_closing = pair[S_[k]][S_[l]];
	if(ptype_closing == 0) return 0;

	pf_t contributions = 0;
	contributions += PR.get(i,j,k+1,l-1)*get_e_stP(k,l);

	for(cand_pos_t d= k+1; d<std::min(l,k+MAXLOOP); ++d){
		for(cand_pos_t dp=l-1; dp > std::max(d+TURN,l-MAXLOOP); --dp){
			cand_pos_t u1 = d-(k)-1; // check these
			cand_pos_t u2 = (l)-dp-1;
			contributions += get_e_intP(k,d,dp,l)*PR.get(i,j,d,dp)*scale[u1 + u2 + 2];
		}
	}
	return contributions;
}

pf_t W_final_pf::get_PRmloop(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	if (!(i <= j && j < k-1 && k <= l)){
		return 0;
	}
	assert(!(i<=0 || l> n));
    pf_t contributions = 0;
	contributions += PRmloop10.get(i,j,k+1,l-1)*expap_penalty*beta2P(l,k);
    contributions += PRmloop01.get(i,j,k+1,l-1)*expap_penalty*beta2P(l,k);

	return contributions;
}

pf_t W_final_pf::get_PMiloop(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	if (!(i <= j && j < k-1 && k <= l)){
		return 0;
	}
	assert(!(i<=0 || l> n));

	const int ptype_closing = pair[S_[j]][S_[k]];
	if(ptype_closing == 0) return 0;

	pf_t contributions = 0;
	contributions += PM.get(i,j-1,k+1,l)*get_e_stP(j-1,k+1);

	for(cand_pos_t d= j-1; d>std::max(i,j-MAXLOOP); --d){
		for (cand_pos_t dp=k+1; dp <std::min(l,k+MAXLOOP); ++dp) {
			cand_pos_t u1 = (j)-d-1; // check these
			cand_pos_t u2 = dp-(k)-1;
			contributions += get_e_intP(d,j,k,dp)*PM.get(i,d,dp,l)*scale[u1 + u2 + 2];
		}
	}
	return contributions;
}

pf_t W_final_pf::get_PMmloop(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	if (!(i <= j && j < k-1 && k <= l)){
		return 0;
	}
	assert(!(i<=0 || l> n));
    pf_t contributions = 0;
	contributions += PMmloop10.get(i,j-1,k+1,l)*expap_penalty*beta2P(j,k);
    contributions += PMmloop01.get(i,j-1,k+1,l)*expap_penalty*beta2P(j,k);

	return contributions;
}

pf_t W_final_pf::get_POiloop(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	if (!(i <= j && j < k-1 && k <= l)){
		return 0;
	}
	assert(!(i<=0 || l> n));

	const int ptype_closing = pair[S_[i]][S_[l]];
	if(ptype_closing == 0) return 0;

	pf_t contributions = 0;
	contributions += PO.get(i+1,j,k,l-1)*get_e_stP(i,l);

	for(cand_pos_t d= i+1; d<std::min(j,i+MAXLOOP); ++d){
		for (cand_pos_t dp=l-1; dp >std::max(l-MAXLOOP,k); --dp) {
			cand_pos_t u1 = d-(i)-1; // check these
			cand_pos_t u2 = (l)-dp-1;
			contributions += get_e_intP(i,d,dp,l)*PO.get(d,j,dp,k)*scale[u1 + u2 + 2];
		}
	}
	return contributions;
}

pf_t W_final_pf::get_POmloop(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	if (!(i <= j && j < k-1 && k <= l)){
		return 0;
	}
	assert(!(i<=0 || l> n));
    pf_t contributions = 0;
	contributions += POmloop10.get(i+1,j,k,l-1)*expap_penalty*beta2P(l,i);
    contributions += POmloop01.get(i+1,j,k,l-1)*expap_penalty*beta2P(l,i);

	return contributions;
}

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


// penalty for closing pair i.l or l.i of an ordinary multiloop
inline pf_t W_final_pf::beta2(cand_pos_t i, cand_pos_t l){
	return expb_penalty;
}

// penalty for closing pair i.l or l.i of a multiloop that spans a band
inline pf_t W_final_pf::beta2P(cand_pos_t i, cand_pos_t l){
	return expbp_penalty;
}