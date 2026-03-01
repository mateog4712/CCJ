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
    cand_pos_t total_length = ((n + 1) * (n + 2)) / 2;
    index[1] = 0;
    for (cand_pos_t i = 2; i <= n; i++)
        index[i] = index[i - 1] + (n + 1) - i + 1;
    // Allocate space
    V.resize(total_length, 0);
    VM.resize(total_length, 0);
    WM.resize(total_length, 0);
    WMv.resize(total_length, 0);
    WMp.resize(total_length, 0);

    // PK
    WBP.resize(total_length,0);
	WPP.resize(total_length,0);
	P.resize(total_length,0);

    std::vector<energy_t> row;
	row.resize(total_length,0);

    // Mateo - Mirroring ModifiedCCJ, changing PL from n^2xn^2 to nxnxn^2
	PL.resize(n+1);
	PO.resize(n+1);
	PfromL.resize(n+1);
	PfromO.resize(n+1);
	PLmloop01.resize(n+1);
	PLmloop10.resize(n+1);
	POmloop01.resize(n+1);
	POmloop10.resize(n+1);
	for(cand_pos_t i = 0; i <= n; ++i) {
		for(cand_pos_t j = 0; j <= n; ++j) {
			PL[i].push_back(row);
			PO[i].push_back(row);
			PfromL[i].push_back(row);
			PfromO[i].push_back(row);

			PLmloop01[i].push_back(row);
			PLmloop10[i].push_back(row);
			
			POmloop01[i].push_back(row);
			POmloop10[i].push_back(row);
		}
	}

	for(cand_pos_t i = 0; i < total_length; ++i) {
		PK.push_back(row);
		PR.push_back(row);
		PM.push_back(row);

		PfromR.push_back(row);
		PfromM.push_back(row);

		PLiloop.push_back(row);
		PLmloop00.push_back(row);

		PRiloop.push_back(row);
		PRmloop00.push_back(row);
		PRmloop01.push_back(row);
		PRmloop10.push_back(row);

		PMiloop.push_back(row);
		PMmloop00.push_back(row);
		PMmloop01.push_back(row);
		PMmloop10.push_back(row);

		POiloop.push_back(row);
		POmloop00.push_back(row);
	}

    rescale_pk_globals();
    exp_params_rescale(MFE_energy);
    W.resize(n + 1, scale[1]);
    // WI.resize(total_length, scale[1]);

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
			energy_t acc = (k>1) ? W[k-1]: 0;
			contributions += acc*get_energy(k, j)*exp_Extloop(k, j);
			contributions += acc*get_P(k, j)*expPS_penalty;
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
    WMp_contributions += (get_P(i, j) * expPSM_penalty * expb_penalty);
    WMv_contributions += (get_energy_WMv(i, j - 1) * expMLbase[1]);
    WMp_contributions += (get_energy_WMp(i, j - 1) * expMLbase[1]);

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
        pf_t qbt2 = get_P(k, j) * expPSM_penalty * expb_penalty;
        contributions += (static_cast<pf_t>(expMLbase[k - i]) * qbt1);
        contributions += (static_cast<pf_t>(expMLbase[k - i]) * qbt2);
        contributions += (get_energy_WM(i, k - 1) * qbt1);
        contributions += (get_energy_WM(i, k - 1) * qbt2);
    }
    contributions += WM[ijminus1] * expMLbase[1];
    WM[ij] = contributions;
}

pf_t W_final_pf::compute_energy_VM(cand_pos_t i, cand_pos_t j) {
    pf_t contributions = 0;
    cand_pos_t ij = index[(i)] + (j) - (i);
    for (cand_pos_t k = i + 1; k <= j - TURN - 1; ++k) {
        contributions += (get_energy_WM(i + 1, k - 1) * get_energy_WMv(k, j - 1) * exp_Mbloop(i, j) * exp_params_->expMLclosing);
        contributions += (get_energy_WM(i + 1, k - 1) * get_energy_WMp(k, j - 1) * exp_Mbloop(i, j) * exp_params_->expMLclosing);
        contributions += (expMLbase[k - i - 1] * get_energy_WMp(k, j - 1) * exp_Mbloop(i, j) * exp_params_->expMLclosing);
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

			compute_PLiloop(i,j,k,l);

			compute_PLmloop00(i,j,k,l);

			compute_PLmloop01(i,j,k,l);

			compute_PLmloop10(i,j,k,l);

			compute_PRiloop(i,j,k,l);

			compute_PRmloop00(i,j,k,l);

			compute_PRmloop01(i,j,k,l);

			compute_PRmloop10(i,j,k,l);

			compute_PMiloop(i,j,k,l);

			compute_PMmloop00(i,j,k,l);

			compute_PMmloop01(i,j,k,l);

			compute_PMmloop10(i,j,k,l);

			compute_POiloop(i,j,k,l);

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

	cand_pos_t il = index[i]+l-i;
	for(cand_pos_t d=i; d< l; ++d){
		for(cand_pos_t e = d+1; e<= l; ++e){
			energy_t common = get_WB(i,d-1)*expcp_pen[l-e]*expPPS_penalty;
			contributions += get_energy(d,e)*beta2P(e,d)*common;
			contributions += get_P(d,e)*expPSM_penalty*common;
		}
	}
	WBP[il] = contributions;
}

void W_final_pf::compute_WPP(cand_pos_t i, cand_pos_t l){
	pf_t contributions = 0;

	cand_pos_t il = index[i]+l-i;
	for(cand_pos_t d=i; d<l; ++d){
		for(cand_pos_t e = d+1; e<= l; ++e){
			energy_t common = get_WP(i,d-1)*expPUP_pen[l-e]*expPPS_penalty;
			contributions += get_energy(d,e)*gamma2(e,d)*common;
			contributions += get_P(d,e)*expPSP_penalty*common;
		}
	}
	WPP[il] = contributions;
}

void W_final_pf::compute_P(cand_pos_t i, cand_pos_t l){
	cand_pos_t il = index[i]+l-i;
    pf_t contributions = 0;
	for(cand_pos_t j=i; j<l; ++j){
		for (cand_pos_t d=j+1; d<l; ++d){
			for (cand_pos_t k=d+1; k<l; ++k){
				contributions += get_PK(i,j,d+1,k)*get_PK(j+1,d,k+1,l);
			}
		}
	}
	P[il]=contributions;
}

void W_final_pf::compute_PK(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	pf_t contributions = 0;

	cand_pos_t ij = index[i]+j-i;
	cand_pos_t kl = index[k]+l-k;

	for(cand_pos_t d=i+1; d< j; ++d){
		contributions += get_PK(i,d,k,l)*get_WP(d+1,j);  //12G1
	}

	for(cand_pos_t d=k+1; d< l; ++d){
		contributions += get_PK(i,j,d,l)*get_WP(k,d-1); //1G21
	}

	contributions += get_PL(i,j,k,l)*gamma2(j,i)*expPB_penalty;
	contributions += get_PM(i,j,k,l)*gamma2(j,k)*expPB_penalty;
	contributions += get_PR(i,j,k,l)*gamma2(l,k)*expPB_penalty;
	contributions += get_PO(i,j,k,l)*gamma2(l,i)*expPB_penalty;
	
    PK[ij][kl]=contributions;
}

void W_final_pf::compute_PL(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	pf_t contributions = 0;

	cand_pos_t kl = index[k]+l-k;
	const int ptype_closing = pair[S_[i]][S_[j]];

	if (ptype_closing>0){

		contributions += get_PLiloop(i,j,k,l);

		contributions += get_PLmloop(i,j,k,l)*expbp_penalty;

		if (j>=(i+TURN+1)){
			contributions += get_PfromL(i+1,j-1,k,l)*gamma2(j,i);
		}
	}
	PL[i][j][kl]=contributions;
}

void W_final_pf::compute_PR(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	pf_t contributions = 0;

	cand_pos_t ij = index[i]+j-i;
	cand_pos_t kl = index[k]+l-k;

	const int ptype_closing = pair[S_[k]][S_[l]];

	if (ptype_closing>0){
		contributions += get_PRiloop(i,j,k,l);

		contributions += get_PRmloop(i,j,k,l)*expbp_penalty;

		if (l>=(k+TURN+1)){
			contributions += get_PfromR(i,j,k+1,l-1)*gamma2(l,k);
		}
	}
	PR[ij][kl]=contributions;
}

void W_final_pf::compute_PM(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	pf_t contributions = 0;

	cand_pos_t ij = index[i]+j-i;
	cand_pos_t kl = index[k]+l-k;

	const int ptype_closing = pair[S_[j]][S_[k]];
	if (ptype_closing>0){
		contributions += get_PMiloop(i,j,k,l);

		contributions += get_PMmloop(i,j,k,l)*expbp_penalty;

		if (k>=(j+TURN-1)){
			contributions += get_PfromM(i,j-1,k+1,l)*gamma2(j,k);
		}

		if(i==j && k==l){
			contributions += gamma2(i,l);
		}
	}
	PM[ij][kl]=contributions;
}

void W_final_pf::compute_PO(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	pf_t contributions = 0;

	cand_pos_t kl = index[k]+l-k;

	const int ptype_closing = pair[S_[i]][S_[l]];

	if (ptype_closing>0){
		contributions += get_POiloop(i,j,k,l);

		contributions += get_POmloop(i,j,k,l)*expbp_penalty;

		// Hosna, July 11, 2014
		// To avoid addition of close base pairs we check for the following here
		if (l>=(i+TURN+1)){
			contributions += get_PfromO(i+1,j,k,l-1)*gamma2(l,i);
		}
	}
	PO[i][j][kl]=contributions;
}

void W_final_pf::compute_PfromL(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	pf_t contributions = 0;

	cand_pos_t kl = index[k]+l-k;

	for(cand_pos_t d=i+1; d< j; ++d){
		contributions += get_PfromL(d,j,k,l)*get_WP(i,d-1);
		contributions += get_PfromL(i,d,k,l)*get_WP(d+1,j);
	}

	contributions += get_PR(i,j,k,l)*gamma2(l,k)*expPB_penalty;

	contributions += get_PM(i,j,k,l)*gamma2(j,k)*expPB_penalty;

	contributions += get_PO(i,j,k,l)*gamma2(l,i)*expPB_penalty;

	PfromL[i][j][kl]=contributions;
}

void W_final_pf::compute_PfromR(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	pf_t contributions = 0;

	cand_pos_t ij = index[i]+j-i;
	cand_pos_t kl = index[k]+l-k;
	for(cand_pos_t d=k+1; d< l; ++d){
		contributions += get_PfromR(i,j,d,l)*get_WP(k,d-1);
		contributions += get_PfromR(i,j,k,d)+get_WP(d+1,l);
	}

	contributions += get_PM(i,j,k,l)*gamma2(j,k)*expPB_penalty;

	contributions += get_PO(i,j,k,l)*gamma2(l,i)*expPB_penalty;

	PfromR[ij][kl]=contributions;
}

void W_final_pf::compute_PfromM(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	pf_t contributions = 0;

	cand_pos_t ij = index[i]+j-i;
	cand_pos_t kl = index[k]+l-k;

	for(cand_pos_t d=i+1; d<j; ++d){
		contributions += get_PfromM(i,d,k,l)*get_WP(d+1,j);
	}
	for(cand_pos_t d=k+1; d<l; ++d){
		contributions += get_PfromM(i,j,d,l)*get_WP(k,d-1);
	}

	contributions += get_PL(i,j,k,l)*gamma2(j,i)*expPB_penalty;

	contributions += get_PR(i,j,k,l)*gamma2(l,k)*expPB_penalty;

	PfromM[ij][kl]=contributions;
}

void W_final_pf::compute_PfromO(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	pf_t contributions = 0;

	cand_pos_t kl = index[k]+l-k;

	for(cand_pos_t d=i+1; d< j; ++d){
		contributions += get_PfromO(d,j,k,l)*get_WP(i,d-1);
	}
	for(cand_pos_t d=k+1; d<l; ++d){
		contributions += get_PfromO(i,j,k,d)*get_WP(d+1,l);
	}

	contributions += get_PL(i,j,k,l)*gamma2(j,i)*expPB_penalty;

	contributions += get_PR(i,j,k,l)*gamma2(l,k)*expPB_penalty;

    PfromO[i][j][kl]=contributions;
}

void W_final_pf::compute_PLiloop(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	pf_t contributions = 0;

	cand_pos_t ij = index[i]+j-i;
	cand_pos_t kl = index[k]+l-k;
	int ptype_closing = pair[S_[i]][S_[j]];
	if (ptype_closing>0){
		contributions += get_PL(i+1,j-1,k,l)*get_e_stP(i,j);

		for(cand_pos_t d= i+1; d<std::min(j,i+MAXLOOP); ++d){
        	for(cand_pos_t dp = j-1; dp > std::max(d+TURN,j-MAXLOOP); --dp){
            	cand_pos_t u1 = d - (i+1) - 1; // check these
                cand_pos_t u2 = (j-1) - dp - 1;
                contributions += get_e_intP(i,d,dp,j)*get_PL(d,dp,k,l)*scale[u1 + u2 + 2];
        	}
    	}
	}
	PLiloop[ij][kl]=contributions;
}

void W_final_pf::compute_PLmloop00(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	pf_t contributions = 0;

	cand_pos_t ij = index[i]+j-i;
	cand_pos_t kl = index[k]+l-k;

	contributions += get_PL(i,j,k,l)*beta2P(j,i);
    for(cand_pos_t d = i; d<=j; ++d){
        if (d>i){
            contributions += get_WB(i,d-1)*get_PLmloop00(d,j,k,l);
        }
        if(d<j){
            contributions += get_PLmloop00(i,d,k,l)*get_WB(d+1,j);
        }

    }
	PLmloop00[ij][kl]=contributions;	
}

void W_final_pf::compute_PLmloop01(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	pf_t contributions = 0;

	cand_pos_t kl = index[k]+l-k;
	for(cand_pos_t d = i; d < j; ++d){
        contributions += get_PLmloop00(i,d,k,l)*get_WBP(d+1,j);
    }
	PLmloop01[i][j][kl] = contributions;

}

void W_final_pf::compute_PLmloop10(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	pf_t contributions = 0;

	cand_pos_t kl = index[k]+l-k;
	for(cand_pos_t d = i+1; d <= j; ++d){
        contributions += get_WBP(i,d-1)*get_PLmloop00(d,j,k,l);
        if(d<j){
            contributions += get_PLmloop10(i,d,k,l)*get_WB(d+1,j);
        }
    }
	PLmloop10[i][j][kl] = contributions;
}

void W_final_pf::compute_PRiloop(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	pf_t contributions = 0;

	cand_pos_t ij = index[i]+j-i;
	cand_pos_t kl = index[k]+l-k;
	int ptype_closing = pair[S_[k]][S_[l]];
	if (ptype_closing>0){
		contributions += get_PR(i,j,k+1,l-1)*get_e_stP(k,l);

		for(cand_pos_t d= k+1; d<std::min(l,k+MAXLOOP); ++d){
        	for(cand_pos_t dp=l-1; dp > std::max(d+TURN,l-MAXLOOP); --dp){
                cand_pos_t u1 = d - (k+1) - 1; // check these
                cand_pos_t u2 = (l-1) - dp - 1;
            	contributions += get_e_intP(k,d,dp,l)*get_PR(i,j,d,dp)*scale[u1 + u2 + 2];
        	}
    	}
	}
	PRiloop[ij][kl]=contributions;
}

void W_final_pf::compute_PRmloop00(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	pf_t contributions = 0;
	cand_pos_t ij = index[i]+j-i;
	cand_pos_t kl = index[k]+l-k;

	contributions += get_PR(i,j,k,l)*beta2P(l,k);
	for(cand_pos_t d=k; d<=l; ++d){
        if(d>k){
            contributions += get_WB(k,d-1)*get_PRmloop00(i,j,d,l);
        }
        if (d<l){
            contributions += get_PRmloop00(i,j,k,d)*get_WB(d+1,l);
        }
    }
	PRmloop00[ij][kl]=contributions;
}


void W_final_pf::compute_PRmloop01(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	pf_t contributions = 0;

	cand_pos_t ij = index[i]+j-i;
	cand_pos_t kl = index[k]+l-k;

	contributions += get_PRmloop01(i,j,k,l-1)*expcp_pen[1];
    for(cand_pos_t d = k; d < l; d++){
        contributions += get_PRmloop00(i,j,k,d)*get_WBP(d+1,l);
    }
	PRmloop01[ij][kl] = contributions;
}

void W_final_pf::compute_PRmloop10(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	pf_t contributions = 0;

	cand_pos_t ij = index[i]+j-i;
	cand_pos_t kl = index[k]+l-k;

	contributions += get_PRmloop10(i,j,k+1,l)*expcp_pen[1];
    for(cand_pos_t d = k+1; d <= l; ++d){
        contributions += get_WBP(k,d-1)*get_PRmloop00(i,j,d,l);
    }
	PRmloop10[ij][kl] = contributions;
}

void W_final_pf::compute_PMiloop(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	pf_t contributions = 0;

	cand_pos_t ij = index[i]+j-i;
	cand_pos_t kl = index[k]+l-k;
	int ptype_closing = pair[S_[j]][S_[k]];

	if (ptype_closing>0){
		contributions += get_PM(i,j-1,k+1,l)*get_e_stP(j-1,k+1);

		for(cand_pos_t d= j-1; d>std::max(i,j-MAXLOOP); --d){
        	for (cand_pos_t dp=k+1; dp <std::min(l,k+MAXLOOP); ++dp) {
                cand_pos_t u1 = d - (j-1) - 1; // check these
                cand_pos_t u2 = (k+1) - dp - 1;
            	contributions += get_e_intP(d,j,k,dp)*get_PM(i,d,dp,l)*scale[u1 + u2 + 2];
        	}
    	}
	}
	PMiloop[ij][kl]=contributions;
}

void W_final_pf::compute_PMmloop00(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	pf_t contributions = 0;

	cand_pos_t ij = index[i]+j-i;
	cand_pos_t kl = index[k]+l-k;

	contributions += get_PM(i,j,k,l)*beta2P(j,k);
    for(cand_pos_t d=i; d<j; ++d){
        contributions=get_PMmloop00(i,d,k,l)*get_WB(d+1,j);
    }
    for(cand_pos_t d=k+1; d<=l; ++d){
        contributions=get_PMmloop00(i,j,d,l)*get_WB(k,d-1);
    }
	PMmloop00[ij][kl]=contributions;
}


void W_final_pf::compute_PMmloop01(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	pf_t contributions = 0;

	cand_pos_t ij = index[i]+j-i;
	cand_pos_t kl = index[k]+l-k;
	contributions += get_PMmloop01(i,j,k+1,l)+expcp_pen[1];
	for(cand_pos_t d = k; d < l; ++d){
        contributions += get_POmloop00(i,j,k,d)*get_WBP(d+1,l);
    }
	PMmloop01[ij][kl] = contributions;
}

void W_final_pf::compute_PMmloop10(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	pf_t contributions = 0;

	cand_pos_t ij = index[i]+j-i;
	cand_pos_t kl = index[k]+l-k;

	contributions += get_PMmloop10(i,j-1,k,l)*expcp_pen[1];
	for(cand_pos_t d=i+1; d<=j;++d){
        contributions += get_WBP(i,d-1)*get_POmloop00(d,j,k,l);
    }
    for(cand_pos_t d = k+1; d < l; ++d){
        contributions += get_POmloop10(i,j,k,d)*get_WB(d+1,l);
    }
	PMmloop10[ij][kl] = contributions;
}

void W_final_pf::compute_POiloop(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	pf_t contributions = 0;

	cand_pos_t ij = index[i]+j-i;
	cand_pos_t kl = index[k]+l-k;
	int ptype_closing = pair[S_[i]][S_[l]];

	if (ptype_closing>0){
		contributions += get_PO(i+1,j,k,l-1)*get_e_stP(i,l);

		for(cand_pos_t d= i+1; d<std::min(j,i+MAXLOOP); ++d){
        	for (cand_pos_t dp=l-1; dp >std::max(l-MAXLOOP,k); --dp) {
                cand_pos_t u1 = d - (i+1) - 1; // check these
                cand_pos_t u2 = (l-1) - dp - 1;
            	contributions += get_e_intP(i,d,dp,l)*get_PO(d,j,dp,k)*scale[u1 + u2 + 2];
        	}
    	}
	}
	POiloop[ij][kl]=contributions;
}

void W_final_pf::compute_POmloop00(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	pf_t contributions = 0;

	cand_pos_t ij = index[i]+j-i;
	cand_pos_t kl = index[k]+l-k;

	contributions += get_PO(i,j,k,l)*beta2P(l,i);
    for(cand_pos_t d=i+1; d<=j; ++d){
        contributions += get_WB(i,d-1)*get_POmloop00(d,j,k,l);
    }
    for(cand_pos_t d=k; d<l; ++d){
        contributions = get_POmloop00(i,j,k,d)*get_WB(d+1,l);
    }
	POmloop00[ij][kl]=contributions;
}


void W_final_pf::compute_POmloop01(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	pf_t contributions = 0;

	cand_pos_t kl = index[k]+l-k;
	for(cand_pos_t d = k; d < l; ++d){
        contributions += get_POmloop00(i,j,k,d)*get_WBP(d+1,l);
    }
	POmloop01[i][j][kl] = contributions;
}

void W_final_pf::compute_POmloop10(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	pf_t contributions = 0;

	cand_pos_t kl = index[k]+l-k;
	for(cand_pos_t d=i+1; d<=j;++d){
        contributions += get_WBP(i,d-1)*get_POmloop00(d,j,k,l);
    }
    for(cand_pos_t d = k+1; d < l; ++d){
        contributions += get_POmloop10(i,j,k,d) + get_WB(d+1,l);
    }

	POmloop10[i][j][kl] = contributions;
}

energy_t W_final_pf::get_WB(cand_pos_t i, cand_pos_t j){
	if (i<=0 || j<=0 || i>n || j>n){
		return 0;
	}
	if (i>j) return 1;
	return (std::min(cp_penalty*(j-i+1),get_WBP(i,j)));
}

energy_t W_final_pf::get_WBP(cand_pos_t i, cand_pos_t j){
	if (i > j || i<=0 || j<=0 || i>n || j>n){
		return 0;
	}
	cand_pos_t ij = index[i]+j-i;
	return WBP[ij];
}

energy_t W_final_pf::get_WP(cand_pos_t i, cand_pos_t j){
	if (i<=0 || j<=0 || i>n || j>n){
		return 0;
	}
	if (i>j) return 1; // needed as this will happen 
	return (std::min(PUP_penalty*(j-i+1),get_WPP(i,j)));
}

energy_t W_final_pf::get_WPP(cand_pos_t i, cand_pos_t j){
	if (i > j  || i<=0 || j<=0 || i>n || j>n){
		return 0;
	}
	cand_pos_t ij = index[i]+j-i;
	return WPP[ij];
}

energy_t W_final_pf::get_P(cand_pos_t i, cand_pos_t j){
	if (i >= j  || i<=0 || j<=0 || i>n || j>n){
		return 0;
	}
	cand_pos_t ij = index[i]+j-i;
	return P[ij];
}

energy_t W_final_pf::get_PK(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	if (!(i <= j && j < k-1 && k <= l)){
		return 0;
	}
	assert(!(i<=0 || l> n));
	
	cand_pos_t ij = index[i]+j-i;
	cand_pos_t kl = index[k]+l-k;

	return PK[ij][kl];
}

energy_t W_final_pf::get_PL(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	if (!(i <= j && j < k-1 && k <= l)){
		return 0;
	}
	assert(!(i<=0 || l> n));

	cand_pos_t kl = index[k]+l-k;

	return PL[i][j][kl];
}

energy_t W_final_pf::get_PR(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	if (!(i <= j && j < k-1 && k <= l)){
		return 0;
	}
	assert(!(i<=0 || l> n));

	cand_pos_t ij = index[i]+j-i;
	cand_pos_t kl = index[k]+l-k;

	return PR[ij][kl];
}

energy_t W_final_pf::get_PM(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	if (!(i <= j && j < k-1 && k <= l)){
		return 0;
	}
	assert(!(i<=0 || l> n));

	int ptype_closing = pair[S_[j]][S_[k]];
	if (ptype_closing == 0){
		return 0;
	}
	cand_pos_t ij = index[i]+j-i;
	cand_pos_t kl = index[k]+l-k;
	if (i ==j && k ==l){
		return gamma2(i,l);
	}

	return PM[ij][kl];
}

energy_t W_final_pf::get_PO(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	if (!(i <= j && j < k-1 && k <= l)){
		return 0;
	}
	assert(!(i<=0 || l> n));

	cand_pos_t kl = index[k]+l-k;

	return PO[i][j][kl];
}

energy_t W_final_pf::get_PfromL(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	if (!(i <= j && j < k-1 && k <= l)){
		return 0;
	}
	assert(!(i<=0 || l> n));

	if (i==j && k==l){
		const int ptype_closing = pair[S_[i]][S_[l]];
		if(ptype_closing == 0) return 0;
		else return gamma2(j,k)*gamma2(k,j);
	}
	cand_pos_t kl = index[k]+l-k;

	return PfromL[i][j][kl];
}

energy_t W_final_pf::get_PfromR(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	if (!(i <= j && j < k-1 && k <= l)){
		return 0;
	}
	assert(!(i<=0 || l> n));

	if (i==j && k==l){
		const int ptype_closing = pair[S_[i]][S_[l]];
		if(ptype_closing == 0) return 0;
		else return gamma2(j,k)*gamma2(k,j);
	}
	cand_pos_t ij = index[i]+j-i;
	cand_pos_t kl = index[k]+l-k;

	return PfromR[ij][kl];
}

energy_t W_final_pf::get_PfromM(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	if (!(i <= j && j < k-1 && k <= l)){
		return 0;
	}
	assert(!(i<=0 || l> n));

	if (i==j && k==l){
		const int ptype_closing = pair[S_[i]][S_[l]];
		if(ptype_closing == 0) return 0;
		else return 1;
	}
	cand_pos_t ij = index[i]+j-i;
	cand_pos_t kl = index[k]+l-k;

	return PfromM[ij][kl];
}

energy_t W_final_pf::get_PfromO(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	if (!(i <= j && j < k-1 && k <= l)){
		return 0;
	}
	assert(!(i<=0 || l> n));

	if (i==j && k==l){
		const int ptype_closing = pair[S_[i]][S_[l]];
		if(ptype_closing == 0) return 0;
		else return 1;
	}
	cand_pos_t kl = index[k]+l-k;

	return PfromO[i][j][kl];
}

energy_t W_final_pf::get_PLiloop(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	if (!(i <= j && j < k-1 && k <= l)){
		return 0;
	}
	assert(!(i<=0 || l> n));

	const int ptype_closing = pair[S_[i]][S_[j]];
	if(ptype_closing == 0) return 0;

	cand_pos_t ij = index[i]+j-i;
	cand_pos_t kl = index[k]+l-k;

	return PLiloop[ij][kl];
}

energy_t W_final_pf::get_PLmloop(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	if (!(i <= j && j < k-1 && k <= l)){
		return 0;
	}
	assert(!(i<=0 || l> n));
    pf_t contributions = 0;
	contributions += get_PLmloop10(i+1,j-1,k,l)*expap_penalty*beta2P(j,i);
    contributions += get_PLmloop01(i+1,j-1,k,l)*ap_penalty*beta2P(j,i);

	return contributions;
}

energy_t W_final_pf::get_PLmloop00(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	if (!(i <= j && j < k-1 && k <= l)){
		return 0;
	}
	assert(!(i<=0 || l> n));

	cand_pos_t ij = index[i]+j-i;
	cand_pos_t kl = index[k]+l-k;

	return PLmloop00[ij][kl];
}

energy_t W_final_pf::get_PLmloop01(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	if (!(i <= j && j < k-1 && k <= l)){
		return 0;
	}
	assert(!(i<=0 || l> n));

	cand_pos_t kl = index[k]+l-k;

	return PLmloop01[i][j][kl];
}

energy_t W_final_pf::get_PLmloop10(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	if (!(i <= j && j < k-1 && k <= l)){
		return 0;
	}
	assert(!(i<=0 || l> n));

	cand_pos_t kl = index[k]+l-k;

	return PLmloop10[i][j][kl];
}

energy_t W_final_pf::get_PRiloop(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	if (!(i <= j && j < k-1 && k <= l)){
		return 0;
	}
	assert(!(i<=0 || l> n));

	const int ptype_closing = pair[S_[k]][S_[l]];
	if(ptype_closing == 0) return 0;

	cand_pos_t ij = index[i]+j-i;
	cand_pos_t kl = index[k]+l-k;

	return PRiloop[ij][kl];
}

energy_t W_final_pf::get_PRmloop(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	if (!(i <= j && j < k-1 && k <= l)){
		return 0;
	}
	assert(!(i<=0 || l> n));
    pf_t contributions = 0;
	contributions += get_PRmloop10(i,j,k+1,l-1)*expap_penalty*beta2P(l,k);
    contributions += get_PRmloop01(i,j,k+1,l-1)*expap_penalty*beta2P(l,k);

	return contributions;
}

energy_t W_final_pf::get_PRmloop00(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	if (!(i <= j && j < k-1 && k <= l)){
		return 0;
	}
	assert(!(i<=0 || l> n));

	cand_pos_t ij = index[i]+j-i;
	cand_pos_t kl = index[k]+l-k;

	return PRmloop00[ij][kl];
}

energy_t W_final_pf::get_PRmloop01(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	if (!(i <= j && j < k-1 && k <= l)){
		return 0;
	}
	assert(!(i<=0 || l> n));

	cand_pos_t ij = index[i]+j-i;
	cand_pos_t kl = index[k]+l-k;

	return PRmloop01[ij][kl];
}

energy_t W_final_pf::get_PRmloop10(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	if (!(i <= j && j < k-1 && k <= l)){
		return 0;
	}
	assert(!(i<=0 || l> n));

	cand_pos_t ij = index[i]+j-i;
	cand_pos_t kl = index[k]+l-k;

	return PRmloop10[ij][kl];
}

energy_t W_final_pf::get_PMiloop(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	if (!(i <= j && j < k-1 && k <= l)){
		return 0;
	}
	assert(!(i<=0 || l> n));

	const int ptype_closing = pair[S_[j]][S_[k]];
	if(ptype_closing == 0) return INF;

	cand_pos_t ij = index[i]+j-i;
	cand_pos_t kl = index[k]+l-k;

	return PMiloop[ij][kl];
}

energy_t W_final_pf::get_PMmloop(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	if (!(i <= j && j < k-1 && k <= l)){
		return 0;
	}
	assert(!(i<=0 || l> n));
    pf_t contributions = 0;
	contributions += get_PMmloop10(i,j-1,k+1,l)*expap_penalty*beta2P(j,k);
    contributions += get_PMmloop01(i,j-1,k+1,l)*expap_penalty*beta2P(j,k);

	return contributions;
}

energy_t W_final_pf::get_PMmloop00(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	if (!(i <= j && j < k-1 && k <= l)){
		return 0;
	}
	assert(!(i<=0 || l> n));

	cand_pos_t ij = index[i]+j-i;
	cand_pos_t kl = index[k]+l-k;

	return PMmloop00[ij][kl];
}

energy_t W_final_pf::get_PMmloop01(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	if (!(i <= j && j < k-1 && k <= l)){
		return 0;
	}
	assert(!(i<=0 || l> n));

	cand_pos_t ij = index[i]+j-i;
	cand_pos_t kl = index[k]+l-k;

	return PMmloop01[ij][kl];
}

energy_t W_final_pf::get_PMmloop10(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	if (!(i <= j && j < k-1 && k <= l)){
		return 0;
	}
	assert(!(i<=0 || l> n));

	cand_pos_t ij = index[i]+j-i;
	cand_pos_t kl = index[k]+l-k;

	return PMmloop10[ij][kl];
}

energy_t W_final_pf::get_POiloop(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	if (!(i <= j && j < k-1 && k <= l)){
		return 0;
	}
	assert(!(i<=0 || l> n));

	const int ptype_closing = pair[S_[i]][S_[l]];
	if(ptype_closing == 0) return 0;

	cand_pos_t ij = index[i]+j-i;
	cand_pos_t kl = index[k]+l-k;

	return POiloop[ij][kl];
}

energy_t W_final_pf::get_POmloop(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	if (!(i <= j && j < k-1 && k <= l)){
		return 0;
	}
	assert(!(i<=0 || l> n));
    pf_t contributions = 0;
	contributions += get_POmloop10(i+1,j,k,l-1)*expap_penalty*beta2P(l,i);
    contributions += get_POmloop01(i+1,j,k,l-1)*expap_penalty*beta2P(l,i);

	return contributions;
}

energy_t W_final_pf::get_POmloop00(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	if (!(i <= j && j < k-1 && k <= l)){
		return 0;
	}
	assert(!(i<=0 || l> n));

	cand_pos_t ij = index[i]+j-i;
	cand_pos_t kl = index[k]+l-k;

	return POmloop00[ij][kl];
}

energy_t W_final_pf::get_POmloop01(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	if (!(i <= j && j < k-1 && k <= l)){
		return 0;
	}
	assert(!(i<=0 || l> n));

	cand_pos_t kl = index[k]+l-k;

	return POmloop01[i][j][kl];
}

energy_t W_final_pf::get_POmloop10(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	if (!(i <= j && j < k-1 && k <= l)){
		return 0;
	}
	assert(!(i<=0 || l> n));

	cand_pos_t kl = index[k]+l-k;

	return POmloop10[i][j][kl];
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
energy_t W_final_pf::beta2(cand_pos_t i, cand_pos_t l){
	return expb_penalty;
}

// penalty for closing pair i.l or l.i of a multiloop that spans a band
energy_t W_final_pf::beta2P(cand_pos_t i, cand_pos_t l){
	return expbp_penalty;
}

// penalty for closing pair i.l or l.i of a pseudoloop
energy_t W_final_pf::gamma2(cand_pos_t i, cand_pos_t l){
	return 1;
}