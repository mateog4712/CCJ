#include "pseudo_loop.hh"
#include "h_externs.hh"
#include <stdio.h>
#include <string>
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <algorithm>

#define debug 0

pseudo_loop::pseudo_loop(std::string seq, s_energy_matrix *V, short *S, short *S1, vrna_param_t *params)
{
	this->seq = seq;
	this->V = V;
	S_ = S;
	S1_ = S1;
	params_ = params;
	make_pair_matrix();
    allocate_space();
}

void pseudo_loop::allocate_space()
{
    n = seq.length();

    index.resize(n+1);
    cand_pos_t total_length = ((n+1) *(n+2))/2;
    index[1] = 0;
    for (cand_pos_t i=2; i <= n; i++)
        index[i] = index[i-1]+(n+1)-i+1;

    WBP.resize(total_length,INF);
	WPP.resize(total_length,INF);
	P.resize(total_length,INF);

	std::vector<energy_t> row;
	row.resize(total_length,INF);

	std::vector<energy_t> row_5;
	row_5.resize(MAXLOOP+1,INF);

	for(cand_pos_t i = 0; i < total_length; ++i) {
		PK.push_back(row);
		PL.push_back(row);
		PR.push_back(row);
		PM.push_back(row);
		PO.push_back(row);

		PfromL.push_back(row);
		PfromR.push_back(row);
		PfromM.push_back(row);
		PfromO.push_back(row);

		PLiloop.push_back(row);
		PLmloop.push_back(row);
		PLmloop0.push_back(row);
		PLmloop1.push_back(row);

		PRiloop.push_back(row);
		PRmloop.push_back(row);
		PRmloop0.push_back(row);
		PRmloop1.push_back(row);

		PMiloop.push_back(row);
		PMmloop.push_back(row);
		PMmloop0.push_back(row);
		PMmloop1.push_back(row);

		POiloop.push_back(row);
		POmloop.push_back(row);
		POmloop0.push_back(row);
		POmloop1.push_back(row);
	}
}

pseudo_loop::~pseudo_loop()
{
}

void pseudo_loop::compute_energies(cand_pos_t i, cand_pos_t l)
{

	// 1) compute all energies over region [i,l]
	compute_P(i,l);

	compute_WBP(i,l);

	compute_WPP(i,l);

	//2) compute all energies over gapped region [i,j]U[k,l]
	for(cand_pos_t j = i; j<l; ++j){
		// Hosna, July 8, 2014
		// in original recurrences we have j< k-1, so I am changing k=j+1 to k=j+2
		for(cand_pos_t k = l; k>=j+2; --k){

			compute_PLiloop(i,j,k,l);

			compute_PLmloop(i,j,k,l);

			compute_PLmloop0(i,j,k,l);

			compute_PLmloop1(i,j,k,l);

			compute_PRiloop(i,j,k,l);

			compute_PRmloop(i,j,k,l);

			compute_PRmloop0(i,j,k,l);

			compute_PRmloop1(i,j,k,l);

			compute_PMiloop(i,j,k,l);

			compute_PMmloop(i,j,k,l);

			compute_PMmloop0(i,j,k,l);

			compute_PMmloop1(i,j,k,l);

			compute_POiloop(i,j,k,l);

			compute_POmloop(i,j,k,l);

			compute_POmloop0(i,j,k,l);

			compute_POmloop1(i,j,k,l);

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

void pseudo_loop::compute_WBP(int i, int l){
	energy_t min_energy= INF, b1 = INF, b2=INF;

	cand_pos_t il = index[i]+l-i;
	// Mateo 2025 Changed i to i+1 as get_WB(i,i-1) is wrong
	for(cand_pos_t d=i+1; d< l; ++d){
		for(cand_pos_t e = d+1; e<= l; ++e){
			// Hosna, August 26, 2014
			// comparing calculation of WI in HFold and WPP in CCJ, I found that
			// in HFold we add another penalty called PPS_penalty for closed regions inside a pseudoloop or multiloop that spans a band
			energy_t common = get_WB(i,d-1) + beta1P*(l-e)+PPS_penalty;
			b1 = V->get_energy(d,e) + beta2P(e,d) + common;
			b2 = get_P(d,e) + gamma0m + common;
			min_energy = std::min({min_energy,b1,b2});
		}
	}
	if (min_energy < INF/2){
		WBP[il] = min_energy;
	}
}

void pseudo_loop::compute_WPP(cand_pos_t i, cand_pos_t l){
	energy_t min_energy = INF, b1 = INF, b2=INF;

	cand_pos_t il = index[i]+l-i;
	for(cand_pos_t d=i+1; d<l; ++d){
		for(cand_pos_t e = d+1; e<= l; ++e){
			// Hosna, August 26, 2014
			// comparing calculation of WI in HFold and WPP in CCJ, I found that
			// in HFold we add another penalty called PPS_penalty for closed regions inside a pseudoloop or multiloop that spans a band
			energy_t common = get_WP(i,d-1) + gamma1*(l-e) +PPS_penalty;
			b1 = V->get_energy(d,e) + gamma2(e,d) + common;
			b2 = get_P(d,e) + gamma0P + common;
			min_energy = std::min({min_energy,b1,b2});
		}
	}
	if (min_energy < INF/2){
		WPP[il] = min_energy;
	}
}

void pseudo_loop::compute_P(cand_pos_t i, cand_pos_t l){
	energy_t min_energy = INF,b1=INF;
	cand_pos_t il = index[i]+l-i;
	for(cand_pos_t j=i; j<l; ++j){
		for (cand_pos_t d=j+1; d<l; ++d){
			for (cand_pos_t k=d+1; k<l; ++k){
				b1 = get_PK(i,j-1,d+1,k-1) + get_PK(j,d,k,l);
				min_energy = std::min(min_energy,b1);
			}
		}
	}
	if (min_energy < INF/2){
		P[il]=min_energy;
	}
}

void pseudo_loop::compute_PK(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	energy_t min_energy = INF,b1=INF,b2=INF,b3=INF,b4=INF,b5=INF,b6=INF;

	cand_pos_t ij = index[i]+j-i;
	cand_pos_t kl = index[k]+l-k;

	// Hosna, july 8, 2014
	// based on original recurrences we should have i<d, and
	// it is not clear to me why we have i<=d here, so I am changing this back to original
	// by changing d=i to d=i+1
	for(cand_pos_t d=i+1; d< j; ++d){
		energy_t tmp = get_PK(i,d,k,l) + get_WP(d+1,j);
		b1=std::min(b1,tmp);
	}

	// Hosna, july 8, 2014
	// based on original recurrences we should have d<l, and
	// it is not clear to me why we have d<=l here, so I am changing this back to original
	// by changing d<=l to d<l
	for(int d=k+1; d< l; ++d){
		energy_t tmp = get_PK(i,j,d,l) + get_WP(k,d-1);
		b2=std::min(b2,tmp);
	}

	b3 = get_PL(i,j,k,l) + gamma2(j,i)+PB_penalty;
	b4 = get_PM(i,j,k,l) + gamma2(j,k)+PB_penalty;
	b5 = get_PR(i,j,k,l) + gamma2(l,k)+PB_penalty;
	b6 = get_PO(i,j,k,l) + gamma2(l,i)+PB_penalty;
	min_energy = std::min({min_energy,b1,b2,b3,b4,b5,b6});
	if (min_energy < INF/2){
		PK[ij][kl]=min_energy;
	}

}

void pseudo_loop::compute_PL(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	energy_t min_energy = INF,b1=INF,b2=INF,b3=INF;

	cand_pos_t ij = index[i]+j-i;
	cand_pos_t kl = index[k]+l-k;
	const int ptype_closing = pair[S_[i]][S_[j]];

	if (ptype_closing>0){

		b1 = get_PLiloop(i,j,k,l);

		// Hosna, April 11, 2014
		// we need to add a branch penalty for multiloop that spans a band
		b2 = get_PLmloop(i,j,k,l) + bp_penalty;

		// Hosna, July 11, 2014
		// To avoid addition of close base pairs we check for the following here
		if (j>=(i+TURN+1)){
			// Hosna April 11, 2014
			// I think we have already added gamma2(j,i) when coming to PL, so there is no need to add it again here.
			// Hosna July 17, 2014
			// I am adding gamma2 back here to avoid single base pair band predictions
			b3 = get_PfromL(i+1,j-1,k,l) + gamma2(j,i);
		}
	}
	min_energy = std::min({min_energy,b1,b2,b3});
	if (min_energy < INF/2){
		PL[ij][kl]=min_energy;
	}
}

void pseudo_loop::compute_PR(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	energy_t min_energy = INF,b1=INF,b2=INF,b3=INF;

	cand_pos_t ij = index[i]+j-i;
	cand_pos_t kl = index[k]+l-k;

	const int ptype_closing = pair[S_[k]][S_[l]];

	if (ptype_closing>0){
		b1 = get_PRiloop(i,j,k,l);

		// Hosna, April 11, 2014
		// we need to add a branch penalty for multiloop that spans a band
		b2 = get_PRmloop(i,j,k,l)+ bp_penalty;

		// Hosna, July 11, 2014
		// To avoid addition of close base pairs we check for the following here
		if (l>=(k+TURN+1)){
			// Hosna April 11, 2014
			// I think we have already added gamma2(l,k) when coming to PR, so there is no need to add it again here.
			// Hosna July 17, 2014
			// I am adding gamma2 back here to avoid single base pair band predictions
			b3 = get_PfromR(i,j,k+1,l-1) + gamma2(l,k);
		}
	}
	min_energy = std::min({min_energy,b1,b2,b3});
	if (min_energy < INF/2){
		PR[ij][kl]=min_energy;
	}
}

void pseudo_loop::compute_PM(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	energy_t min_energy = INF,b1=INF,b2=INF,b3=INF,b4=INF;

	cand_pos_t ij = index[i]+j-i;
	cand_pos_t kl = index[k]+l-k;

	const int ptype_closing = pair[S_[j]][S_[k]];
	if (ptype_closing>0){
		b1 = get_PMiloop(i,j,k,l);
		// Hosna, April 11, 2014
		// we need to add a branch penalty for multiloop that spans a band
		b2 = get_PMmloop(i,j,k,l) + bp_penalty;

		// Hosna, July 11, 2014
		// To avoid addition of close base pairs we check for the following here
		if (k>=(j+TURN-1)){
			// Hosna April 11, 2014
			// I think we have already added gamma2(j,k) when coming to PM, so there is no need to add it again here.
			// Hosna July 17, 2014
			// I am adding gamma2 back here to avoid single base pair band predictions
			b3 = get_PfromM(i,j-1,k+1,l) + gamma2(j,k);
		}

		// Hosna April 11, 2014
		// adding calculation of branch 4 here too
		if(i==j && k==l){
			b4=gamma2(i,l);
		}
	}
	min_energy = std::min({min_energy,b1,b2,b3,b4});
	if (min_energy < INF/2){
		PM[ij][kl]=min_energy;
	}
}

void pseudo_loop::compute_PO(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	energy_t min_energy = INF,b1=INF,b2=INF,b3=INF;

	cand_pos_t ij = index[i]+j-i;
	cand_pos_t kl = index[k]+l-k;

	const int ptype_closing = pair[S_[i]][S_[l]];

	if (ptype_closing>0){
		b1 = get_POiloop(i,j,k,l);

		// Hosna, April 11, 2014
		// we need to add a branch penalty for multiloop that spans a band
		b2 = get_POmloop(i,j,k,l)+ bp_penalty;

		// Hosna, July 11, 2014
		// To avoid addition of close base pairs we check for the following here
		if (l>=(i+TURN+1)){
			// Hosna April 11, 2014
			// I think we have already added gamma2(l,i) when coming to PO, so there is no need to add it again here.
			// Hosna July 17, 2014
			// I am adding gamma2 back here to avoid single base pair band predictions
			b3 = get_PfromO(i+1,j,k,l-1) + gamma2(l,i);
		}
	}
	min_energy = std::min({min_energy,b1,b2,b3});
	if (min_energy < INF/2){
		PO[ij][kl]=min_energy;
	}
}

void pseudo_loop::compute_PfromL(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	energy_t min_energy = INF,b1=INF,b2=INF,b3=INF,b4=INF,b5=INF;

	cand_pos_t ij = index[i]+j-i;
	cand_pos_t kl = index[k]+l-k;

	for(cand_pos_t d=i+1; d< j; ++d){
		energy_t tmp = get_PfromL(d,j,k,l)+get_WP(i,d-1);
		b1 = std::min(b1,tmp);
		tmp = get_PfromL(i,d,k,l)+get_WP(d+1,j);
		b2 = std::min(b2,tmp);
	}
	//Hosna, July 28, 2014
	// I think going from PfromX to PX we are changing bands and so should be paying a band penalty
	b3 = get_PR(i,j,k,l) + gamma2(l,k) + PB_penalty; //;

	//Hosna, July 28, 2014
	// I think going from PfromX to PX we are changing bands and so should be paying a band penalty
	b4 = get_PM(i,j,k,l) + gamma2(j,k)+ PB_penalty;//;

	//Hosna, July 28, 2014
	// I think going from PfromX to PX we are changing bands and so should be paying a band penalty
	b5 = get_PO(i,j,k,l) + gamma2(l,i)+ PB_penalty;//;

	min_energy = std::min({min_energy,b1,b2,b3,b4,b5});
	if (min_energy < INF/2){
		PfromL[ij][kl]=min_energy;
	}
}

void pseudo_loop::compute_PfromR(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	energy_t min_energy = INF,b1=INF,b2=INF,b3=INF,b4=INF;

	cand_pos_t ij = index[i]+j-i;
	cand_pos_t kl = index[k]+l-k;
	for(cand_pos_t d=k+1; d< l; ++d){
		energy_t tmp = get_PfromR(i,j,d,l)+get_WP(k,d-1);
		b1 = std::min(b1,tmp);
		tmp = get_PfromR(i,j,k,d)+get_WP(d+1,l);
		b2 = std::min(b2,tmp);
	}
	//Hosna, July 28, 2014
	// I think going from PfromX to PX we are changing bands and so should be paying a band penalty
	b3 = get_PM(i,j,k,l) + gamma2(j,k)+ PB_penalty;//;

	//Hosna, July 28, 2014
	// I think going from PfromX to PX we are changing bands and so should be paying a band penalty
	b4 = get_PO(i,j,k,l) + gamma2(l,i)+ PB_penalty;//;

	min_energy = std::min({min_energy,b1,b2,b3,b4});
	if (min_energy < INF/2){
		PfromR[ij][kl]=min_energy;
	}
}

void pseudo_loop::compute_PfromM(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	energy_t min_energy = INF,b1=INF,b2=INF,b3=INF,b4=INF; //b5=INF;

	cand_pos_t ij = index[i]+j-i;
	cand_pos_t kl = index[k]+l-k;

	for(cand_pos_t d=i+1; d<j; ++d){
		energy_t tmp = get_PfromM(i,d,k,l)+get_WP(d+1,j);
		b1 = std::min(b1,tmp);
	}
	for(cand_pos_t d=k+1; d<l; ++d){
		energy_t tmp=get_PfromM(i,j,d,l)+get_WP(k,d-1);
		b2 = std::min(b2,tmp);
	}
	//Hosna, July 28, 2014
	// I think going from PfromX to PX we are changing bands and so should be paying a band penalty
	b3 = get_PL(i,j,k,l) + gamma2(j,i)+ PB_penalty;//;
	//Hosna, July 28, 2014
	// I think going from PfromX to PX we are changing bands and so should be paying a band penalty
	b4 = get_PR(i,j,k,l) + gamma2(l,k)+PB_penalty;

	min_energy = std::min({min_energy,b1,b2,b3,b4});
	if (min_energy < INF/2){
		PfromM[ij][kl]=min_energy;
	}
}

void pseudo_loop::compute_PfromO(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	energy_t min_energy = INF,b1=INF,b2=INF,b3=INF,b4=INF;

	cand_pos_t ij = index[i]+j-i;
	cand_pos_t kl = index[k]+l-k;

	for(cand_pos_t d=i+1; d< j; ++d){
		energy_t tmp=get_PfromO(d,j,k,l)+get_WP(i,d-1);
		b1 = std::min(b1,tmp);
	}
	for(cand_pos_t d=k+1; d<l; ++d){
		energy_t tmp=get_PfromO(i,j,k,d)+get_WP(d+1,l);
		b2 = std::min(b2,tmp);
	}

	//Hosna, July 28, 2014
	// I think going from PfromX to PX we are changing bands and so should be paying a band penalty
	b3 = get_PL(i,j,k,l) + gamma2(j,i) + PB_penalty;

	//Hosna, July 28, 2014
	// I think going from PfromX to PX we are changing bands and so should be paying a band penalty
	b4 = get_PR(i,j,k,l) + gamma2(l,k)+PB_penalty;

	min_energy = std::min({min_energy,b1,b2,b3,b4});
	if (min_energy < INF/2){
		PfromO[ij][kl]=min_energy;
	}

}

void pseudo_loop::compute_PLiloop(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	energy_t min_energy = INF,b1=INF,b2=INF,tmp=INF;

	cand_pos_t ij = index[i]+j-i;
	cand_pos_t kl = index[k]+l-k;
	int ptype_closing = pair[S_[i]][S_[j]];
	if (ptype_closing>0){
		b1 = get_PL(i+1,j-1,k,l) + get_e_stP(i,j);

		for(cand_pos_t d= i+1; d<std::min(j,i+MAXLOOP); ++d){
        	for(cand_pos_t dp = j-1; dp > std::max(d+TURN,j-MAXLOOP); --dp){
            	tmp = get_e_intP(i,d,dp,j) + get_PL(d,dp,k,l);
            	b2 = std::min(b2,tmp);
        	}
    	}
		min_energy = std::min(b1,b2);
	}
	if (min_energy < INF/2){
		PLiloop[ij][kl]=min_energy;
	}
}

void pseudo_loop::compute_PLmloop(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	energy_t min_energy = INF,b1=INF,b2=INF;

	cand_pos_t ij = index[i]+j-i;
	cand_pos_t kl = index[k]+l-k;

	for(cand_pos_t d = i+1; d < j-1; ++d){
		// Hosna, Feb 23, 2014
		// Since PLmloop comes from PL and in PL i and j must pair, so we must have the same thing here
		// therefore we should have PLmloop0(d,j-1,k,l)
		energy_t tmp = get_PLmloop0(d,j-1,k,l) + get_WBP(i+1,d-1) + beta0P + beta2P(j,i);
		b1 = std::min(b1,tmp);
		energy_t tmp2 = get_PLmloop1(d,j-1,k,l) + get_WB(i+1,d-1) + beta0P + beta2P(j,i);
		b2 = std::min(b2,tmp2);
	}
	min_energy = std::min(b1,b2);
	if (min_energy < INF/2){
		PLmloop[ij][kl]=min_energy;
	}
}

void pseudo_loop::compute_PLmloop0(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	energy_t min_energy = INF;

	cand_pos_t ij = index[i]+j-i;
	cand_pos_t kl = index[k]+l-k;
	for(cand_pos_t d = i+1; d < j; ++d){
		// Hosna, feb 23, 2014
		// changed the recurrences so that j-1 is accounted for in PLmloop
		energy_t tmp = get_PL(i,d,k,l) + get_WB(d+1,j) + beta2P(d,i);
		min_energy = std::min(min_energy,tmp);
	}
	if (min_energy < INF/2){
		PLmloop0[ij][kl] = min_energy;
	}

}

void pseudo_loop::compute_PLmloop1(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	cand_pos_t min_energy = INF;

	cand_pos_t ij = index[i]+j-i;
	cand_pos_t kl = index[k]+l-k;
	for(cand_pos_t d = i+1; d < j; ++d){
		// Hosna, feb 23, 2014
		// changed the recurrences so that j-1 is accounted for in PLmloop
		energy_t tmp = get_PL(i,d,k,l) + get_WBP(d+1,j) + beta2P(d,i);
		min_energy = std::min(min_energy,tmp);
	}
	if (min_energy < INF/2){
		PLmloop1[ij][kl] = min_energy;
	}

}

void pseudo_loop::compute_PRiloop(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	energy_t min_energy = INF,b1=INF,b2=INF,tmp=INF;

	cand_pos_t ij = index[i]+j-i;
	cand_pos_t kl = index[k]+l-k;
	int ptype_closing = pair[S_[k]][S_[l]];
	if (ptype_closing>0){
		b1 = get_PR(i,j,k+1,l-1) + get_e_stP(k,l);

		for(cand_pos_t d= k+1; d<std::min(l,k+MAXLOOP); ++d){
        	for(cand_pos_t dp=l-1; dp > std::max(d+TURN,l-MAXLOOP); --dp){
            	tmp = get_e_intP(k,d,dp,l) + get_PR(i,j,d,dp);
            	b2 = std::min(b2,tmp);
        	}
    	}
		min_energy = std::min(b1,b2);
	}
	if (min_energy < INF/2){
		PRiloop[ij][kl]=min_energy;
	}

}

void pseudo_loop::compute_PRmloop(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){

	energy_t min_energy = INF,b1=INF,b2=INF;
	cand_pos_t ij = index[i]+j-i;
	cand_pos_t kl = index[k]+l-k;

	//Hosna Feb 23, 2014
	// changed the recurrences to have l-1 instead of l in PRmloop and removed l-1 from PRmloop0,1, as we are accounting for k.l here
	for(cand_pos_t d = k+1; d < l-1; ++d){
		energy_t tmp = get_PRmloop0(i,j,d,l-1) + get_WBP(k+1,d-1) + beta0P + beta2P(l,k);
		b1 = std::min(b1,tmp);
		tmp = get_PRmloop1(i,j,d,l-1) + get_WB(k+1,d-1) + beta0P + beta2P(l,k);
		b2 = std::min(b2,tmp);
	}
	min_energy = std::min(b1,b2);
	if (min_energy < INF/2){
		PRmloop[ij][kl]=min_energy;
	}

}


void pseudo_loop::compute_PRmloop0(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	energy_t min_energy = INF;

	cand_pos_t ij = index[i]+j-i;
	cand_pos_t kl = index[k]+l-k;
	for(cand_pos_t d = k+1; d < l; ++d){
		energy_t tmp = get_PR(i,j,k,d) + get_WB(d+1,l) + beta2P(k,d);
		min_energy = std::min(min_energy,tmp);
	}
	if (min_energy < INF/2){
		PRmloop0[ij][kl] = min_energy;
	}

}

void pseudo_loop::compute_PRmloop1(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	energy_t min_energy = INF;

	cand_pos_t ij = index[i]+j-i;
	cand_pos_t kl = index[k]+l-k;
	for(cand_pos_t d = k+1; d < l; ++d){
		energy_t tmp = get_PR(i,j,k,d) + get_WBP(d+1,l) + beta2P(k,d);
		min_energy = std::min(min_energy,tmp);
	}
	if (min_energy < INF/2){
		PRmloop1[ij][kl] = min_energy;
	}

}

void pseudo_loop::compute_PMiloop(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	energy_t min_energy = INF,b1=INF,b2=INF,tmp=INF;

	cand_pos_t ij = index[i]+j-i;
	cand_pos_t kl = index[k]+l-k;
	int ptype_closing = pair[S_[j]][S_[k]];

	if (ptype_closing>0){
		b1 = get_PM(i,j-1,k+1,l) + get_e_stP(j-1,k+1); // Why is this j-1,k+1, it makes no sense. Everything else is like j,k

		for(cand_pos_t d= j-1; d>std::max(i,j-MAXLOOP); --d){
        	for (cand_pos_t dp=k+1; dp <std::min(l,k+MAXLOOP); ++dp) {
            	tmp = get_e_intP(d,j,k,dp) + get_PM(i,d,dp,l);
            	b2= std::min(b2,tmp);
        	}
    	}
		min_energy = std::min(b1,b2);
	}
	if (min_energy < INF/2){
		PMiloop[ij][kl]=min_energy;
	}

}

void pseudo_loop::compute_PMmloop(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	energy_t min_energy = INF,b1=INF,b2=INF;

	cand_pos_t ij = index[i]+j-i;
	cand_pos_t kl = index[k]+l-k;

	// Hosna Feb 23, 2014
	// changed the recurrence of PMmloop, to have k+1 instead of k when calling PMmloop0,1, and removed k+1 from PMmloop0,1
	// as j.k pair is accounted for in this recurrence
	for(cand_pos_t d = i+1; d < j; ++d){
		energy_t tmp = get_PMmloop0(i,d,k+1,l) + get_WBP(d+1,j-1) + beta0P + beta2P(j,k);
		b1 = std::min(b1,tmp);
		energy_t tmp2 = get_PMmloop1(i,d,k+1,l) + get_WB(d+1,j-1) + beta0P + beta2P(j,k);
		b2 = std::min(b2,tmp2);
	}
	min_energy = std::min(b1,b2);

	if (min_energy < INF/2){
		PMmloop[ij][kl]=min_energy;
	}
}


void pseudo_loop::compute_PMmloop0(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	energy_t min_energy = INF;

	cand_pos_t ij = index[i]+j-i;
	cand_pos_t kl = index[k]+l-k;
	for(cand_pos_t d = k+1; d < l; ++d){
		energy_t tmp = get_PM(i,j,d,l) + get_WB(k,d-1) + beta2P(j,d);
		min_energy = std::min(min_energy,tmp);
	}

	if (min_energy < INF/2){
		PMmloop0[ij][kl] = min_energy;
	}
}

void pseudo_loop::compute_PMmloop1(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	energy_t min_energy = INF;

	cand_pos_t ij = index[i]+j-i;
	cand_pos_t kl = index[k]+l-k;
	for(cand_pos_t d = k+1; d < l; ++d){
		energy_t tmp = get_PM(i,j,d,l) + get_WBP(k,d-1) + beta2P(j,d);
		min_energy = std::min(min_energy,tmp);
	}

	if (min_energy < INF/2){
		PMmloop1[ij][kl] = min_energy;
	}
}

void pseudo_loop::compute_POiloop(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	energy_t min_energy = INF,b1=INF,b2=INF,tmp=INF;

	cand_pos_t ij = index[i]+j-i;
	cand_pos_t kl = index[k]+l-k;
	int ptype_closing = pair[S_[i]][S_[l]];

	if (ptype_closing>0){
		b1 = get_PO(i+1,j,k,l-1) + get_e_stP(i,l);

		for(cand_pos_t d= i+1; d<std::min(j,i+MAXLOOP); ++d){
        	for (cand_pos_t dp=l-1; dp >std::max(l-MAXLOOP,k); --dp) {
            	tmp = get_e_intP(i,d,dp,l) + get_PO(d,j,dp,k);
            	b2 = std::min(b2,tmp);
        	}
    	}
		min_energy = std::min(b1,b2);
	}
	if (min_energy < INF/2){
		POiloop[ij][kl]=min_energy;
	}

}

void pseudo_loop::compute_POmloop(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	energy_t min_energy = INF,b1=INF,b2=INF;

	cand_pos_t ij = index[i]+j-i;
	cand_pos_t kl = index[k]+l-k;

	// Hosna Feb 23, 2014
	// changed recurrences for POmloop to have l-1 instead of l and removed l-1 from POmloop0,1 as i.l is accounted for here in this recurrence
	for(cand_pos_t d = i+1; d < j; ++d){
		energy_t tmp = get_POmloop0(d,j,k,l-1) + get_WBP(i+1,d-1) + beta0P + beta2P(l,i);
		b1 = std::min(b1,tmp);
		energy_t tmp2 = get_POmloop1(d,j,k,l-1) + get_WB(i+1,d-1) + beta0P + beta2P(l,i);
		b2 = std::min(b2,tmp2);
	}
	min_energy = std::min(b1,b2);

	if (min_energy < INF/2){
		POmloop[ij][kl]=min_energy;
	}

}


void pseudo_loop::compute_POmloop0(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	energy_t min_energy = INF;

	cand_pos_t ij = index[i]+j-i;
	cand_pos_t kl = index[k]+l-k;
	for(cand_pos_t d = k+1; d < l; ++d){
		energy_t tmp = get_PO(i,j,k,d) + get_WB(d+1,l) + beta2P(d,i);
		min_energy = std::min(min_energy,tmp);
	}

	if (min_energy < INF/2){
		POmloop0[ij][kl] = min_energy;
	}

}

void pseudo_loop::compute_POmloop1(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	energy_t min_energy = INF;

	cand_pos_t ij = index[i]+j-i;
	cand_pos_t kl = index[k]+l-k;
	for(cand_pos_t d = k+1; d < l; ++d){
		energy_t tmp = get_PO(i,j,k,d) + get_WBP(d+1,l) + beta2P(d,i);
		min_energy = std::min(min_energy,tmp);
	}

	if (min_energy < INF/2){
		POmloop1[ij][kl] = min_energy;
	}
}


energy_t pseudo_loop::get_WB(cand_pos_t i, cand_pos_t j){
	if (i<=0 || j<=0 || i>n || j>n){
		return INF;
	}
	return (std::min(beta1P*(j-i+1),get_WBP(i,j)));
}

energy_t pseudo_loop::get_WBP(cand_pos_t i, cand_pos_t j){
	if (i > j || i<=0 || j<=0 || i>n || j>n){
		return INF;
	}
	cand_pos_t ij = index[i]+j-i;
	return WBP[ij];
}

energy_t pseudo_loop::get_WP(cand_pos_t i, cand_pos_t j){
	if (i<=0 || j<=0 || i>n || j>n){
		return INF;
	}
	if (i>j) return 0; // needed as this will happen 
	return (std::min(gamma1*(j-i+1),get_WPP(i,j)));
}

energy_t pseudo_loop::get_WPP(cand_pos_t i, cand_pos_t j){
	if (i > j  || i<=0 || j<=0 || i>n || j>n){
		return INF;
	}
	cand_pos_t ij = index[i]+j-i;
	return WPP[ij];
}

energy_t pseudo_loop::get_P(cand_pos_t i, cand_pos_t j){
	if (i >= j  || i<=0 || j<=0 || i>n || j>n){
		return INF;
	}
	cand_pos_t ij = index[i]+j-i;
	return P[ij];
}

energy_t pseudo_loop::get_PK(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	if (!(i <= j && j < k-1 && k <= l)){
		return INF;
	}
	// Hosna, April 3, 2014
	// adding impossible cases
	// Mateo -These are needed as k > l due to some of the function doing k+1 when k==l. I would prefer to make it so that won't happen so that these ifs 
	// are not needed as they are most likely expensive time-wise
	if (i<=0 ||j<=0 || k<=0 || l<=0 || i>n || j>n || k> n || l> n){
		return INF;
	}
	cand_pos_t ij = index[i]+j-i;
	cand_pos_t kl = index[k]+l-k;

	return PK[ij][kl];
}

energy_t pseudo_loop::get_PL(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	if (!(i <= j && j < k-1 && k <= l)){
		return INF;
	}
	// Hosna, April 3, 2014
	// adding impossible cases
	// Mateo -These are needed as k > l due to some of the function doing k+1 when k==l. I would prefer to make it so that won't happen so that these ifs 
	// are not needed as they are most likely expensive time-wise
	if (i<=0 ||j<=0 || k<=0 || l<=0 || i>n || j>n || k> n || l> n){
		return INF;
	}
	cand_pos_t ij = index[i]+j-i;
	cand_pos_t kl = index[k]+l-k;

	return PL[ij][kl];
}

energy_t pseudo_loop::get_PR(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	if (!(i <= j && j < k-1 && k <= l)){
		return INF;
	}
	// Hosna, April 3, 2014
	// adding impossible cases
	// Mateo -These are needed as k > l due to some of the function doing k+1 when k==l. I would prefer to make it so that won't happen so that these ifs 
	// are not needed as they are most likely expensive time-wise
	if (i<=0 ||j<=0 || k<=0 || l<=0 || i>n || j>n || k> n || l> n){
		return INF;
	}
	cand_pos_t ij = index[i]+j-i;
	cand_pos_t kl = index[k]+l-k;

	return PR[ij][kl];
}

energy_t pseudo_loop::get_PM(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	if (!(i <= j && j < k-1 && k <= l)){
		return INF;
	}
	// Hosna, April 3, 2014
	// adding impossible cases
	// Mateo -These are needed as k > l due to some of the function doing k+1 when k==l. I would prefer to make it so that won't happen so that these ifs 
	// are not needed as they are most likely expensive time-wise
	if (i<=0 ||j<=0 || k<=0 || l<=0 || i>n || j>n || k> n || l> n){
		return INF;
	}
	int ptype_closing = pair[S_[j]][S_[k]];
	if (ptype_closing == 0){
		return INF;
	}
	cand_pos_t ij = index[i]+j-i;
	cand_pos_t kl = index[k]+l-k;
	if (i ==j && k ==l){
		return gamma2(i,l);
	}

	return PM[ij][kl];
}

energy_t pseudo_loop::get_PO(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	if (!(i <= j && j < k-1 && k <= l)){
		return INF;
	}
	// Hosna, April 3, 2014
	// adding impossible cases
	// Mateo -These are needed as k > l due to some of the function doing k+1 when k==l. I would prefer to make it so that won't happen so that these ifs 
	// are not needed as they are most likely expensive time-wise
	if (i<=0 ||j<=0 || k<=0 || l<=0 || i>n || j>n || k> n || l> n){
		return INF;
	}
	cand_pos_t ij = index[i]+j-i;
	cand_pos_t kl = index[k]+l-k;

	return PO[ij][kl];
}

energy_t pseudo_loop::get_PfromL(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	if (!(i <= j && j < k-1 && k <= l)){
		return INF;
	}
	// Hosna, April 3, 2014
	// adding impossible cases
	// Mateo -These are needed as k > l due to some of the function doing k+1 when k==l. I would prefer to make it so that won't happen so that these ifs 
	// are not needed as they are most likely expensive time-wise
	if (i<=0 ||j<=0 || k<=0 || l<=0 || i>n || j>n || k> n || l> n){
		return INF;
	}
	if (i==j && k==l){
		// Hosna, August 13, 2014
		// in some cases it used P_fromM to exit from a case with no band!
		const int ptype_closing = pair[S_[i]][S_[l]];
		if(ptype_closing == 0) return INF;
		else return gamma2(j,k) + gamma2(k,j);
	}
	cand_pos_t ij = index[i]+j-i;
	cand_pos_t kl = index[k]+l-k;

	return PfromL[ij][kl];
}

energy_t pseudo_loop::get_PfromR(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	if (!(i <= j && j < k-1 && k <= l)){
		return INF;
	}
	// Hosna, April 3, 2014
	// adding impossible cases
	// Mateo -These are needed as k > l due to some of the function doing k+1 when k==l. I would prefer to make it so that won't happen so that these ifs 
	// are not needed as they are most likely expensive time-wise
	if (i<=0 ||j<=0 || k<=0 || l<=0 || i>n || j>n || k> n || l> n){
		return INF;
	}
	if (i==j && k==l){
		// Hosna, August 13, 2014
		// in some cases it used P_fromM to exit from a case with no band!
		const int ptype_closing = pair[S_[i]][S_[l]];
		if(ptype_closing == 0) return INF;
		else return gamma2(j,k) + gamma2(k,j);
	}
	cand_pos_t ij = index[i]+j-i;
	cand_pos_t kl = index[k]+l-k;

	return PfromR[ij][kl];
}

energy_t pseudo_loop::get_PfromM(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	if (!(i <= j && j < k-1 && k <= l)){
		return INF;
	}
	// Hosna, April 3, 2014
	// adding impossible cases
	// Mateo -These are needed as k > l due to some of the function doing k+1 when k==l. I would prefer to make it so that won't happen so that these ifs 
	// are not needed as they are most likely expensive time-wise
	if (i<=0 ||j<=0 || k<=0 || l<=0 || i>n || j>n || k> n || l> n){
		return INF;
	}
	if (i==j && k==l){
		// Hosna, August 13, 2014
		// in some cases it used P_fromM to exit from a case with no band!
		const int ptype_closing = pair[S_[i]][S_[l]];
		if(ptype_closing == 0) return INF;
		else return 0;
	}
	cand_pos_t ij = index[i]+j-i;
	cand_pos_t kl = index[k]+l-k;

	return PfromM[ij][kl];
}

energy_t pseudo_loop::get_PfromO(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	if (!(i <= j && j < k-1 && k <= l)){
		return INF;
	}
	// Hosna, April 3, 2014
	// adding impossible cases
	// Mateo -These are needed as k > l due to some of the function doing k+1 when k==l. I would prefer to make it so that won't happen so that these ifs 
	// are not needed as they are most likely expensive time-wise
	if (i<=0 ||j<=0 || k<=0 || l<=0 || i>n || j>n || k> n || l> n){
		return INF;
	}
	if (i==j && k==l){
		// Hosna, August 13, 2014
		// in some cases it used P_fromM to exit from a case with no band!
		const int ptype_closing = pair[S_[i]][S_[l]];
		if(ptype_closing == 0) return INF;
		else return 0;
	}
	cand_pos_t ij = index[i]+j-i;
	cand_pos_t kl = index[k]+l-k;

	return PfromO[ij][kl];
}

energy_t pseudo_loop::get_PLiloop(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	if (!(i <= j && j < k-1 && k <= l)){
		return INF;
	}
	// Hosna, April 3, 2014
	// adding impossible cases
	// Mateo -These are needed as k > l due to some of the function doing k+1 when k==l. I would prefer to make it so that won't happen so that these ifs 
	// are not needed as they are most likely expensive time-wise
	if (i<=0 ||j<=0 || k<=0 || l<=0 || i>n || j>n || k> n || l> n){
		return INF;
	}
	const int ptype_closing = pair[S_[i]][S_[j]];
	if(ptype_closing == 0) return INF;

	cand_pos_t ij = index[i]+j-i;
	cand_pos_t kl = index[k]+l-k;

	return PLiloop[ij][kl];
}

energy_t pseudo_loop::get_PLmloop(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	if (!(i <= j && j < k-1 && k <= l)){
		return INF;
	}
	// Hosna, April 3, 2014
	// adding impossible cases
	// Mateo -These are needed as k > l due to some of the function doing k+1 when k==l. I would prefer to make it so that won't happen so that these ifs 
	// are not needed as they are most likely expensive time-wise
	if (i<=0 ||j<=0 || k<=0 || l<=0 || i>n || j>n || k> n || l> n){
		return INF;
	}
	cand_pos_t ij = index[i]+j-i;
	cand_pos_t kl = index[k]+l-k;

	return PLmloop[ij][kl];
}

energy_t pseudo_loop::get_PLmloop0(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	if (!(i <= j && j < k-1 && k <= l)){
		return INF;
	}
	// Hosna, April 3, 2014
	// adding impossible cases
	// Mateo -These are needed as k > l due to some of the function doing k+1 when k==l. I would prefer to make it so that won't happen so that these ifs 
	// are not needed as they are most likely expensive time-wise
	if (i<=0 ||j<=0 || k<=0 || l<=0 || i>n || j>n || k> n || l> n){
		return INF;
	}
	cand_pos_t ij = index[i]+j-i;
	cand_pos_t kl = index[k]+l-k;

	return PLmloop0[ij][kl];
}

energy_t pseudo_loop::get_PLmloop1(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	if (!(i <= j && j < k-1 && k <= l)){
		return INF;
	}
	// Hosna, April 3, 2014
	// adding impossible cases
	// Mateo -These are needed as k > l due to some of the function doing k+1 when k==l. I would prefer to make it so that won't happen so that these ifs 
	// are not needed as they are most likely expensive time-wise
	if (i<=0 ||j<=0 || k<=0 || l<=0 || i>n || j>n || k> n || l> n){
		return INF;
	}
	cand_pos_t ij = index[i]+j-i;
	cand_pos_t kl = index[k]+l-k;

	return PLmloop1[ij][kl];
}

energy_t pseudo_loop::get_PRiloop(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	if (!(i <= j && j < k-1 && k <= l)){
		return INF;
	}
	// Hosna, April 3, 2014
	// adding impossible cases
	// Mateo -These are needed as k > l due to some of the function doing k+1 when k==l. I would prefer to make it so that won't happen so that these ifs 
	// are not needed as they are most likely expensive time-wise
	if (i<=0 ||j<=0 || k<=0 || l<=0 || i>n || j>n || k> n || l> n){
		return INF;
	}
	const int ptype_closing = pair[S_[k]][S_[l]];
	if(ptype_closing == 0) return INF;

	cand_pos_t ij = index[i]+j-i;
	cand_pos_t kl = index[k]+l-k;

	return PRiloop[ij][kl];
}

energy_t pseudo_loop::get_PRmloop(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	if (!(i <= j && j < k-1 && k <= l)){
		return INF;
	}
	// Hosna, April 3, 2014
	// adding impossible cases
	// Mateo -These are needed as k > l due to some of the function doing k+1 when k==l. I would prefer to make it so that won't happen so that these ifs 
	// are not needed as they are most likely expensive time-wise
	if (i<=0 ||j<=0 || k<=0 || l<=0 || i>n || j>n || k> n || l> n){
		return INF;
	}
	cand_pos_t ij = index[i]+j-i;
	cand_pos_t kl = index[k]+l-k;

	return PRmloop[ij][kl];
}

energy_t pseudo_loop::get_PRmloop0(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	if (!(i <= j && j < k-1 && k <= l)){
		return INF;
	}
	// Hosna, April 3, 2014
	// adding impossible cases
	// Mateo -These are needed as k > l due to some of the function doing k+1 when k==l. I would prefer to make it so that won't happen so that these ifs 
	// are not needed as they are most likely expensive time-wise
	if (i<=0 ||j<=0 || k<=0 || l<=0 || i>n || j>n || k> n || l> n){
		return INF;
	}
	cand_pos_t ij = index[i]+j-i;
	cand_pos_t kl = index[k]+l-k;

	return PRmloop0[ij][kl];
}

energy_t pseudo_loop::get_PRmloop1(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	if (!(i <= j && j < k-1 && k <= l)){
		return INF;
	}
	// Hosna, April 3, 2014
	// adding impossible cases
	// Mateo -These are needed as k > l due to some of the function doing k+1 when k==l. I would prefer to make it so that won't happen so that these ifs 
	// are not needed as they are most likely expensive time-wise
	if (i<=0 ||j<=0 || k<=0 || l<=0 || i>n || j>n || k> n || l> n){
		return INF;
	}
	cand_pos_t ij = index[i]+j-i;
	cand_pos_t kl = index[k]+l-k;

	return PRmloop1[ij][kl];
}

energy_t pseudo_loop::get_PMiloop(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	if (!(i <= j && j < k-1 && k <= l)){
		return INF;
	}
	// Hosna, April 3, 2014
	// adding impossible cases
	// Mateo -These are needed as k > l due to some of the function doing k+1 when k==l. I would prefer to make it so that won't happen so that these ifs 
	// are not needed as they are most likely expensive time-wise
	if (i<=0 ||j<=0 || k<=0 || l<=0 || i>n || j>n || k> n || l> n){
		return INF;
	}
	const int ptype_closing = pair[S_[j]][S_[k]];
	if(ptype_closing == 0) return INF;

	cand_pos_t ij = index[i]+j-i;
	cand_pos_t kl = index[k]+l-k;

	return PMiloop[ij][kl];
}

energy_t pseudo_loop::get_PMmloop(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	if (!(i <= j && j < k-1 && k <= l)){
		return INF;
	}
	// Hosna, April 3, 2014
	// adding impossible cases
	// Mateo -These are needed as k > l due to some of the function doing k+1 when k==l. I would prefer to make it so that won't happen so that these ifs 
	// are not needed as they are most likely expensive time-wise
	if (i<=0 ||j<=0 || k<=0 || l<=0 || i>n || j>n || k> n || l> n){
		return INF;
	}
	cand_pos_t ij = index[i]+j-i;
	cand_pos_t kl = index[k]+l-k;

	return PMmloop[ij][kl];
}

energy_t pseudo_loop::get_PMmloop0(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	if (!(i <= j && j < k-1 && k <= l)){
		return INF;
	}
	// Hosna, April 3, 2014
	// adding impossible cases
	// Mateo -These are needed as k > l due to some of the function doing k+1 when k==l. I would prefer to make it so that won't happen so that these ifs 
	// are not needed as they are most likely expensive time-wise
	if (i<=0 ||j<=0 || k<=0 || l<=0 || i>n || j>n || k> n || l> n){
		return INF;
	}

	cand_pos_t ij = index[i]+j-i;
	cand_pos_t kl = index[k]+l-k;

	return PMmloop0[ij][kl];
}

energy_t pseudo_loop::get_PMmloop1(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	if (!(i <= j && j < k-1 && k <= l)){
		return INF;
	}
	// Hosna, April 3, 2014
	// adding impossible cases
	// Mateo -These are needed as k > l due to some of the function doing k+1 when k==l. I would prefer to make it so that won't happen so that these ifs 
	// are not needed as they are most likely expensive time-wise
	if (i<=0 ||j<=0 || k<=0 || l<=0 || i>n || j>n || k> n || l> n){
		return INF;
	}
	cand_pos_t ij = index[i]+j-i;
	cand_pos_t kl = index[k]+l-k;

	return PMmloop1[ij][kl];
}

energy_t pseudo_loop::get_POiloop(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	if (!(i <= j && j < k-1 && k <= l)){
		return INF;
	}
	// Hosna, April 3, 2014
	// adding impossible cases
	// Mateo -These are needed as k > l due to some of the function doing k+1 when k==l. I would prefer to make it so that won't happen so that these ifs 
	// are not needed as they are most likely expensive time-wise
	if (i<=0 ||j<=0 || k<=0 || l<=0 || i>n || j>n || k> n || l> n){
		return INF;
	}
	const int ptype_closing = pair[S_[i]][S_[l]];
	if(ptype_closing == 0) return INF;

	cand_pos_t ij = index[i]+j-i;
	cand_pos_t kl = index[k]+l-k;

	return POiloop[ij][kl];
}

energy_t pseudo_loop::get_POmloop(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	if (!(i <= j && j < k-1 && k <= l)){
		return INF;
	}
	// Hosna, April 3, 2014
	// adding impossible cases
	// Mateo -These are needed as k > l due to some of the function doing k+1 when k==l. I would prefer to make it so that won't happen so that these ifs 
	// are not needed as they are most likely expensive time-wise
	if (i<=0 ||j<=0 || k<=0 || l<=0 || i>n || j>n || k> n || l> n){
		return INF;
	}
	cand_pos_t ij = index[i]+j-i;
	cand_pos_t kl = index[k]+l-k;

	return POmloop[ij][kl];
}

energy_t pseudo_loop::get_POmloop0(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	if (!(i <= j && j < k-1 && k <= l)){
		return INF;
	}
	// Hosna, April 3, 2014
	// adding impossible cases
	// Mateo -These are needed as k > l due to some of the function doing k+1 when k==l. I would prefer to make it so that won't happen so that these ifs 
	// are not needed as they are most likely expensive time-wise
	if (i<=0 ||j<=0 || k<=0 || l<=0 || i>n || j>n || k> n || l> n){
		return INF;
	}
	cand_pos_t ij = index[i]+j-i;
	cand_pos_t kl = index[k]+l-k;

	return POmloop0[ij][kl];
}

energy_t pseudo_loop::get_POmloop1(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	if (!(i <= j && j < k-1 && k <= l)){
		return INF;
	}
	// Hosna, April 3, 2014
	// adding impossible cases
	// Mateo -These are needed as k > l due to some of the function doing k+1 when k==l. I would prefer to make it so that won't happen so that these ifs 
	// are not needed as they are most likely expensive time-wise
	if (i<=0 ||j<=0 || k<=0 || l<=0 || i>n || j>n || k> n || l> n){
		return INF;
	}
	cand_pos_t ij = index[i]+j-i;
	cand_pos_t kl = index[k]+l-k;

	return POmloop1[ij][kl];
}

energy_t pseudo_loop::compute_int(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){

	const pair_type ptype_closing = pair[S_[i]][S_[j]];
    return E_IntLoop(k-i-1,j-l-1,ptype_closing,rtype[pair[S_[k]][S_[l]]],S1_[i+1],S1_[j-1],S1_[k-1],S1_[l+1],const_cast<paramT *>(params_));
}

energy_t pseudo_loop::get_e_stP(cand_pos_t i, cand_pos_t j){
	if (i+1 == j-1){ // TODO: do I need something like that or stack is taking care of this?
		return INF;
	}
	energy_t ss = compute_int(i,j,i+1,j-1);
	return lrint(e_stP_penalty * ss);
}

energy_t pseudo_loop::get_e_intP(cand_pos_t i, cand_pos_t ip, cand_pos_t jp, cand_pos_t j){
	energy_t e_int = compute_int(i,j,ip,jp);
	energy_t energy = lrint(e_intP_penalty * e_int);
	return energy;
}

//penalty for z unpaired bases in an internal loop that spans a band
energy_t pseudo_loop::alpha1P(cand_pos_t z){
	// Hosna, April 2nd, 2014
	// for simplicity I am not checking to see what type of internal loop (i.e. I, B, or H) I have that spans a band, and I am assuming it is an intetnal loop of type I
	// Hosna, August 21, 2014
	// penalty by size returns inf for size 1
	if (z== 0) return 0;
	if (z == 1 || z==2){
		return (cand_pos_t) lrint(params_->bulge[z]*e_intP_penalty);
	}
	return (cand_pos_t) lrint(params_->internal_loop[z]*e_intP_penalty);
	
}

// penalty for closing pair i.l of an internal loop that spans a band
energy_t pseudo_loop::alpha2P(cand_pos_t i, cand_pos_t l){
	// Hosna April 2, 2014
	// I added this similar to calculation of closing base pair contribution in a general internal loop in simfold and
	// multiplied that to the intP_penalty for an internal loop that spans a band
	int ptype_closing = pair[S_[i]][S_[l]];
	return (cand_pos_t) lrint(params_->mismatchI[ptype_closing][S_[i+1]][S_[l-1]]*e_intP_penalty);
}

//penalty for asymmetry of z in an internal loop that spans a band
energy_t pseudo_loop::alpha3P(cand_pos_t z){
	if (z== 0) return 0;

	return (energy_t) lrint(std::min(MAX_NINIO,params_->ninio[2]*z)*e_intP_penalty);
}

// penalty for closing pair i.l or l.i of an ordinary multiloop
energy_t pseudo_loop::beta2(cand_pos_t i, cand_pos_t l){
	// Hosna, April 2, 2014
	// I don't think this is the complete value, but since HFold's WM and CCJ's Vmloop recurrences are the same I am not changing this value here, unless I find out it is needed
	// the correct value should be: Non-GC-penalty(i,l)+b_penalty
	return b_penalty;
}

// penalty for closing pair i.l or l.i of a multiloop that spans a band
energy_t pseudo_loop::beta2P(cand_pos_t i, cand_pos_t l){
	// Hosna, April 2, 2014
	// I don't think this is the complete value, but since HFold's WM and CCJ's Vmloop recurrences are the same I am not changing this value here, unless I find out it is needed
	// the correct value should be: Non-GC-penalty(i,l) *0.74 + bp_penalty
	return bp_penalty;
}

// penalty for closing pair i.l or l.i of a pseudoloop
energy_t pseudo_loop::gamma2(cand_pos_t i, cand_pos_t l){
	// Hosna, April 2, 2014
	// I changed this value to be 0 as I can't find its correct value
	// the correct value for this penalty should be similar to what we have in case of an internal loop or a multiloop, but the value is missing here
	return 0;
	// Hosna July 17, 2014
	// To avoid addition of single base pair bands I am giving a very small non-zero value to gamma2
	// return 1;
}

// Hosna, Feb 18, 2014
// I am changing the backtrack function such that it does not deal with structure
// instead it only fills the minimum_fold array, f, and passes it to W_final
// then in W_final one pass over f, will create the structure in dot bracket format
// This is the solution I found for the problem of not knowing what kind of brackets and
// how many different brackets to use when fillinf f and structure at the same time in pseudoloop.cpp
//void pseudo_loop::back_track(char *structure, minimum_fold *f, seq_interval *cur_interval)
void pseudo_loop::backtrack(minimum_fold *f, seq_interval *cur_interval){

	// this->structure = structure;
	this->f = f;
	switch (cur_interval->type)
	{
		case P_P:
		{
			cand_pos_t i = cur_interval->i;
			cand_pos_t l = cur_interval->j;
			if(debug) printf("At P_P with i=%d,l=%d\n",i,l);

			if (i >= l){
				std::cerr << "border case: This should not have happened!, P_P" << std::endl;
				exit(EXIT_FAILURE);
			}
			energy_t min_energy = INF,b1=INF;

			cand_pos_t best_d=0, best_j=0,best_k=0;
			for(cand_pos_t j=i; j< l; j++){
				for (cand_pos_t d=j+1; d<l; d++){
					for (cand_pos_t k=d+1; k<l; k++){
						b1 = get_PK(i,j,d+1,k) +get_PK(j+1,d,k+1,l);
						if(b1 < min_energy){
							min_energy = b1;
							best_d = d;
							best_j = j;
							best_k= k;
						}
					}
				}
			}

			insert_node(i,best_k,best_j,best_d+1,P_PK);
			insert_node(best_j+1,l,best_d,best_k+1,P_PK);
		}
		break;

		case P_PK:
		{
			cand_pos_t i = cur_interval->i;
			cand_pos_t l = cur_interval->j;
			cand_pos_t j = cur_interval->k;
			cand_pos_t k = cur_interval->l;
			if(debug) printf("At P_PK with i=%d,l=%d,j=%d,k=%d\n",i,l,j,k);
		
			if (!(i <= j && j < k-1 && k <= l)){
				std::cerr << "border cases: This should not have happened!, P_PK" << std::endl;
				exit(EXIT_FAILURE);
			}

			if (i<=0 ||j<=0 || k<=0 || l<=0 || i>n || j>n || k> n || l> n){
				std::cerr << "impossible cases: This should not have happened!, P_PK" << std::endl;
				exit(EXIT_FAILURE);
			}
			energy_t min_energy = INF,temp=INF;
			cand_pos_t best_row = -1,best_d=-1;

			// branch 1
			// Hosna, july 8, 2014
			// based on original recurrences we should have i<d, and
			// it is not clear to me why we have d=i here, so I am changing this back to original
			// by changing d=i to d=i+1
			for(cand_pos_t  d=i+1; d< j; d++){
				temp = get_PK(i,d,k,l) + get_WP(d+1,j);
				if (temp < min_energy){
					min_energy=temp;
					best_row = 1;
					best_d=d;
				}
			}

			// branch 2
			// Hosna, july 8, 2014
			// based on original recurrences we should have d<l, and
			// it is not clear to me why we have d<=l here, so I am changing this back to original
			// by changing d<=l to d<l
			for(cand_pos_t d=k+1; d< l; d++){
				temp = get_PK(i,j,d,l) + get_WP(k,d-1);
				if (temp < min_energy){
					min_energy=temp;
					best_row = 2;
					best_d = d;
				}
			}

			// branch 3
			temp = get_PL(i,j,k,l) + gamma2(j,i)+PB_penalty;
			if(temp < min_energy){
				min_energy = temp;
				best_row = 3;
				best_d = -1;
			}

			//branch 4
			temp = get_PM(i,j,k,l) + gamma2(j,k)+PB_penalty;
			if(temp < min_energy){
				min_energy = temp;
				best_row = 4;
				best_d = -1;
			}

			// branch 5
			temp = get_PR(i,j,k,l) + gamma2(l,k)+PB_penalty;
			if(temp < min_energy){
				min_energy = temp;
				best_row = 5;
				best_d = -1;
			}

			// branch 6
			temp = get_PO(i,j,k,l) + gamma2(l,i)+PB_penalty;
			if(temp < min_energy){
				min_energy = temp;
				best_row = 6;
				best_d = -1;
			}

			switch (best_row)
			{
				case 1:
					if (best_d > -1){
						insert_node(i,l,best_d,k,P_PK);
						insert_node(best_d+1,j,P_WP);
					}
					break;
				case 2:
					if (best_d > -1){
						insert_node(i,l,j,best_d,P_PK);
						insert_node(k,best_d-1,P_WP);
					}
					break;
				case 3:
					insert_node(i,l,j,k,P_PL);
					break;
				case 4:
					insert_node(i,l,j,k,P_PM);
					break;
				case 5:
					insert_node(i,l,j,k,P_PR);
					break;
				case 6:
					insert_node(i,l,j,k,P_PO);
					break;
			}
		}
		break;

		case P_PL:
		{
			cand_pos_t i = cur_interval->i;
			cand_pos_t l = cur_interval->j;
			cand_pos_t j = cur_interval->k;
			cand_pos_t k = cur_interval->l;
			if(debug) printf("At P_PL with i=%d,l=%d,j=%d,k=%d\n",i,l,j,k);

			if (!(i <= j && j < k-1 && k <= l)){
				std::cerr << "border cases: This should not have happened!, P_PL" << std::endl;
				exit(EXIT_FAILURE);
			}

			if (i<=0 ||j<=0 || k<=0 || l<=0 || i>n || j>n || k> n || l> n){
				std::cerr << "impossible cases: This should not have happened!, P_PL" << std::endl;
				exit(EXIT_FAILURE);
			}

			energy_t min_energy = INF,temp=INF;
			cand_pos_t best_row = -1;

			if (pair[S_[i]][S_[j]]> 0){

				//branch 1
				temp = get_PLiloop(i,j,k,l);
				if(temp < min_energy){
					min_energy = temp;
					best_row = 1;
				}

				//branch 2
				// Hosna, April 11, 2014
				// we need to add a branch penalty for multiloop that spans a band
				temp = get_PLmloop(i,j,k,l) + bp_penalty;
				if(temp < min_energy){
					min_energy = temp;
					best_row=2;
				}

				//branch 3
				// Hosna, July 11, 2014
				// To avoid addition of close base pairs we check for the following here
				if (j>=(i+TURN+1)){
					// Hosna April 11, 2014
					// I think we have already added gamma2(j,i) when coming to PL, so there is no need to add it again here.
					// Hosna July 17, 2014
					// I am adding gamma2 back here to avoid single base pair band predictions
					temp = get_PfromL(i+1,j-1,k,l) + gamma2(j,i);
					if(temp < min_energy){
						min_energy = temp;
						best_row= 3;
					}
				}
			}

			switch (best_row)
			{
				case 1:
					insert_node(i,l,j,k,P_PLiloop);
					break;
				case 2:
					insert_node(i,l,j,k,P_PLmloop);
					break;
				case 3:
					insert_node(i+1,l,j-1,k,P_PfromL);
					// Hosna, Feb 18, 2014
					// filling the structure
					f[i].pair = j;
					f[j].pair = i;
					f[i].type = P_PL;
					f[j].type = P_PL;
					break;
			}
		}
		break;

		case P_PR:
		{
			cand_pos_t i = cur_interval->i;
			cand_pos_t l = cur_interval->j;
			cand_pos_t j = cur_interval->k;
			cand_pos_t k = cur_interval->l;
			if(debug) printf("At P_PR with i=%d,l=%d,j=%d,k=%d\n",i,l,j,k);

			if (!(i <= j && j < k-1 && k <= l)){
				std::cerr << "boder cases: This should not have happened!, P_PR" <<std::endl;
				exit(EXIT_FAILURE);
			}

			if (i<0 ||j<0 || k<0 || l<0 || i>=n || j>=n || k>= n || l>= n){
				std::cerr << "impossible cases: This should not have happened!, P_PR" << std::endl;
				exit(EXIT_FAILURE);
			}

			energy_t min_energy = INF,temp=INF;
			cand_pos_t best_row = -1;
			if (pair[S_[k]][S_[l]]>0){
				//branch 1
				temp = get_PRiloop(i,j,k,l);
				if(temp < min_energy){
					min_energy = temp;
					best_row = 1;
				}

				//branch 2
				// Hosna, April 11, 2014
				// we need to add a branch penalty for multiloop that spans a band
				temp = get_PRmloop(i,j,k,l)+ bp_penalty;
				if(temp < min_energy){
					min_energy = temp;
					best_row=2;
				}

				//branch 3
				// Hosna, July 11, 2014
				// To avoid addition of close base pairs we check for the following here
				if (l>=(k+TURN+1)){
					// Hosna April 11, 2014
					// I think we have already added gamma2(l,k) when coming to PR, so there is no need to add it again here.
					// Hosna July 17, 2014
					// I am adding gamma2 back here to avoid single base pair band predictions
					temp = get_PfromR(i,j,k+1,l-1) + gamma2(l,k);
					if(temp < min_energy){
						min_energy = temp;
						best_row= 3;
					}
				}
			}

			switch (best_row)
			{
				case 1:
					insert_node(i,l,j,k,P_PRiloop);
					break;
				case 2:
					insert_node(i,l,j,k,P_PRmloop);
					break;
				case 3:
					insert_node(i,l-1,j,k+1,P_PfromR);
					// Hosna, Feb 18, 2014
					f[k].pair = l;
					f[l].pair = k;
					f[k].type = P_PR;
					f[l].type = P_PR;

					break;
			}

		}
		break;

		case P_PM:
		{
			cand_pos_t i = cur_interval->i;
			cand_pos_t l = cur_interval->j;
			cand_pos_t j = cur_interval->k;
			cand_pos_t k = cur_interval->l;
			if(debug) printf("At P_PM with i=%d,l=%d,j=%d,k=%d\n",i,l,j,k);

			if (!(i <= j && j < k-1 && k <= l)){
				std::cerr << "border cases: This should not have happened!, P_PM" << std::endl;
				exit(EXIT_FAILURE);
			}

			if (i<=0 ||j<=0 || k<=0 || l<=0 || i>n || j>n || k> n || l> n){
				std::cerr << "impossible cases: This should not have happened!, P_PM" << std::endl;
				exit(EXIT_FAILURE);
			}

			if (i==j && k ==l){
				// Hosna, Feb 25, 2014
				f[j].pair = k;
				f[k].pair = j;
				f[j].type = P_PM;
				f[k].type = P_PM;
				return;
			}

			energy_t min_energy = INF,temp=INF;
			cand_pos_t best_row = -1;

			if (pair[S_[j]][S_[k]]>0){
				//branch 1
				temp = get_PMiloop(i,j,k,l);
				if(temp < min_energy){
					min_energy = temp;
					best_row = 1;
				}
				//branch 2
				// Hosna, April 11, 2014
				// we need to add a branch penalty for multiloop that spans a band
				temp = get_PMmloop(i,j,k,l) + bp_penalty;
				if(temp < min_energy){
					min_energy = temp;
					best_row=2;
				}
				//branch 3
				// Hosna, July 11, 2014
				// To avoid addition of close base pairs we check for the following here
				if (k>=(j+TURN-1)){
					// Hosna April 11, 2014
					// I think we have already added gamma2(j,k) when coming to PM, so there is no need to add it again here.
					// Hosna July 17, 2014
					// I am adding gamma2 back here to avoid single base pair band predictions
					temp = get_PfromM(i,j-1,k+1,l) + gamma2(j,k);
					if(temp < min_energy){
						min_energy = temp;
						best_row= 3;
					}
				}
			}

			switch (best_row)
			{
				case 1:
					insert_node(i,l,j,k,P_PMiloop);
					break;
				case 2:
					insert_node(i,l,j,k,P_PMmloop);
					break;
				case 3:
					insert_node(i,l,j-1,k+1,P_PfromM);
					// Hosna, Feb 18, 2014
					f[j].pair = k;
					f[k].pair = j;
					f[j].type = P_PM;
					f[k].type = P_PM;
					break;
			}
		}
		break;

		case P_PO:
		{
			cand_pos_t i = cur_interval->i;
			cand_pos_t l = cur_interval->j;
			cand_pos_t j = cur_interval->k;
			cand_pos_t k = cur_interval->l;
			if(debug) printf("At P_PO with i=%d,l=%d,j=%d,k=%d\n",i,l,j,k);

			if (!(i <= j && j < k-1 && k <= l)){
				std::cerr << "border cases: This should not have happened!, P_PO" << std::endl;
				exit(EXIT_FAILURE);
			}
			// Hosna, April 3, 2014
			// adding impossible cases
			if (i<=0 ||j<=0 || k<=0 || l<=0 || i>n || j>n || k> n || l> n){
				std::cerr << "impossible cases: This should not have happened!, P_PO" << std::endl;
				exit(EXIT_FAILURE);
			}

			energy_t min_energy = INF,temp=INF;
			cand_pos_t best_row = -1;
			if (pair[S_[i]][S_[l]] > 0){
				//branch 1
				temp = get_POiloop(i,j,k,l);
				if(temp < min_energy){
					min_energy = temp;
					best_row = 1;
				}
				//branch 2
				// Hosna, April 11, 2014
				// we need to add a branch penalty for multiloop that spans a band
				temp = get_POmloop(i,j,k,l)+bp_penalty;
				if(temp < min_energy){
					min_energy = temp;
					best_row=2;
				}
				//branch 3
				// Hosna, July 11, 2014
				// To avoid addition of close base pairs we check for the following here
				if (l>=(i+TURN+1)){
					// Hosna April 11, 2014
					// I think we have already added gamma2(l,i) when coming to PO, so there is no need to add it again here.
					// Hosna July 17, 2014
					// I am adding gamma2 back here to avoid single base pair band predictions
					temp = get_PfromO(i+1,j,k,l-1) + gamma2(l,i);
					if(temp < min_energy){
						min_energy = temp;
						best_row= 3;
					}
				}
			}
			switch (best_row)
			{
				case 1:
					insert_node(i,l,j,k,P_POiloop);
					break;
				case 2:
					insert_node(i,l,j,k,P_POmloop);
					break;
				case 3:
					insert_node(i+1,l-1,j,k,P_PfromO);
					// Hosna, Feb 18, 2014
					f[i].pair = l;
					f[l].pair = i;
					f[i].type = P_PO;
					f[l].type = P_PO;
					break;
			}

		}
		break;

		case P_PfromL:
		{
			cand_pos_t i = cur_interval->i;
			cand_pos_t l = cur_interval->j;
			cand_pos_t j = cur_interval->k;
			cand_pos_t k = cur_interval->l;
			if(debug) printf("At P_PfromL with i=%d,l=%d,j=%d,k=%d\n",i,l,j,k);

			if (!(i <= j && j < k-1 && k <= l)){
				std::cerr << "This should not have happened!, P_PfromL" << std::endl;
				exit(EXIT_FAILURE);
			}

			if (i<=0 ||j<=0 || k<=0 || l<=0 || i>n || j>n || k> n || l> n){
				std::cerr << "This should not have happened!, P_PfromL" << std::endl;
				exit(EXIT_FAILURE);
			}

			if (i==j && k==l){
				return;
			}

			energy_t min_energy = INF,temp=INF;
			cand_pos_t best_row = -1,best_d=-1;

			for(cand_pos_t d=i+1; d< j; d++){
				//branch 1
				temp=get_PfromL(d,j,k,l)+get_WP(i,d-1);
				if(temp < min_energy){
					min_energy=temp;
					best_row=1;
					best_d = d;

				}
				//branch 2
				temp=get_PfromL(i,d,k,l)+get_WP(d+1,j);
				if(temp < min_energy){
					min_energy=temp;
					best_row=2;
					best_d = d;
				}
			}
			// branch 3
			//Hosna, July 28, 2014
			// I think going from PfromX to PX we are changing bands and so should be paying a band penalty
			temp = get_PR(i,j,k,l) + gamma2(l,k)+PB_penalty;
			if(temp < min_energy){
				min_energy = temp;
				best_row=3;
				best_d = -1;
			}

			//branch 4
			//Hosna, July 28, 2014
			// I think going from PfromX to PX we are changing bands and so should be paying a band penalty
			temp = get_PM(i,j,k,l) + gamma2(j,k)+PB_penalty;
			if(temp < min_energy){
				min_energy = temp;
				best_row=4;
				best_d = -1;
			}

			// branch 5
			//Hosna, July 28, 2014
			// I think going from PfromX to PX we are changing bands and so should be paying a band penalty
			temp = get_PO(i,j,k,l) + gamma2(l,i)+PB_penalty;
			if(temp < min_energy){
				min_energy = temp;
				best_row=5;
				best_d = -1;
			}
			switch (best_row)
			{
				case 1:
					if (best_d > -1){
						insert_node(best_d,l,j,k,P_PfromL);
						insert_node(i,best_d-1,P_WP);
					}
					break;
				case 2:
					if (best_d > -1){
						insert_node(i,l,best_d,k,P_PfromL);
						insert_node(best_d+1,j,P_WP);
					}
					break;
				case 3:
					insert_node(i,l,j,k,P_PR);
					break;
				case 4:
					insert_node(i,l,j,k,P_PM);
					break;
				case 5:
					insert_node(i,l,j,k,P_PO);
					break;
			}


		}
			break;

		case P_PfromR:
		{
			cand_pos_t i = cur_interval->i;
			cand_pos_t l = cur_interval->j;
			cand_pos_t j = cur_interval->k;
			cand_pos_t k = cur_interval->l;
			if(debug) printf("At P_PfromR with i=%d,l=%d,j=%d,k=%d\n",i,l,j,k);

			if (!(i <= j && j < k-1 && k <= l)){
				std::cerr << "This should not have happened!, P_PfromR" << std::endl;
				exit(EXIT_FAILURE);
			}

			if (i<=0 ||j<=0 || k<=0 || l<=0 || i>n || j>n || k> n || l> n){
				std::cerr << "impossible case: This should not have happened!, P_PfromR" << std::endl;
				exit(EXIT_FAILURE);
			}

			if (i==j && k==l){
				return;
			}

			energy_t min_energy = INF,temp=INF;
			cand_pos_t best_row = -1,best_d=-1;

			for(cand_pos_t d=k+1; d< l; d++){
				//branch 1
				temp=get_PfromR(i,j,d,l)+get_WP(k,d-1);
				if(temp < min_energy){
					min_energy = temp;
					best_row=1;
					best_d = d;

				}
				//branch 2
				temp=get_PfromR(i,j,k,d)+get_WP(d+1,l);
				if(temp < min_energy){
					min_energy=temp;
					best_row=2;
					best_d = d;
				}
			}

			//branch 3
			//Hosna, July 28, 2014
			// I think going from PfromX to PX we are changing bands and so should be paying a band penalty
			temp = get_PM(i,j,k,l) + gamma2(j,k)+PB_penalty;
			if(temp < min_energy){
				min_energy = temp;
				best_row=3;
				best_d = -1;
			}

			// branch 4
			//Hosna, July 28, 2014
			// I think going from PfromX to PX we are changing bands and so should be paying a band penalty
			temp = get_PO(i,j,k,l) + gamma2(l,i)+PB_penalty;
			if(temp < min_energy){
				min_energy = temp;
				best_row=4;
				best_d = -1;
			}
			switch (best_row)
			{
				case 1:
					if (best_d > -1){
						insert_node(i,l,j,best_d,P_PfromR);
						insert_node(k,best_d-1,P_WP);
					}
					break;
				case 2:
					if (best_d > -1){
						insert_node(i,best_d,j,k,P_PfromR);
						insert_node(best_d+1,l,P_WP);
					}
					break;
				case 3:
					insert_node(i,l,j,k,P_PM);
					break;
				case 4:
					insert_node(i,l,j,k,P_PO);
					break;
			}


		}
		break;

		case P_PfromM:
		{
			cand_pos_t i = cur_interval->i;
			cand_pos_t l = cur_interval->j;
			cand_pos_t j = cur_interval->k;
			cand_pos_t k = cur_interval->l;
			if(debug) printf("At P_PfromM with i=%d,l=%d,j=%d,k=%d\n",i,l,j,k);

			if (!(i <= j && j < k-1 && k <= l)){
				std::cerr << "This should not have happened!, P_PfromM" << std::endl;
				exit(EXIT_FAILURE);
			}

			if (i<=0 ||j<=0 || k<=0 || l<=0 || i>n || j>n || k> n || l> n){
				std::cerr << "This should not have happened!, P_PfromM" << std::endl;
				exit(EXIT_FAILURE);
			}

			if (i==j && k==l){
				return;
			}

			energy_t min_energy = INF,temp=INF;
			cand_pos_t best_row = -1,best_d=-1;

			for(cand_pos_t d=i+1; d< j; d++){
				//branch 1
				temp=get_PfromM(i,d,k,l)+get_WP(d+1,j);
				if(temp < min_energy){
					min_energy=temp;
					best_row=1;
					best_d = d;

				}
			}
			for(cand_pos_t d=k+1; d< l; d++){
				//branch 2
				temp=get_PfromM(i,j,d,l)+get_WP(k,d-1);
				if(temp < min_energy){
					min_energy=temp;
					best_row=2;
					best_d = d;
				}
			}
			// branch 3
			//Hosna, July 28, 2014
			// I think going from PfromX to PX we are changing bands and so should be paying a band penalty
			temp = get_PL(i,j,k,l) + gamma2(j,i)+PB_penalty;
			if(temp < min_energy){
				min_energy = temp;
				best_row=3;
				best_d = -1;
			}

			//branch 4
			//Hosna, July 28, 2014
			// I think going from PfromX to PX we are changing bands and so should be paying a band penalty
			temp = get_PR(i,j,k,l) + gamma2(l,k)+PB_penalty;
			if(temp < min_energy){
				min_energy = temp;
				best_row=4;
				best_d = -1;
			}

			switch (best_row)
			{
				case 1:
					if (best_d > -1){
						insert_node(i,l,best_d,k,P_PfromM);
						insert_node(best_d+1,j,P_WP);
					}
					break;
				case 2:
					if (best_d > -1){
						insert_node(i,l,j,best_d,P_PfromM);
						insert_node(k,best_d-1,P_WP);
					}
					break;
				case 3:
					insert_node(i,l,j,k,P_PL);
					break;
				case 4:
					insert_node(i,l,j,k,P_PR);
					break;
			}


		}
			break;

		case P_PfromO:
		{
			cand_pos_t i = cur_interval->i;
			cand_pos_t l = cur_interval->j;
			cand_pos_t j = cur_interval->k;
			cand_pos_t k = cur_interval->l;
			if(debug) printf("At P_PfromO with i=%d,l=%d,j=%d,k=%d\n",i,l,j,k);

			if (!(i <= j && j < k-1 && k <= l)){
				std::cerr << "border cases: This should not have happened!, P_PfromO" << std::endl;
				exit(EXIT_FAILURE);
			}

			if (i<=0 ||j<=0 || k<=0 || l<=0 || i>n || j>n || k> n || l> n){
				std::cerr << "impossible case: This should not have happened!, P_PfromO" << std::endl;
				exit(EXIT_FAILURE);
			}

			if (i==j && k==l){
				return;
			}

			energy_t min_energy = INF,temp=INF;
			cand_pos_t best_row = -1,best_d=-1;

			for(cand_pos_t d=i+1; d< j; d++){
				//branch 1
				temp=get_PfromO(d,j,k,l)+get_WP(i,d-1);
				if(temp < min_energy){
					min_energy=temp;
					best_row=1;
					best_d = d;

				}
			}
			for(cand_pos_t d=k+1; d< l; d++){
				//branch 2
				temp=get_PfromO(i,j,k,d)+get_WP(d+1,l);
				if(temp < min_energy){
					min_energy=temp;
					best_row=2;
					best_d = d;
				}
			}
			// branch 3
			//Hosna, July 28, 2014
			// I think going from PfromX to PX we are changing bands and so should be paying a band penalty
			temp = get_PL(i,j,k,l) + gamma2(j,i)+PB_penalty;
			if(temp < min_energy){
				min_energy = temp;
				best_row=3;
				best_d = -1;
			}

			//branch 4
			//Hosna, July 28, 2014
			// I think going from PfromX to PX we are changing bands and so should be paying a band penalty
			temp = get_PR(i,j,k,l) + gamma2(l,k)+PB_penalty;
			if(temp < min_energy){
				min_energy = temp;
				best_row=4;
				best_d = -1;
			}

			switch (best_row)
			{
				case 1:
					if (best_d > -1){
						insert_node(best_d,l,j,k,P_PfromO);
						insert_node(i,best_d-1,P_WP);
					}
					break;
				case 2:
					if (best_d > -1){
						insert_node(i,best_d,j,k,P_PfromO);
						insert_node(best_d+1,l,P_WP);
					}
					break;
				case 3:
					insert_node(i,l,j,k,P_PL);
					break;
				case 4:
					insert_node(i,l,j,k,P_PR);
					break;
			}


		}
			break;
		case P_WB:
		{
			cand_pos_t i = cur_interval->i;
			cand_pos_t l = cur_interval->j;
			if(debug) printf("At P_WM with i=%d,l=%d\n",i,l);

			if (i<=0 ||l<=0 || i>n || l >n){
				std::cerr << "impossible cases: This should not have happened!, P_WB" << std::endl;
				exit(EXIT_FAILURE);
			}

			if (i>l){
				return;
			}

			energy_t min_energy = INF,temp=INF;
			cand_pos_t best_row = -1;
			//branch 1
			temp = get_WBP(i,l);
			if (temp < min_energy){
				min_energy = temp;
				best_row=1;
			}

			temp = beta1P*(l-i+1);
			if (temp < min_energy){
				min_energy = temp;
				best_row=2;
			}

			switch (best_row)
			{
				case 1:
					insert_node(i,l,P_WBP);
					break;
				case 2:
					// do nothing.
					break;
			}

		}
			break;
		case P_WBP:
		{
			cand_pos_t i = cur_interval->i;
			cand_pos_t l = cur_interval->j;
			if(debug) printf("At P_WBP with i=%d,l=%d\n",i,l);
			if (i>l){
				std::cerr << "border case: This should not have happened!, P_WBP" << std::endl;
				exit(EXIT_FAILURE);
			}

			if (i<=0 ||l<=0 || i>n || l >n){
				std::cerr << "impossible cases: This should not have happened!, P_WBP" << std::endl;
				exit(EXIT_FAILURE);
			}

			energy_t min_energy = INF,temp=INF;
			cand_pos_t best_row = -1,best_d=-1,best_e=-1;

			for(cand_pos_t d=i; d< l; d++){
				for(cand_pos_t e = d+1; e<= l; e++){
					//branch 1
					// Hosna, August 26, 2014
					// comparing calculation of WI in HFold and WPP in CCJ, I found that
					// in HFold we add another penalty called PPS_penalty for closed regions inside a pseudoloop or multiloop that spans a band
					//int common = get_WB(i,d-1) + beta1P*(l-e);
					energy_t common = get_WB(i,d-1) + beta1P*(l-e)+PPS_penalty;
					temp = common + V->get_energy(d,e) +beta2P(e,d);
					if (temp < min_energy){
						min_energy = temp;
						best_row = 1;
						best_d = d;
						best_e = e;
					}

					//branch 2
					temp = common + get_P(d,e) + gamma0m;
					if (temp < min_energy){
						min_energy = temp;
						best_row = 2;
						best_d = d;
						best_e = e;
					}
				}
			}

			switch (best_row)
			{
				case 1:
					insert_node(i,best_d-1,P_WB);
					insert_node(best_d,best_e,LOOP);
					break;
				case 2:
					insert_node(i,best_d-1,P_WB);
					insert_node(best_d,best_e,P_P);

					break;
			}
		}
			break;

		case P_WP:
		{
			cand_pos_t i = cur_interval->i;
			cand_pos_t l = cur_interval->j;
			if(debug) printf("At P_WP with i=%d,l=%d\n",i,l);

			if (i<=0 ||l<=0 || i>n || l >n){
				std::cerr << "impossible cases: This should not have happened!, P_WP" << std::endl;
				exit(EXIT_FAILURE);
			}

			if (i>l){
				return;
			}

			energy_t min_energy = INF,temp=INF;
			cand_pos_t best_row = -1;
			//branch 1
			temp = get_WPP(i,l);
			if (temp < min_energy){
				min_energy = temp;
				best_row=1;
			}

			temp = gamma1*(l-i+1);
			if (temp < min_energy){
				min_energy = temp;
				best_row=2;
			}

			switch (best_row)
			{
				case 1:
					insert_node(i,l,P_WPP);
					break;
				case 2:
					// do nothing.
					break;
			}

		}
			break;
		case P_WPP:
		{//TODO: change WPP
			cand_pos_t i = cur_interval->i;
			cand_pos_t l = cur_interval->j;
			if(debug) printf("At P_WPP with i=%d,l=%d\n",i,l);
			if (i>l){
				std::cerr << "border case: This should not have happened!, P_WPP" << std::endl;
				exit(EXIT_FAILURE);
			}

			if (i<=0 ||l<=0 || i>n || l >n){
				std::cerr << "impossible cases: This should not have happened!, P_WPP" << std::endl;
				exit(EXIT_FAILURE);
			}

			energy_t min_energy = INF,temp=INF;
			cand_pos_t best_row = -1,best_d=-1,best_e=-1;

			for(cand_pos_t d=i; d< l; d++){
				for(cand_pos_t e = d+1; e<= l; e++){
					// Hosna, August 26, 2014
					// comparing calculation of WI in HFold and WPP in CCJ, I found that
					// in HFold we add another penalty called PPS_penalty for closed regions inside a pseudoloop or multiloop that spans a band
					//int common = get_WP(i,d-1) + gamma1*(l-e);
					energy_t common = get_WP(i,d-1) + gamma1*(l-e)+PPS_penalty;
					//branch 1
					temp = V->get_energy(d,e) + gamma2(e,d) + common;
					if (temp < min_energy){
						min_energy = temp;
						best_row = 1;
						best_d = d;
						best_e = e;
					}

					//branch 2
					temp = get_P(d,e) + gamma0P + common;
					if (temp < min_energy){
						min_energy = temp;
						best_row = 2;
						best_d = d;
						best_e = e;
					}
				}
			}

			switch (best_row)
			{
				case 1:
					insert_node(i,best_d-1,P_WP);
					insert_node(best_d,best_e,LOOP);
					break;
				case 2:
					insert_node(i,best_d-1,P_WP);
					insert_node(best_d,best_e,P_P);
					break;
			}
		}
			break;
		case P_PLiloop:
		{
			// changing gapped region borders to what we mean in recurrences
			cand_pos_t i = cur_interval->i;
			cand_pos_t l = cur_interval->j;
			cand_pos_t j = cur_interval->k;
			cand_pos_t k = cur_interval->l;
			if(debug) printf("At P_PLiloop with i=%d,l=%d,j=%d,k=%d\n",i,l,j,k);

			if (!(i < j && j < k-1 && k < l)){
				std::cerr << "border cases: This should not have happened!, P_PLiloop" << std::endl;
				exit(EXIT_FAILURE);
			}

			if (i<=0 || j<=0 || k<=0 || l<=0 || i>n || j>n || k> n || l> n){
				std::cerr << "impossbible cases: This should not have happened!, P_PLiloop" << std::endl;
				exit(EXIT_FAILURE);
			}

			// Hosna, Feb 25, 2014
			f[i].pair = j;
			f[j].pair = i;
			f[i].type = P_PLiloop;
			f[j].type = P_PLiloop;

			energy_t min_energy = INF,temp=INF;
			cand_pos_t best_row = -1,best_d=-1,best_dp=-1;
			if (pair[S_[i]][S_[j]] > 0){
				//branch 1
				temp = get_PL(i+1,j-1,k,l) + get_e_stP(i,j);
				if (temp < min_energy){
					min_energy = temp;
					best_row=1;
				}
				//branch 2
				for(cand_pos_t d= i+1; d<std::min(j,i+MAXLOOP); ++d){
					for(cand_pos_t dp = j-1; dp > std::max(d+TURN,j-MAXLOOP); --dp){
						temp = get_e_intP(i,d,dp,j) + get_PL(d,dp,k,l);
						if(temp < min_energy){
							min_energy = temp;
							best_d = d;
							best_dp = dp;
							best_row = 2;
						}
					}
            }	
				
			}

			switch (best_row)
			{
				case 1:
					insert_node(i+1,l,j-1,k,P_PL);
					break;
				case 2:
					insert_node(best_d,l,best_dp,k,P_PL);
					break;
			}
		}
			break;

		case P_PLmloop:
		{
			// changing gapped region borders to what we mean in recurrences
			cand_pos_t i = cur_interval->i;
			cand_pos_t l = cur_interval->j;
			cand_pos_t j = cur_interval->k;
			cand_pos_t k = cur_interval->l;
			if(debug) printf("At P_PLmloop with i=%d,l=%d,j=%d,k=%d\n",i,l,j,k);

			if (!(i <= j && j < k-1 && k <= l)){
				std::cerr << "border cases: This should not have happened!, P_PLmloop" << std::endl;
				exit(EXIT_FAILURE);
			}

			if (i<=0 || j<=0 || k<=0 || l<=0 || i>n || j>n || k> n || l> n){
				std::cerr << "impossible cases: This should not have happened!, P_PLmloop" << std::endl;
				exit(EXIT_FAILURE);
			}

			// Hosna, Feb 25, 2014
			f[i].pair = j;
			f[j].pair = i;
			f[i].type = P_PLmloop;
			f[j].type = P_PLmloop;


			energy_t min_energy = INF,temp=INF;
			cand_pos_t best_row = -1, best_d=-1;
			for(cand_pos_t d = i+1; d < j-1; d++){
				//branch 1
				// Hosna, Feb 23, 2014
				// Since PLmloop comes from PL and in PL i and j must pair, so we must have the same thing here
				// therefore we should have PLmloop0(d,j-1,k,l)
				temp = get_PLmloop0(d,j-1,k,l) + get_WBP(i+1,d-1) + beta0P + beta2P(j,i);
				if (temp < min_energy){
					min_energy = temp;
					best_row = 1;
					best_d = d;
				}
				//branch 2
				temp = get_PLmloop1(d,j-1,k,l) + get_WB(i+1,d-1) + beta0P + beta2P(j,i);
				if (temp < min_energy){
					min_energy = temp;
					best_row = 2;
					best_d = d;
				}
			}
			switch (best_row)
			{
				case 1:
					insert_node(best_d,l,j-1,k,P_PLmloop0);
					insert_node(i+1,best_d-1,P_WBP);
					break;
				case 2:
					insert_node(best_d,l,j-1,k,P_PLmloop1);
					insert_node(i+1,best_d-1,P_WB);
					break;
			}

		}
			break;
		case P_PLmloop0:
		{
			// changing gapped region borders to what we mean in recurrences
			cand_pos_t i = cur_interval->i;
			cand_pos_t l = cur_interval->j;
			cand_pos_t j = cur_interval->k;
			cand_pos_t k = cur_interval->l;
			if(debug) printf("At P_PLmloop0 with i=%d,l=%d,j=%d,k=%d\n",i,l,j,k);

			if (!(i <= j && j < k-1 && k <= l)){
				std::cerr << "border cases: This should not have happened!, P_PLmloop0" << std::endl;
				exit(EXIT_FAILURE);
			}

			if (i<=0 || j<=0 || k<=0 || l<=0 || i>n || j>n || k> n || l> n){
				std::cerr << "impossible cases: This should not have happened!, P_PLmloop0" << std::endl;
				exit(EXIT_FAILURE);
			}

			energy_t min_energy = INF,temp=INF;
			cand_pos_t best_d=-1;
			for(cand_pos_t d = i+1; d < j; d++){
				// Hosna, feb 23, 2014
				// changed the recurrences so that j-1 is accounted for in PLmloop
				temp = get_PL(i,d,k,l) + get_WB(d+1,j) + beta2P(d,i);
				if (temp < min_energy){
					min_energy = temp;
					best_d = d;
				}
			}
			insert_node(i,l,best_d,k,P_PL);
			insert_node(best_d+1,j,P_WB);

		}
			break;
		case P_PLmloop1:
		{
			// changing gapped region borders to what we mean in recurrences
			cand_pos_t i = cur_interval->i;
			cand_pos_t l = cur_interval->j;
			cand_pos_t j = cur_interval->k;
			cand_pos_t k = cur_interval->l;
			if(debug) printf("At P_PLmloop1 with i=%d,l=%d,j=%d,k=%d\n",i,l,j,k);

			if (!(i <= j && j < k-1 && k <= l)){
				std::cerr<< "border cases: This should not have happened!, P_PLmloop1" << std::endl;
				exit(EXIT_FAILURE);
			}

			if (i<=0 || j<=0 || k<=0 || l<=0 || i>n || j>n || k> n || l> n){
				std::cerr<< "impossible cases: This should not have happened!, P_PLmloop1" << std::endl;
				exit(EXIT_FAILURE);
			}

			energy_t min_energy = INF,temp=INF;
			cand_pos_t best_d=-1;
			for(cand_pos_t d = i+1; d < j; d++){
				// Hosna, feb 23, 2014
				// changed the recurrences so that j-1 is accounted for in PLmloop
				temp = get_PL(i,d,k,l) + get_WBP(d+1,j) + beta2P(d,i);
				if (temp < min_energy){
					min_energy = temp;
					best_d = d;
				}
			}
			insert_node(i,l,best_d,k,P_PL);
			insert_node(best_d+1,j,P_WBP);
		}
			break;

		case P_PRiloop:
		{
			// changing gapped region borders to what we mean in recurrences
			cand_pos_t i = cur_interval->i;
			cand_pos_t l = cur_interval->j;
			cand_pos_t j = cur_interval->k;
			cand_pos_t k = cur_interval->l;
			if(debug) printf("At P_PRiloop with i=%d,l=%d,j=%d,k=%d\n",i,l,j,k);

			if (!(i <= j && j < k-1 && k <= l)){
				std::cerr << "border cases: This should not have happened!, P_PRiloop" << std::endl;
				exit(EXIT_FAILURE);
			}

			if (i<=0 || j<=0 || k<=0 || l<=0 || i>n || j>n || k> n || l> n){
				std::cerr << "impossible cases: This should not have happened!, P_PRiloop" << std::endl;
				exit(EXIT_FAILURE);
			}

			// Hosna, Feb 25, 2014
			f[k].pair = l;
			f[l].pair = k;
			f[k].type = P_PRiloop;
			f[l].type = P_PRiloop;

			energy_t min_energy = INF,temp=INF;
			cand_pos_t best_row=-1,best_d=-1,best_dp=-1;
			if (pair[S_[k]][S_[l]]>0){
				//branch 1
				temp = get_PR(i,j,k+1,l-1) + get_e_stP(k,l);
				if(temp < min_energy){
					min_energy = temp;
					best_row = 1;
				}
				//branch 2
				for(cand_pos_t d= k+1; d<std::min(l,k+MAXLOOP); ++d){
					for(cand_pos_t dp=l-1; dp > std::max(d+TURN,l-MAXLOOP); --dp){
						temp = get_e_intP(k,d,dp,l) + get_PR(i,j,d,dp);
						if(temp < min_energy){
							min_energy = temp;
							best_d = d;
							best_dp = dp;
							best_row = 2;
						}
					}
				}
			}

			switch (best_row)
			{
				case 1:
					insert_node(i,l-1,j,k+1,P_PR);
					break;
				case 2:
					insert_node(i,best_dp,j,best_d,P_PR);
					break;
			}

		}
			break;

		case P_PRmloop:
		{
			// changing gapped region borders to what we mean in recurrences
			cand_pos_t i = cur_interval->i;
			cand_pos_t l = cur_interval->j;
			cand_pos_t j = cur_interval->k;
			cand_pos_t k = cur_interval->l;
			if(debug) printf("At P_PRmloop with i=%d,l=%d,j=%d,k=%d\n",i,l,j,k);

			if (!(i <= j && j < k-1 && k <= l)){
				std::cerr << "border cases: This should not have happened!, P_PRmloop" << std::endl;
				exit(EXIT_FAILURE);
			}
			if (i<=0 || j<=0 || k<=0 || l<=0 || i>n || j>n || k> n || l> n){
				std::cerr << "impossible cases: This should not have happened!, P_PRmloop" << std::endl;
				exit(EXIT_FAILURE);
			}

			// Hosna, Feb 25, 2014
			f[k].pair = l;
			f[l].pair = k;
			f[k].type = P_PRmloop;
			f[l].type = P_PRmloop;

			energy_t min_energy = INF,temp=INF;
			cand_pos_t best_row=-1, best_d=-1;
			//Hosna Feb 23, 2014
			// changed the recurrences to have l-1 instead of l in PRmloop and removed l-1 from PRmloop0,1, as we are accounting for k.l here
			for(cand_pos_t d = k+1; d < l-1; d++){
				// branch 1
				temp = get_PRmloop0(i,j,d,l-1) + get_WBP(k+1,d-1) + beta0P + beta2P(l,k);
				if (temp < min_energy){
					min_energy = temp;
					best_row = 1;
					best_d=d;
				}
				// branch 2
				temp = get_PRmloop1(i,j,d,l-1) + get_WB(k+1,d-1) + beta0P + beta2P(l,k);
				if (temp < min_energy){
					min_energy = temp;
					best_row=2;
					best_d = d;
				}
			}
			switch (best_row)
			{
				case 1:
					insert_node(i,l-1,j,best_d,P_PRmloop0);
					insert_node(k+1,best_d-1,P_WBP);
					break;
				case 2:
					insert_node(i,l-1,j,best_d,P_PRmloop1);
					insert_node(k+1,best_d-1,P_WB);
					break;
			}

		}
			break;

		case P_PRmloop0:
		{
			// changing gapped region borders to what we mean in recurrences
			cand_pos_t i = cur_interval->i;
			cand_pos_t l = cur_interval->j;
			cand_pos_t j = cur_interval->k;
			cand_pos_t k = cur_interval->l;
			if(debug) printf("At P_PRmloop0 with i=%d,l=%d,j=%d,k=%d\n",i,l,j,k);

			if (!(i <= j && j < k-1 && k <= l)){
				std::cerr << "border cases: This should not have happened!, P_PRmloop0" << std::endl;
				exit(EXIT_FAILURE);
			}

			if (i<=0 || j<=0 || k<=0 || l<=0 || i>n || j>n || k> n || l> n){
				std::cerr << "impossible cases: This should not have happened!, P_PRmloop0" << std::endl;
				exit(EXIT_FAILURE);
			}

			energy_t min_energy = INF,temp=INF;
			cand_pos_t best_d=-1;
			for(cand_pos_t d = k+1; d < l; d++){
				temp = get_PR(i,j,k,d) + get_WB(d+1,l) + beta2P(k,d);
				if (temp < min_energy){
					min_energy = temp;
					best_d=d;
				}
			}
			insert_node(i,best_d,j,k,P_PR);
			insert_node(best_d+1,l,P_WB);

		}
			break;
		case P_PRmloop1:
		{
			// changing gapped region borders to what we mean in recurrences
			cand_pos_t i = cur_interval->i;
			cand_pos_t l = cur_interval->j;
			cand_pos_t j = cur_interval->k;
			cand_pos_t k = cur_interval->l;
			if(debug) printf("At P_PRmloop1 with i=%d,l=%d,j=%d,k=%d\n",i,l,j,k);

			if (!(i <= j && j < k-1 && k <= l)){
				std::cerr << "border cases: This should not have happened!, P_PRmloop1" << std::endl;
				exit(EXIT_FAILURE);
			}

			if (i<=0 || j<=0 || k<=0 || l<=0 || i>n || j>n || k> n || l> n){
				std::cerr << "impossible cases: This should not have happened!, P_PRmloop1" << std::endl;
				exit(EXIT_FAILURE);
			}

			energy_t min_energy = INF,temp=INF; 
			cand_pos_t best_d=-1;
			for(cand_pos_t d = k+1; d < l; d++){
				temp = get_PR(i,j,k,d) + get_WBP(d+1,l) + beta2P(k,d);
				if (temp < min_energy){
					min_energy = temp;
					best_d=d;
				}
			}
			insert_node(i,best_d,j,k,P_PR);
			insert_node(best_d+1,l,P_WBP);

		}
			break;

		case P_PMiloop:
		{
			// changing gapped region borders to what we mean in recurrences
			cand_pos_t i = cur_interval->i;
			cand_pos_t l = cur_interval->j;
			cand_pos_t j = cur_interval->k;
			cand_pos_t k = cur_interval->l;
			if(debug) printf("At P_PMiloop with i=%d,l=%d,j=%d,k=%d\n",i,l,j,k);

			if (!(i <= j && j < k-1 && k <= l)){
				std::cerr << "border cases: This should not have happened!, P_PMiloop" << std::endl;
				exit(EXIT_FAILURE);
			}

			if (i<=0 || j<=0 || k<=0 || l<=0 || i>n || j>n || k> n || l> n){
				std::cerr << "impossible cases: This should not have happened!, P_PMiloop" << std::endl;
				exit(EXIT_FAILURE);
			}

			// Hosna, Feb 25, 2014
			f[j].pair = k;
			f[k].pair = j;
			f[j].type = P_PMiloop;
			f[k].type = P_PMiloop;
			
			energy_t min_energy = INF,temp=INF;
			cand_pos_t best_d=-1,best_dp=-1,best_row=-1;

			if (pair[S_[j]][S_[k]]>0){
					// branch 1
				temp = get_PM(i,j-1,k+1,l) + get_e_stP(j-1,k+1);
				if (temp < min_energy){
					min_energy = temp;
					best_row=1;
				}

				for(cand_pos_t d= j-1; d>std::max(i,j-MAXLOOP); --d){
					for (cand_pos_t dp=k+1; dp <std::min(l,k+MAXLOOP); ++dp) {
						temp = get_e_intP(d,j,k,dp) + get_PM(i,d,dp,l);
						if(temp < min_energy){
							min_energy = temp;
							best_d = d;
							best_dp = dp;
							best_row = 2;
						}
					}
				}
			}

			switch (best_row)
			{
				case 1:
					insert_node(i,l,j-1,k+1,P_PM);
					break;
				case 2:
					insert_node(i,l,best_d,best_dp,P_PM);
					break;
			}
		}
			break;

		case P_PMmloop:
		{
			// changing gapped region borders to what we mean in recurrences
			cand_pos_t i = cur_interval->i;
			cand_pos_t l = cur_interval->j;
			cand_pos_t j = cur_interval->k;
			cand_pos_t k = cur_interval->l;
			if(debug) printf("At P_PMmloop with i=%d,l=%d,j=%d,k=%d\n",i,l,j,k);

			if (!(i <= j && j < k-1 && k <= l)){
				std::cerr << "border cases: This should not have happened!, P_PMmloop" << std::endl;
				exit(EXIT_FAILURE);
			}

			if (i<=0 || j<=0 || k<=0 || l<=0 || i>n || j>n || k> n || l> n){
				std::cerr << "impossible cases: This should not have happened!, P_PMmloop" << std::endl;
				exit(EXIT_FAILURE);
			}

			// Hosna, Feb 25, 2014
			f[j].pair = k;
			f[k].pair = j;
			f[j].type = P_PMmloop;
			f[k].type = P_PMmloop;

			energy_t min_energy = INF,temp=INF;
			cand_pos_t best_row=-1, best_d=-1;
			for(cand_pos_t d = i+1; d < j; d++){
				// Hosna Feb 23, 2014
				// changed the recurrence of PMmloop, to have k+1 instead of k when calling PMmloop0,1, and removed k+1 from PMmloop0,1
				// as j.k pair is accounted for in this recurrence

				// branch 1
				temp = get_PMmloop0(i,d,k+1,l) + get_WBP(d+1,j-1) + beta0P + beta2P(j,k);
				if (temp < min_energy){
					min_energy = temp;
					best_row=1;
					best_d=d;
				}
				// branch 2
				temp = get_PMmloop1(i,d,k+1,l) + get_WB(d+1,j-1) + beta0P + beta2P(j,k);
				if (temp < min_energy){
					min_energy = temp;
					best_row=2;
					best_d=d;
				}
			}

			switch (best_row)
			{
				case 1:
					insert_node(i,l,best_d,k+1,P_PMmloop0);
					insert_node(best_d+1,j-1,P_WBP);
					break;
				case 2:
					insert_node(i,l,best_d,k+1,P_PMmloop1);
					insert_node(best_d+1,j-1,P_WB);
					break;
			}

		}
			break;

		case P_PMmloop0:
		{
			// changing gapped region borders to what we mean in recurrences
			cand_pos_t i = cur_interval->i;
			cand_pos_t l = cur_interval->j;
			cand_pos_t j = cur_interval->k;
			cand_pos_t k = cur_interval->l;
			if(debug) printf("At P_PMmloop0 with i=%d,l=%d,j=%d,k=%d\n",i,l,j,k);

			if (!(i <= j && j < k-1 && k <= l)){
				std::cerr << "border cases: This should not have happened!, P_PMmloop0" << std::endl;
				exit(EXIT_FAILURE);
			}

			if (i<=0 || j<=0 || k<=0 || l<=0 || i>n || j>n || k> n || l> n){
				std::cerr << "impossible cases: This should not have happened!, P_PMmloop0" << std::endl;
				exit(EXIT_FAILURE);
			}

			energy_t min_energy = INF,temp=INF;
			cand_pos_t best_d=-1;
			for(int d = k+1; d < l; d++){
				//temp = get_PM(i,j,d,l) + get_WB(k+1,d-1) + beta2P(j,d);
				temp = get_PM(i,j,d,l) + get_WB(k,d-1) + beta2P(j,d);
				if (temp < min_energy){
					min_energy = temp;
					best_d =d;
				}
			}
			insert_node(i,l,j,best_d,P_PM);
			insert_node(k,best_d-1,P_WB);

		}
			break;
		case P_PMmloop1:
		{
			// changing gapped region borders to what we mean in recurrences
			cand_pos_t i = cur_interval->i;
			cand_pos_t l = cur_interval->j;
			cand_pos_t j = cur_interval->k;
			cand_pos_t k = cur_interval->l;
			if(debug) printf("At P_PMmloop1 with i=%d,l=%d,j=%d,k=%d\n",i,l,j,k);

			if (!(i <= j && j < k-1 && k <= l)){
				std::cerr << "border cases: This should not have happened!, P_PMmloop1" << std::endl;
				exit(EXIT_FAILURE);
			}

			if (i<=0 || j<=0 || k<=0 || l<=0 || i>n || j>n || k> n || l> n){
				std::cerr << "impossible cases: This should not have happened!, P_PMmloop1" << std::endl;
				exit(EXIT_FAILURE);
			}

			energy_t min_energy = INF,temp=INF;
			cand_pos_t best_d=-1;
			for(cand_pos_t d = k+1; d < l; d++){
				temp = get_PM(i,j,d,l) + get_WBP(k,d-1) + beta2P(j,d);
				if (temp < min_energy){
					min_energy = temp;
					best_d=d;
				}
			}
			insert_node(i,l,j,best_d,P_PM);
			insert_node(k,best_d-1,P_WBP);

		}
			break;
		case P_POiloop:
		{
			// changing gapped region borders to what we mean in recurrences
			cand_pos_t i = cur_interval->i;
			cand_pos_t l = cur_interval->j;
			cand_pos_t j = cur_interval->k;
			cand_pos_t k = cur_interval->l;
			if(debug) printf("At P_POiloop with i=%d,l=%d,j=%d,k=%d\n",i,l,j,k);

			if (i<=0 || j<=0 || k<=0 || l<=0 || i>n || j>n || k> n || l> n){
				std::cerr << "impossible cases: This should not have happened!, P_POiloop" << std::endl;
				exit(EXIT_FAILURE);
			}

			if (!(i <= j && j < k-1 && k <= l)){
				std::cerr << "border cases: This should not have happened!, P_POiloop" << std::endl;
				exit(EXIT_FAILURE);
			}

			// Hosna, Feb 25, 2014
			f[i].pair = l;
			f[l].pair = i;
			f[i].type = P_POiloop;
			f[l].type = P_POiloop;


			energy_t min_energy = INF;
			cand_pos_t best_d=-1,best_dp=-1,best_row=-1;
			if (pair[S_[i]][S_[l]]>0){
				//branch 1
				energy_t temp = get_PO(i+1,j,k,l-1) + get_e_stP(i,l);
				if(temp < min_energy){
					min_energy = temp;
					best_row=1;
				}
				for(cand_pos_t d= i+1; d<std::min(j,i+MAXLOOP); ++d){
					for (cand_pos_t dp=l-1; dp >std::max(l-MAXLOOP,k); --dp) {
						energy_t branch2 = get_e_intP(i,d,dp,l) + get_PO(d,j,dp,k);
						if(branch2 < min_energy){
							min_energy = branch2;
							best_row = 2;
							best_d = d;
							best_dp = dp;
						}
					}

    			}
			}

			switch (best_row)
			{
				case 1:
					insert_node(i+1,l-1,j,k,P_PO);
					break;
				case 2:
					insert_node(best_d,k,j,best_dp,P_PO);
					break;
			}

		}
			break;

		case P_POmloop:
		{
			// changing gapped region borders to match the recurrences
			cand_pos_t i = cur_interval->i;
			cand_pos_t l = cur_interval->j;
			cand_pos_t j = cur_interval->k;
			cand_pos_t k = cur_interval->l;
			if(debug) printf("At P_POmloop with i=%d,l=%d,j=%d,k=%d\n",i,l,j,k);

			if (!(i <= j && j < k-1 && k <= l)){
				std::cerr << "border cases: This should not have happened!, P_POmloop" << std::endl;
				exit(EXIT_FAILURE);
			}

			if (i<=0 || j<=0 || k<=0 || l<=0 || i>n || j>n || k> n || l> n){
				std::cerr << "impossible cases: This should not have happened!, P_POmloop" << std::endl;
				exit(EXIT_FAILURE);
			}
			// Hosna, Feb 25, 2014
			f[i].pair = l;
			f[l].pair = i;
			f[i].type = P_POmloop;
			f[l].type = P_POmloop;


			energy_t min_energy = INF,temp=INF;
			cand_pos_t best_row=-1,best_d=-1;
			// Hosna Feb 23, 2014
			// changed recurrences for POmloop to have l-1 instead of l and removed l-1 from POmloop0,1 as i.l is accounted for here in this recurrence
			for(cand_pos_t d = i+1; d < j; d++){
				//branch 1
				temp = get_POmloop0(d,j,k,l-1) + get_WBP(i+1,d-1) + beta0P + beta2P(l,i);
				if (temp < min_energy){
					min_energy = temp;
					best_row=1;
					best_d=d;
				}
				//branch 2
				temp = get_POmloop1(d,j,k,l-1) + get_WB(i+1,d-1) + beta0P + beta2P(l,i);
				if (temp < min_energy){
					min_energy = temp;
					best_row=2;
					best_d=d;
				}
			}

			switch (best_row)
			{
				case 1:
					insert_node(best_d,l-1,j,k,P_POmloop0);
					insert_node(i+1,best_d-1,P_WBP);
					break;
				case 2:
					insert_node(best_d,l-1,j,k,P_POmloop1);
					insert_node(i+1,best_d-1,P_WB);
					break;
			}
		}
			break;

		case P_POmloop0:
		{
			// changing gapped region borders to match the recurrences
			cand_pos_t i = cur_interval->i;
			cand_pos_t l = cur_interval->j;
			cand_pos_t j = cur_interval->k;
			cand_pos_t k = cur_interval->l;
			if(debug) printf("At P_POmloop0 with i=%d,l=%d,j=%d,k=%d\n",i,l,j,k);

			if (!(i <= j && j < k-1 && k <= l)){
				std::cerr << "border cases: This should not have happened!, P_POmloop0" << std::endl;
				exit(EXIT_FAILURE);
			}

			if (i<=0 || j<=0 || k<=0 || l<=0 || i>n || j>n || k> n || l> n){
				std::cerr << "impossible cases: This should not have happened!, P_POmloop0" << std::endl;
				exit(EXIT_FAILURE);
			}

			energy_t min_energy = INF,temp=INF;
			cand_pos_t best_d=-1;
			for(int d = k+1; d < l; d++){
				temp = get_PO(i,j,k,d) + get_WB(d+1,l) + beta2P(d,i);
				if (temp < min_energy){
					min_energy = temp;
					best_d=d;
				}
			}
			insert_node(i,best_d,j,k,P_PO);
			insert_node(best_d+1,l,P_WB);
		}
			break;
		case P_POmloop1:
		{
			// changing gapped region borders to match the recurrences
			cand_pos_t i = cur_interval->i;
			cand_pos_t l = cur_interval->j;
			cand_pos_t j = cur_interval->k;
			cand_pos_t k = cur_interval->l;
			if(debug) printf("At P_POmloop1 with i=%d,l=%d,j=%d,k=%d\n",i,l,j,k);

			if (!(i <= j && j < k-1 && k <= l)){
				std::cerr << "border cases: This should not have happened!, P_POmloop1" << std::endl;
				exit(EXIT_FAILURE);
			}

			if (i<=0 || j<=0 || k<=0 || l<=0 || i>n || j>n || k> n || l> n){
				std::cerr << "impossible cases: This should not have happened!, P_POmloop1" << std::endl;
				exit(EXIT_FAILURE);
			}

			energy_t min_energy = INF,temp=INF;
			cand_pos_t best_d=-1;
			for(int d = k+1; d < l; d++){
				temp = get_PO(i,j,k,d) + get_WBP(d+1,l) + beta2P(d,i);
				if (temp < min_energy){
					min_energy = temp;
					best_d=d;
				}
			}
			insert_node(i,best_d,j,k,P_PO);
			insert_node(best_d+1,l,P_WBP);
		}
			break;
	}
}


void pseudo_loop::insert_node(int i, int j, char type)
{
	seq_interval *tmp;
    tmp = new seq_interval;
    tmp->i = i;
    tmp->j = j;
    tmp->type = type;
    tmp->next = stack_interval;
    stack_interval = tmp;

}

// corresponds to region [i,k]U[l,j]
void pseudo_loop::insert_node(int i, int j, int k, int l, char type){
	seq_interval *tmp;
    tmp = new seq_interval;
    tmp->i = i;
    tmp->j = j;
	tmp->k = k;
	tmp->l = l;
	tmp->asym = -1;
    tmp->type = type;
    tmp->next = stack_interval;
    stack_interval = tmp;
}

// corresponds to region [i,k]U[l,j] with asymmetry s (used for Pxiloop5, where x=L,R,M,O)
void pseudo_loop::insert_node(int i, int j, int k, int l, int s, char type){
	seq_interval *tmp;
    tmp = new seq_interval;
    tmp->i = i;
    tmp->j = j;
	tmp->k = k;
	tmp->l = l;
	tmp->asym = s;
    tmp->type = type;
    tmp->next = stack_interval;
    stack_interval = tmp;
}

void pseudo_loop::set_stack_interval(seq_interval *stack_interval){
	this->stack_interval = stack_interval;
}
