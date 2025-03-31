#include "pseudo_loop.hh"
#include "h_externs.hh"
#include "common.hh"
#include <stdio.h>
#include <string>
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <algorithm>

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

	for(cand_pos_t i = 1; i <= total_length; ++i) {
		std::vector<energy_t> row;
		row.resize(total_length,INF);
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

		std::vector< std::vector< energy_t> > row_row;
		for (cand_pos_t j=1; j<=MAXLOOP; ++j){
			std::vector<energy_t> row_5;
			row_5.resize(MAXLOOP,INF);
			row_row.push_back(row_5);
		}

		PLiloop5.push_back(row_row);
		PRiloop5.push_back(row_row);
		PMiloop5.push_back(row_row);
		POiloop5.push_back(row_row);

	}
	WMB.resize(total_length,INF);
}

pseudo_loop::~pseudo_loop()
{
}

void pseudo_loop::compute_energies(int i, int l)
{

	// 1) compute all energies over region [i,l]
	compute_P(i,l);

	compute_WBP(i,l);

	compute_WPP(i,l);

	//2) compute all energies over gapped region [i,j]U[k,l]
	for(int j = i; j<l; j++){
		// Hosna, July 8, 2014
		// in original recurrences we have j< k-1, so I am changing k=j+1 to k=j+2
		for(int k = l; k>=j+2; k--){//for(int k = j+2; k<=l; k++){

			// 3) compute all energies PXiloop5 for all s
			for(int s=0; s<MAXLOOP; s++){

				compute_PLiloop5(i,j,k,l,s);

				compute_PRiloop5(i,j,k,l,s);

				compute_PMiloop5(i,j,k,l,s);

				compute_POiloop5(i,j,k,l,s);
			}

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
	for(cand_pos_t d=i; d< l; d++){
		for(cand_pos_t e = d+1; e<= l; e++){
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
	for(cand_pos_t d=i; d< l; d++){
		for(cand_pos_t e = d+1; e<= l; e++){
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
	for(cand_pos_t j=i; j< l; j++){
		for (cand_pos_t d=j+1; d<l; d++){
			for (cand_pos_t k=d+1; k<l; k++){
				b1 = get_PK(i,j,d+1,k) +get_PK(j+1,d,k+1,l);
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
	for(cand_pos_t d=i+1; d< j; d++){
		energy_t tmp = get_PK(i,d,k,l) + get_WP(d+1,j);
		b1=std::min(b1,tmp);
	}

	// Hosna, july 8, 2014
	// based on original recurrences we should have d<l, and
	// it is not clear to me why we have d<=l here, so I am changing this back to original
	// by changing d<=l to d<l
	for(int d=k+1; d< l; d++){
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
	//	printf (">>>>>>>>> PR(%d,%d,%d,%d) branch %d energy %d\n", i, j, k,l,best_branch, min_energy);

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

	for(cand_pos_t d=i+1; d< j; d++){
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
	for(cand_pos_t d=k+1; d< l; d++){
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

	for(cand_pos_t d=i+1; d< j; d++){
		energy_t tmp = get_PfromM(i,d,k,l)+get_WP(d+1,j);
		b1 = std::min(b1,tmp);
	}
	for(cand_pos_t d=k+1; d<l; d++){
		energy_t tmp=get_PfromM(i,j,d,l)+get_WP(k,d-1);
		b2 = std::min(b2,tmp);
	}
	//Hosna, July 28, 2014
	// I think going from PfromX to PX we are changing bands and so should be paying a band penalty
	b3 = get_PL(i,j,k,l) + gamma2(j,i)+ PB_penalty;//;
	//Hosna, July 28, 2014
	// I think going from PfromX to PX we are changing bands and so should be paying a band penalty
	b4 = get_PR(i,j,k,l) + gamma2(l,k)+PB_penalty;

	//Hosna, May 2, 2014
	// I think I should block going from PM to PO as it will solve the same band from 2 sides jumping over possible internal loops
	/*
	b5 = get_PO(i,j,k,l) + gamma2(l,i);
	if(b5 < min_energy){
		min_energy = b5;
	}
	*/
	min_energy = std::min({min_energy,b1,b2,b3,b4});
	if (min_energy < INF/2){
		PfromM[ij][kl]=min_energy;
	}
}

void pseudo_loop::compute_PfromO(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	energy_t min_energy = INF,b1=INF,b2=INF,b3=INF,b4=INF;

	cand_pos_t ij = index[i]+j-i;
	cand_pos_t kl = index[k]+l-k;

	for(cand_pos_t d=i+1; d< j; d++){
		energy_t tmp=get_PfromO(d,j,k,l)+get_WP(i,d-1);
		b1 = std::min(b1,tmp);
	}
	for(cand_pos_t d=k+1; d<l; d++){
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
	energy_t min_energy = INF,b1=INF,b2=INF;

	cand_pos_t ij = index[i]+j-i;
	cand_pos_t kl = index[k]+l-k;
	int ptype_closing = pair[S_[i]][S_[j]];
	if (ptype_closing>0){

		if (i+1 < n && j-1 >= 0){
			b1 = get_PL(i+1,j-1,k,l) + get_e_stP(i,j);
			// Hosna, August 21, 2014
			// revising the max_s value
			// there are j-i+1 bases between i and j, from which we need an ip and jp and at least 3 bases between ip and jp=> j-i+1-2-3
			//int max_s = MIN(MAX(j-i-5,0),MAXLOOP-1);
			cand_pos_t max_s = std::min(std::max(j-i-4,0),MAXLOOP-1);
			// Hosna April 11, 2014 changed s=0 to s=1 to avoid considering stack as internal loop
			for(int s = 1; s <= max_s; s++){
				// Hosna, April 2, 2014
				// since in h_common.cpp I don't have access to int_sequence, I am changing alpha2P definition to accept 4 values
				//int temp = get_PLiloop5(i,j,k,l,s) + alpha0P + alpha2P(j,i);
				// Hosna, August 21, 2014 changing indexes in alpha2P so the first index is always less than the second.
				// so the order of alpha2P(i,l,i+1,l-1)
				//int temp = get_PLiloop5(i,j,k,l,s) + alpha0P + alpha2P(int_sequence[j],int_sequence[i],int_sequence[j-1],int_sequence[i+1]);
				// energy_t tmp = get_PLiloop5(i,j,k,l,s) + alpha0P + alpha2P(int_sequence[i],int_sequence[j],int_sequence[i+1],int_sequence[j-1]);
				energy_t tmp = get_PLiloop5(i,j,k,l,s) + get_e_stP(i,j);
				b2 = std::min(b2,tmp);
			}
		}
		min_energy = std::min(b1,b2);
	}
	if (min_energy < INF/2){
		PLiloop[ij][kl]=min_energy;
	}
}

void pseudo_loop::compute_PLiloop5(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l, cand_pos_t s){
	energy_t min_energy = INF,b1=INF,b2=INF,b3=INF;

	cand_pos_t ij = index[i]+j-i;
	cand_pos_t kl = index[k]+l-k;

	if (s >= 2 && (i+1) < n && (j-1)>=0){
		b1 = get_PLiloop5(i+1,j-1,k,l,s-2) + alpha1P(2); // Fix here Mateo

	}
	// Hosna, April 2, 2014
	// since in h_common.cpp I don't have access to int_sequence, I am changing alpha2P definition to accept 4 values
	//b2 = get_PL(i+s+1,j-1,k,l) + alpha1P(s) + alpha3P(3) + alpha2P(j-1,i+s+1);
	if((j-2) >= 0 && (i+s+2) < n){
		// Hosna, August 21, 2014 changing indexes in alpha2P so the first index is always less than the second.
		// so the order of alpha2P(i,l,i+1,l-1)
		//b2 = get_PL(i+s+1,j-1,k,l) + alpha1P(s) + alpha3P(s) + alpha2P(int_sequence[j-1],int_sequence[i+s+1],int_sequence[j-2],int_sequence[i+s+2]);
		// b2 = get_PL(i+s+1,j-1,k,l) + alpha1P(s) + alpha3P(s) + alpha2P(int_sequence[i+s+1],int_sequence[j-1],int_sequence[i+s+2],int_sequence[j-2]);
		b2 = get_PL(i+s+1,j-1,k,l) + get_e_intP(i+s+1,i+s+2,j-2,j-1);
		//Hosna, August 21, 2014, adding V here is wrong! since a PX will handle its energy in other recurrences
		/*
		// Hosna, April 11, 2014
		// This branch needs to add the energy of the hairpin loop to the internal loop as well, so adding V->get_energy()*penalty
		// because we have internal loop here, we use the internal loop penalty
		int hairpin_energy = (int)round(e_intP_penalty * (double)V->get_energy(i+s+1,j-1));
		b2 += hairpin_energy;
		 */
	}


	// Hosna, April 2, 2014
	// since in h_common.cpp I don't have access to int_sequence, I am changing alpha2P definition to accept 4 values
	//b3 = get_PL(i+1,j-s-1,k,l) + alpha1P(s) + alpha3P(3) + alpha2P(j-s-1,i+1);
	if ((j-s-2) >= 0 && (i+2) < n){
		// Hosna, August 21, 2014 changing indexes in alpha2P so the first index is always less than the second.
		// so the order of alpha2P(i,l,i+1,l-1)
		//b3 = get_PL(i+1,j-s-1,k,l) + alpha1P(s) + alpha3P(s) + alpha2P(int_sequence[j-s-1],int_sequence[i+1],int_sequence[j-s-2],int_sequence[i+2]);
		// b3 = get_PL(i+1,j-s-1,k,l) + alpha1P(s) + alpha3P(s) + alpha2P(int_sequence[i+1],int_sequence[j-s-1],int_sequence[i+2],int_sequence[j-s-1]);
		b3 = get_PL(i+1,j-s-1,k,l) + get_e_intP(i+1,i+2,j-s-1,j-s-1);
		//Hosna, August 21, 2014, adding V here is wrong! since a PX will handle its energy in other recurrences
		/*
		// Hosna, April 11, 2014
		// This branch needs to add the energy of the hairpin loop to the internal loop as well, so adding V->get_energy()*penalty
		// because we have internal loop here, we use the internal loop penalty
		int hairpin_energy = (int)round(e_intP_penalty * (double)V->get_energy(i+1,j-s-1));
		b3 += hairpin_energy;
		 */

	}

	min_energy = std::min({min_energy,b1,b2,b3});
	if (min_energy < INF/2){
		PLiloop5[ij][kl][s] = min_energy;
	}

}

void pseudo_loop::compute_PLmloop(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	energy_t min_energy = INF,b1=INF,b2=INF;

	cand_pos_t ij = index[i]+j-i;
	cand_pos_t kl = index[k]+l-k;

	for(cand_pos_t d = i+1; d < j-1; d++){
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
	int min_energy = INF;

	cand_pos_t ij = index[i]+j-i;
	cand_pos_t kl = index[k]+l-k;
	for(cand_pos_t d = i+1; d < j; d++){
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
	for(cand_pos_t d = i+1; d < j; d++){
		// Hosna, feb 23, 2014
		// changed the recurrences so that j-1 is accounted for in PLmloop
		energy_t tmp = get_PL(i,d,k,l) + get_WBP(d+1,j) + beta2P(d,i);
		min_energy = std::min(min_energy,tmp);
	}
	if (min_energy < INF/2){
		PLmloop1[ij][kl] = min_energy;
	}

}


energy_t pseudo_loop::get_WB(cand_pos_t i, cand_pos_t j){
	return (std::min(beta1P*(j-i+1),get_WBP(i,j)));
}

energy_t pseudo_loop::get_WBP(cand_pos_t i, cand_pos_t j){
	cand_pos_t ij = index[i]+j-i;
	return WBP[ij];
}

energy_t pseudo_loop::get_WP(cand_pos_t i, cand_pos_t j){
	return (std::min(gamma1*(j-i+1),get_WPP(i,j)));
}

energy_t pseudo_loop::get_WPP(cand_pos_t i, cand_pos_t j){
	cand_pos_t ij = index[i]+j-i;
	return WPP[ij];
}

energy_t pseudo_loop::get_P(cand_pos_t i, cand_pos_t j){
	cand_pos_t ij = index[i]+j-i;
	return P[ij];
}

energy_t pseudo_loop::get_PK(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){

	cand_pos_t ij = index[i]+j-i;
	cand_pos_t kl = index[k]+l-k;

	return PK[ij][kl];
}

energy_t pseudo_loop::get_PL(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){

	cand_pos_t ij = index[i]+j-i;
	cand_pos_t kl = index[k]+l-k;

	return PL[ij][kl];
}

energy_t pseudo_loop::get_PR(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){

	cand_pos_t ij = index[i]+j-i;
	cand_pos_t kl = index[k]+l-k;

	return PR[ij][kl];
}

energy_t pseudo_loop::get_PM(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){

	cand_pos_t ij = index[i]+j-i;
	cand_pos_t kl = index[k]+l-k;
	if (i ==j && k ==l){
		return gamma2(i,l);
	}

	return PM[ij][kl];
}

energy_t pseudo_loop::get_PO(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){

	cand_pos_t ij = index[i]+j-i;
	cand_pos_t kl = index[k]+l-k;

	return PO[ij][kl];
}

energy_t pseudo_loop::get_PfromL(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){

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
	const int ptype_closing = pair[S_[i]][S_[j]];
	if(ptype_closing == 0) return INF;

	cand_pos_t ij = index[i]+j-i;
	cand_pos_t kl = index[k]+l-k;

	return PLiloop[ij][kl];
}

energy_t pseudo_loop::get_PLiloop5(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l, cand_pos_t s){

	cand_pos_t ij = index[i]+j-i;
	cand_pos_t kl = index[k]+l-k;

	return PLiloop5[ij][kl][s];
}

energy_t pseudo_loop::get_PLmloop(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){

	cand_pos_t ij = index[i]+j-i;
	cand_pos_t kl = index[k]+l-k;

	return PLmloop[ij][kl];
}

energy_t pseudo_loop::get_PLmloop0(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){

	cand_pos_t ij = index[i]+j-i;
	cand_pos_t kl = index[k]+l-k;

	return PLmloop0[ij][kl];
}

energy_t pseudo_loop::get_PLmloop1(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){

	cand_pos_t ij = index[i]+j-i;
	cand_pos_t kl = index[k]+l-k;

	return PLmloop1[ij][kl];
}

energy_t pseudo_loop::get_PRiloop(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	const int ptype_closing = pair[S_[i]][S_[j]];
	if(ptype_closing == 0) return INF;

	cand_pos_t ij = index[i]+j-i;
	cand_pos_t kl = index[k]+l-k;

	return PRiloop[ij][kl];
}

energy_t pseudo_loop::get_PRiloop5(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l, cand_pos_t s){

	cand_pos_t ij = index[i]+j-i;
	cand_pos_t kl = index[k]+l-k;

	return PRiloop5[ij][kl][s];
}

energy_t pseudo_loop::get_PRmloop(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){

	cand_pos_t ij = index[i]+j-i;
	cand_pos_t kl = index[k]+l-k;

	return PRmloop[ij][kl];
}

energy_t pseudo_loop::get_PRmloop0(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){

	cand_pos_t ij = index[i]+j-i;
	cand_pos_t kl = index[k]+l-k;

	return PRmloop0[ij][kl];
}

energy_t pseudo_loop::get_PRmloop1(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){

	cand_pos_t ij = index[i]+j-i;
	cand_pos_t kl = index[k]+l-k;

	return PRmloop1[ij][kl];
}

energy_t pseudo_loop::get_PMiloop(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	const int ptype_closing = pair[S_[i]][S_[j]];
	if(ptype_closing == 0) return INF;

	cand_pos_t ij = index[i]+j-i;
	cand_pos_t kl = index[k]+l-k;

	return PMiloop[ij][kl];
}

energy_t pseudo_loop::get_PMiloop5(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l, cand_pos_t s){

	cand_pos_t ij = index[i]+j-i;
	cand_pos_t kl = index[k]+l-k;

	return PMiloop5[ij][kl][s];
}

energy_t pseudo_loop::get_PMmloop(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){

	cand_pos_t ij = index[i]+j-i;
	cand_pos_t kl = index[k]+l-k;

	return PMmloop[ij][kl];
}

energy_t pseudo_loop::get_PMmloop0(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){

	cand_pos_t ij = index[i]+j-i;
	cand_pos_t kl = index[k]+l-k;

	return PMmloop0[ij][kl];
}

energy_t pseudo_loop::get_PMmloop1(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){

	cand_pos_t ij = index[i]+j-i;
	cand_pos_t kl = index[k]+l-k;

	return PMmloop1[ij][kl];
}

energy_t pseudo_loop::get_POiloop(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	const int ptype_closing = pair[S_[i]][S_[j]];
	if(ptype_closing == 0) return INF;

	cand_pos_t ij = index[i]+j-i;
	cand_pos_t kl = index[k]+l-k;

	return POiloop[ij][kl];
}

energy_t pseudo_loop::get_POiloop5(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l, cand_pos_t s){

	cand_pos_t ij = index[i]+j-i;
	cand_pos_t kl = index[k]+l-k;

	return POiloop5[ij][kl][s];
}

energy_t pseudo_loop::get_POmloop(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){

	cand_pos_t ij = index[i]+j-i;
	cand_pos_t kl = index[k]+l-k;

	return POmloop[ij][kl];
}

energy_t pseudo_loop::get_POmloop0(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){

	cand_pos_t ij = index[i]+j-i;
	cand_pos_t kl = index[k]+l-k;

	return POmloop0[ij][kl];
}

energy_t pseudo_loop::get_POmloop1(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){

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
energy_t pseudo_loop::get_WMB(cand_pos_t i, cand_pos_t j){
	return INF+i+j;
}

void pseudo_loop::set_stack_interval(seq_interval *stack_interval){
	this->stack_interval = stack_interval;
}
