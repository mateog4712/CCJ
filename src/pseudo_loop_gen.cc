#include "pseudo_loop.hh"
#include "h_externs.hh"
#include <stdio.h>
#include <string>
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <algorithm>
#include <cassert>


// void pseudo_loop::compute_PL(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
// 	energy_t min_energy = INF,b1=INF,b2=INF,b3=INF;

// 	const int ptype_closing = pair[S_[i]][S_[j]];

// 	if (ptype_closing>0){
// 		Index4D x(i,j,k,l);
// 		b1 = get_PLiloop(i,j,k,l);
// 		b2 = calc_PXmloop(x,MType::L);

// 		// Hosna, July 11, 2014
// 		// To avoid addition of close base pairs we check for the following here
// 		if (j>(i+TURN)){
// 			b3 = PfromL.get(i+1,j-1,k,l) + gamma2(j,i);
// 		}
// 	}
// 	min_energy = std::min({b1,b2,b3});
// 	if (min_energy < INF/2){
// 		PL.set(i,j,k,l,min_energy);
// 	}
// }

// void pseudo_loop::compute_PR(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
// 	energy_t min_energy = INF,b1=INF,b2=INF,b3=INF;

// 	const int ptype_closing = pair[S_[k]][S_[l]];

// 	if (ptype_closing>0){
// 		Index4D x(i,j,k,l);
// 		b1 = get_PRiloop(i,j,k,l);

// 		b2 = calc_PXmloop(x,MType::M);

// 		// Hosna, July 11, 2014
// 		// To avoid addition of close base pairs we check for the following here
// 		if (l>=(k+TURN+1)){
// 			b3 = PfromR.get(i,j,k+1,l-1) + gamma2(l,k);
// 		}
// 	}
// 	min_energy = std::min({b1,b2,b3});
// 	if (min_energy < INF/2){
// 		PR.set(i,j,k,l,min_energy);
// 	}
// }

// void pseudo_loop::compute_PM(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
// 	energy_t min_energy = INF,b1=INF,b2=INF,b3=INF,b4=INF;

// 	const int ptype_closing = pair[S_[j]][S_[k]];
// 	if (ptype_closing>0){
// 		Index4D x(i,j,k,l);
// 		b1 = get_PMiloop(i,j,k,l);

// 		b2 = calc_PXmloop(x,MType::M);

// 		// Hosna, July 11, 2014
// 		// To avoid addition of close base pairs we check for the following here
// 		if (k>=(j+TURN-1)){
// 			b3 = PfromM.get(i,j-1,k+1,l) + gamma2(j,k);
// 		}
// 	}
// 	min_energy = std::min({b1,b2,b3,b4});
// 	if (min_energy < INF/2){
// 		PM.set(i,j,k,l,min_energy);
// 	}
// }

// void pseudo_loop::compute_POm(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
// 	energy_t min_energy = INF,b1=INF,b2=INF,b3=INF;

// 	const int ptype_closing = pair[S_[i]][S_[l]];

// 	if (ptype_closing>0){
// 		Index4D x(i,j,k,l);
// 		b1 = get_POmiloop(i,j,k,l);

// 		b2 = calc_PXmloop(x,MType::Om);

// 		// Hosna, July 11, 2014
// 		// To avoid addition of close base pairs we check for the following here
// 		if (l>=(i+TURN+1)){
// 			b3 = PfromO.get(i+1,j,k,l-1) + gamma2(l,i);
// 		}
// 	}
// 	min_energy = std::min({b1,b2,b3});
// 	if (min_energy < INF/2){
// 		POm.set(i,j,k,l,min_energy);
// 	}
// }

// void pseudo_loop::compute_POs(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
// 	energy_t min_energy = INF,b1=INF,b2=INF,b3=INF;

// 	const int ptype_closing = pair[S_[i]][S_[l]];

// 	if (ptype_closing>0){
// 		Index4D x(i,j,k,l);
// 		b1 = get_POsiloop(i,j,k,l);

// 		b2 = calc_PXmloop(x,MType::Os);

// 		if(i==j && k==l){
// 			b3=gamma2(l,i);
// 		}
// 	}
// 	min_energy = std::min({b1,b2,b3});
// 	if (min_energy < INF/2){
// 		POs.set(i,j,k,l,min_energy);
// 	}
// }

// void pseudo_loop::compute_PLreO(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
// 	energy_t min_energy = INF,b1=INF,b2=INF,b3=INF;

// 	const int ptype_closing = pair[S_[i]][S_[j]];

// 	if (ptype_closing>0){
// 		Index4D x(i,j,k,l);
// 		b1 = get_PLreOiloop(i,j,k,l);
// 		b2 = calc_PXmloop(x,MType::LreO);

// 		// Hosna, July 11, 2014
// 		// To avoid addition of close base pairs we check for the following here
// 		if (j>(i+TURN)){
// 			b3 = PfromLreO.get(i+1,j-1,k,l) + gamma2(j,i);
// 		}
// 	}
// 	min_energy = std::min({b1,b2,b3});
// 	if (min_energy < INF/2){
// 		PLreO.set(i,j,k,l,min_energy);
// 	}
// }

// void pseudo_loop::compute_PLreR(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
// 	energy_t min_energy = INF,b1=INF,b2=INF,b3=INF;

// 	const int ptype_closing = pair[S_[i]][S_[j]];

// 	if (ptype_closing>0){
// 		Index4D x(i,j,k,l);
// 		b1 = get_PLreRiloop(i,j,k,l);
// 		b2 = calc_PXmloop(x,MType::LreR);

// 		// Hosna, July 11, 2014
// 		// To avoid addition of close base pairs we check for the following here
// 		if (j>(i+TURN)){
// 			b3 = PfromLreR.get(i+1,j-1,k,l) + gamma2(j,i);
// 		}
// 	}
// 	min_energy = std::min({b1,b2,b3});
// 	if (min_energy < INF/2){
// 		PLreR.set(i,j,k,l,min_energy);
// 	}
// }

// void pseudo_loop::compute_PMreO(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
// 	energy_t min_energy = INF,b1=INF,b2=INF,b3=INF,b4=INF;

// 	const int ptype_closing = pair[S_[j]][S_[k]];
// 	if (ptype_closing>0){
// 		Index4D x(i,j,k,l);
// 		b1 = get_PMreOiloop(i,j,k,l);
// 		b2 = calc_PXmloop(x,MType::MreO);

// 		// Hosna, July 11, 2014
// 		// To avoid addition of close base pairs we check for the following here
// 		if (k>=(j+TURN-1)){
// 			b3 = PfromMreO.get(i,j-1,k+1,l) + gamma2(j,k);
// 		}
// 	}
// 	min_energy = std::min({b1,b2,b3,b4});
// 	if (min_energy < INF/2){
// 		PMreO.set(i,j,k,l,min_energy);
// 	}
// }

// void pseudo_loop::compute_PMreR(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
// 	energy_t min_energy = INF,b1=INF,b2=INF,b3=INF,b4=INF;

// 	const int ptype_closing = pair[S_[j]][S_[k]];
// 	if (ptype_closing>0){
// 		Index4D x(i,j,k,l);
// 		b1 = get_PMreRiloop(i,j,k,l);
// 		b2 = calc_PXmloop(x,MType::MreR);

// 		// Hosna, July 11, 2014
// 		// To avoid addition of close base pairs we check for the following here
// 		if (k>=(j+TURN-1)){
// 			b3 = PfromMreR.get(i,j-1,k+1,l) + gamma2(j,k);
// 		}
// 	}
// 	min_energy = std::min({b1,b2,b3,b4});
// 	if (min_energy < INF/2){
// 		PMreR.set(i,j,k,l,min_energy);
// 	}
// }

// void pseudo_loop::compute_PLMreR(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
// 	energy_t min_energy = INF,b1=INF,b2=INF,b3=INF;

// 	const int ptype_closing = pair[S_[i]][S_[j]];

// 	if (ptype_closing>0){
// 		Index4D x(i,j,k,l);
// 		b1 = get_PLMreRiloop(i,j,k,l);
// 		b2 = calc_PXmloop(x,MType::LMreR);

// 		// Hosna, July 11, 2014
// 		// To avoid addition of close base pairs we check for the following here
// 		if (j>(i+TURN)){
// 			b3 = PfromLMreR.get(i+1,j-1,k,l) + gamma2(j,i);
// 		}
// 	}
// 	min_energy = std::min({b1,b2,b3});
// 	if (min_energy < INF/2){
// 		PLMreR.set(i,j,k,l,min_energy);
// 	}
// }

// void pseudo_loop::compute_PLMorO(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
// 	energy_t min_energy = INF,b1=INF,b2=INF,b3=INF;

// 	const int ptype_closing = pair[S_[i]][S_[j]];

// 	if (ptype_closing>0){
// 		Index4D x(i,j,k,l);
// 		b1 = get_PLMorOiloop(i,j,k,l);
// 		b2 = calc_PXmloop(x,MType::LMorO);

// 		// Hosna, July 11, 2014
// 		// To avoid addition of close base pairs we check for the following here
// 		if (j>(i+TURN)){
// 			b3 = PfromLMorO.get(i+1,j-1,k,l) + gamma2(j,i);
// 		}
// 	}
// 	min_energy = std::min({b1,b2,b3});
// 	if (min_energy < INF/2){
// 		PLMorO.set(i,j,k,l,min_energy);
// 	}
// }

// I will decide at a later time if this should be grabbing from a matrix or calculated every time it is called
energy_t pseudo_loop::calc_PLiloop(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	if (!(i <= j && j < k-1 && k <= l)){
		return INF;
	}
	assert(!(i<=0 || l> n));

	if(!can_pair(i,j)) return INF;

	energy_t min_energy = INF;
	if (i+TURN+2<j) { 
		min_energy = PL.get(i+1,j-1,k,l) + get_e_stP(i,j);
	}
	cand_pos_t max_d = std::min(j,i+MAXLOOP);
	for(cand_pos_t d= i+1; d<max_d; ++d){
		cand_pos_t min_dp = std::max(d+TURN,j-MAXLOOP);
		for(cand_pos_t dp = j-1; dp > min_dp; --dp){
			if (!can_pair(d,dp)) continue;
			min_energy = std::min(min_energy,get_e_intP(i,d,dp,j) + PL.get(d,dp,k,l));
		}
	}
	return min_energy;
}

energy_t pseudo_loop::calc_PRiloop(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	if (!(i <= j && j < k-1 && k <= l)){
		return INF;
	}
	assert(!(i<=0 || l> n));

	if(!can_pair(k,l)) return INF;

	energy_t min_energy = INF;
	if (k+TURN+2<l) { 
		min_energy = PR.get(i,j,k+1,l-1) + get_e_stP(k,l);
	}
	cand_pos_t max_d = std::min(l,k+MAXLOOP);
	for(cand_pos_t d= k+1; d<max_d; ++d){
		cand_pos_t min_dp = std::max(d+TURN,l-MAXLOOP);
		for(cand_pos_t dp=l-1; dp > min_dp; --dp){
			if (!can_pair(d,dp)) continue;
			min_energy = std::min(min_energy,get_e_intP(k,d,dp,l) + PR.get(i,j,d,dp));
		}
	}
	return min_energy;
}

energy_t pseudo_loop::calc_PMiloop(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	if (!(i <= j && j < k-1 && k <= l)){
		return INF;
	}
	assert(!(i<=0 || l> n));

	if(!can_pair(j,k)) return INF;

	energy_t min_energy = INF;
	if ( i<j && k<l ) {
		min_energy = PM.get(i,j-1,k+1,l) + get_e_stP(j-1,k+1);
	}
	cand_pos_t max_d = std::max(i,j-MAXLOOP);
	for(cand_pos_t d= j-1; d>max_d; --d){
		cand_pos_t min_dp = std::min(l,k+MAXLOOP);
		for (cand_pos_t dp=k+1; dp <min_dp; ++dp) {
			if (!can_pair(d,dp)) continue;
			min_energy = std::min(min_energy,get_e_intP(d,j,k,dp) + PM.get(i,d,dp,l));
		}
	}
	return min_energy;
}

energy_t pseudo_loop::calc_POmiloop(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	if (!(i <= j && j < k-1 && k <= l)){
		return INF;
	}
	assert(!(i<=0 || l> n));

	if(!can_pair(i,l)) return INF;

	energy_t min_energy = INF;
	if (i<j && k<l ) { 
		min_energy = POm.get(i+1,j,k,l-1) + get_e_stP(i,l);
	}
	cand_pos_t max_d = std::min(j,i+MAXLOOP);
	for(cand_pos_t d= i+1; d<max_d; ++d){
		cand_pos_t min_dp = std::max(l-MAXLOOP,k);
		for (cand_pos_t dp=l-1; dp >min_dp; --dp) {
			if (!can_pair(d,dp)) continue;
			min_energy = std::min(min_energy,get_e_intP(i,d,dp,l) + POm.get(d,j,dp,k));
		}
	}
	return min_energy;
}

energy_t pseudo_loop::calc_POsiloop(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	if (!(i <= j && j < k-1 && k <= l)){
		return INF;
	}
	assert(!(i<=0 || l> n));

	if(!can_pair(i,l)) return INF;

	energy_t min_energy = INF;
	if (i<j && k<l ) { 
		min_energy = POs.get(i+1,j,k,l-1) + get_e_stP(i,l);
	}
	cand_pos_t max_d = std::min(j,i+MAXLOOP);
	for(cand_pos_t d= i+1; d<max_d; ++d){
		cand_pos_t min_dp = std::max(l-MAXLOOP,k);
		for (cand_pos_t dp=l-1; dp >min_dp; --dp) {
			if (!can_pair(d,dp)) continue;
			min_energy = std::min(min_energy,get_e_intP(i,d,dp,l) + POs.get(d,j,dp,k));
		}
	}
	return min_energy;
}

energy_t pseudo_loop::calc_PLreOiloop(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	if (!(i <= j && j < k-1 && k <= l)){
		return INF;
	}
	assert(!(i<=0 || l> n));

	if(!can_pair(i,j)) return INF;

	energy_t min_energy = INF;
	if (i+TURN+2<j) { 
		min_energy = PLreO.get(i+1,j-1,k,l) + get_e_stP(i,j);
	}
	cand_pos_t max_d = std::min(j,i+MAXLOOP);
	for(cand_pos_t d= i+1; d<max_d; ++d){
		cand_pos_t min_dp = std::max(d+TURN,j-MAXLOOP);
		for(cand_pos_t dp = j-1; dp > min_dp; --dp){
			if (!can_pair(d,dp)) continue;
			min_energy = std::min(min_energy,get_e_intP(i,d,dp,j) + PLreO.get(d,dp,k,l));
		}
	}
	return min_energy;
}

energy_t pseudo_loop::calc_PLreRiloop(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	if (!(i <= j && j < k-1 && k <= l)){
		return INF;
	}
	assert(!(i<=0 || l> n));

	if(!can_pair(i,j)) return INF;

	energy_t min_energy = INF;
	if (i+TURN+2<j) { 
		min_energy = PLreR.get(i+1,j-1,k,l) + get_e_stP(i,j);
	}
	cand_pos_t max_d = std::min(j,i+MAXLOOP);
	for(cand_pos_t d= i+1; d<max_d; ++d){
		cand_pos_t min_dp = std::max(d+TURN,j-MAXLOOP);
		for(cand_pos_t dp = j-1; dp > min_dp; --dp){
			if (!can_pair(d,dp)) continue;
			min_energy = std::min(min_energy,get_e_intP(i,d,dp,j) + PLreR.get(d,dp,k,l));
		}
	}
	return min_energy;
}

energy_t pseudo_loop::calc_PMreOiloop(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	if (!(i <= j && j < k-1 && k <= l)){
		return INF;
	}
	assert(!(i<=0 || l> n));

	if(!can_pair(j,k)) return INF;

	energy_t min_energy = INF;
	if ( i<j && k<l ) {
		min_energy = PMreO.get(i,j-1,k+1,l) + get_e_stP(j-1,k+1);
	}
	cand_pos_t max_d = std::max(i,j-MAXLOOP);
	for(cand_pos_t d= j-1; d>max_d; --d){
		cand_pos_t min_dp = std::min(l,k+MAXLOOP);
		for (cand_pos_t dp=k+1; dp <min_dp; ++dp) {
			if (!can_pair(d,dp)) continue;
			min_energy = std::min(min_energy,get_e_intP(d,j,k,dp) + PMreO.get(i,d,dp,l));
		}
	}
	return min_energy;
}

energy_t pseudo_loop::calc_PMreRiloop(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	if (!(i <= j && j < k-1 && k <= l)){
		return INF;
	}
	assert(!(i<=0 || l> n));

	if(!can_pair(j,k)) return INF;

	energy_t min_energy = INF;
	if ( i<j && k<l ) {
		min_energy = PMreR.get(i,j-1,k+1,l) + get_e_stP(j-1,k+1);
	}
	cand_pos_t max_d = std::max(i,j-MAXLOOP);
	for(cand_pos_t d= j-1; d>max_d; --d){
		cand_pos_t min_dp = std::min(l,k+MAXLOOP);
		for (cand_pos_t dp=k+1; dp <min_dp; ++dp) {
			if (!can_pair(d,dp)) continue;
			min_energy = std::min(min_energy,get_e_intP(d,j,k,dp) + PMreR.get(i,d,dp,l));
		}
	}
	return min_energy;
}

energy_t pseudo_loop::calc_PLMreRiloop(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	if (!(i <= j && j < k-1 && k <= l)){
		return INF;
	}
	assert(!(i<=0 || l> n));

	if(!can_pair(i,j)) return INF;

	energy_t min_energy = INF;
	if (i+TURN+2<j) { 
		min_energy = PLMreR.get(i+1,j-1,k,l) + get_e_stP(i,j);
	}
	cand_pos_t max_d = std::min(j,i+MAXLOOP);
	for(cand_pos_t d= i+1; d<max_d; ++d){
		cand_pos_t min_dp = std::max(d+TURN,j-MAXLOOP);
		for(cand_pos_t dp = j-1; dp > min_dp; --dp){
			if (!can_pair(d,dp)) continue;
			min_energy = std::min(min_energy,get_e_intP(i,d,dp,j) + PLMreR.get(d,dp,k,l));
		}
	}
	return min_energy;
}

energy_t pseudo_loop::calc_PLMorOiloop(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	if (!(i <= j && j < k-1 && k <= l)){
		return INF;
	}
	assert(!(i<=0 || l> n));

	if(!can_pair(i,j)) return INF;

	energy_t min_energy = INF;
	if (i+TURN+2<j) { 
		min_energy = PLMorO.get(i+1,j-1,k,l) + get_e_stP(i,j);
	}
	cand_pos_t max_d = std::min(j,i+MAXLOOP);
	for(cand_pos_t d= i+1; d<max_d; ++d){
		cand_pos_t min_dp = std::max(d+TURN,j-MAXLOOP);
		for(cand_pos_t dp = j-1; dp > min_dp; --dp){
			if (!can_pair(d,dp)) continue;
			min_energy = std::min(min_energy,get_e_intP(i,d,dp,j) + PLMorO.get(d,dp,k,l));
		}
	}
	return min_energy;
}