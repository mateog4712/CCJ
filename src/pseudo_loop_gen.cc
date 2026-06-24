#include "pseudo_loop.hh"
#include "h_externs.hh"
#include "W_final.hh"
#include <stdio.h>
#include <string>
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <algorithm>
#include <cassert>
#define debug 0

/**
 * These should probably just be in a matrix at this point; there are so many arrays, why have four less and recompute every time; Remember to probably switch this!!
 */
energy_t pseudo_loop::calc_PLiloop(const Index4D &x, MType type){
	if(impossible_case(x)) return INF;
	if(!can_pair(x.lend(type),x.rend(type))) return INF;
	const cand_pos_t i = x.i(), j = x.j(), k = x.k(), l = x.l();

	Matrix4D &PX = PX_by_mtype(type);
	energy_t min_energy = INF;
	if (i+TURN+2<j) { 
		min_energy = PX.get(i+1,j-1,k,l) + get_e_stP(i,j);
	}
	cand_pos_t max_d = std::min(j,i+MAXLOOP);
	for(cand_pos_t d= i+1; d<max_d; ++d){
		cand_pos_t min_dp = std::max(d+TURN,j-MAXLOOP);
		for(cand_pos_t dp = j-1; dp > min_dp; --dp){
			if (!can_pair(d,dp)) continue;
			min_energy = std::min(min_energy,get_e_intP(i,d,dp,j) + PX.get(d,dp,k,l));
		}
	}
	return min_energy;
}

// I know that stacks flip the second stack basepair, do I need to consider that here and when is k.j better done than j.k
// Would it make it easier for understanding since we are moving outward instead of inward?
energy_t pseudo_loop::calc_PMiloop(const Index4D &x, MType type){
	if(impossible_case(x)) return INF;
	if(!can_pair(x.lend(type),x.rend(type))) return INF;
	const cand_pos_t i = x.i(), j = x.j(), k = x.k(), l = x.l();

	Matrix4D &PX = PX_by_mtype(type);
	energy_t min_energy = INF;
	if (i<j && k<l) {
		min_energy = PX.get(i,j-1,k+1,l) + get_e_stP(j-1,k+1);
	}
	cand_pos_t max_d = std::max(i,j-MAXLOOP);
	for(cand_pos_t d= j-1; d>max_d; --d){
		cand_pos_t min_dp = std::min(l,k+MAXLOOP); // could switch these here so that we are increasing in the first for like all the others
		for (cand_pos_t dp=k+1; dp <min_dp; ++dp) {
			if (!can_pair(d,dp)) continue;
			min_energy = std::min(min_energy,get_e_intP(d,j,k,dp) + PX.get(i,d,dp,l));
		}
	}
	return min_energy;
}

energy_t pseudo_loop::calc_PRiloop(const Index4D &x, MType type){
	if(impossible_case(x)) return INF;
	if(!can_pair(x.lend(type),x.rend(type))) return INF;
	const cand_pos_t i = x.i(), j = x.j(), k = x.k(), l = x.l();

	Matrix4D &PX = PX_by_mtype(type);
	energy_t min_energy = INF;
	if (k+TURN+2<l) { 
		min_energy = PX.get(i,j,k+1,l-1) + get_e_stP(k,l);
	}
	cand_pos_t max_d = std::min(l,k+MAXLOOP);
	for(cand_pos_t d= k+1; d<max_d; ++d){
		cand_pos_t min_dp = std::max(d+TURN,l-MAXLOOP);
		for(cand_pos_t dp=l-1; dp > min_dp; --dp){
			if (!can_pair(d,dp)) continue;
			min_energy = std::min(min_energy,get_e_intP(k,d,dp,l) + PX.get(i,j,d,dp));
		}
	}
	return min_energy;
}

energy_t pseudo_loop::calc_POiloop(const Index4D &x, MType type){
	if(impossible_case(x)) return INF;
	if(!can_pair(x.lend(type),x.rend(type))) return INF;
	const cand_pos_t i = x.i(), j = x.j(), k = x.k(), l = x.l();

	Matrix4D &PX = PX_by_mtype(type);
	energy_t min_energy = INF;
	if (i<j && k<l ) { 
		min_energy = PX.get(i+1,j,k,l-1) + get_e_stP(i,l);
	}
	cand_pos_t max_d = std::min(j,i+MAXLOOP);
	for(cand_pos_t d= i+1; d<max_d; ++d){
		cand_pos_t min_dp = std::max(l-MAXLOOP,k);
		for (cand_pos_t dp=l-1; dp >min_dp; --dp) {
			if (!can_pair(d,dp)) continue;
			min_energy = std::min(min_energy,get_e_intP(i,d,dp,l) + PX.get(d,j,k,dp));
		}
	}
	return min_energy;
}

void pseudo_loop::compute_PfromL(const Index4D &x, MType type){
	energy_t min_energy = INF;

	Matrix4D &PfromXprime = PfromXprime_by_mtype(type);
	for(cand_pos_t d=x.i()+1; d<=x.j(); ++d){
		energy_t tmp = calc_WP(x.i(),d-1) + PfromXprime.get(d,x.j(),x.k(),x.l());
		min_energy = std::min(min_energy,tmp);
	}

	if (min_energy < INF/2){
		Matrix4D &PfromX = PfromX_by_mtype(type);
		PfromX.set(x.i(),x.j(),x.k(),x.l(),min_energy);
	}
}

void pseudo_loop::compute_PfromM(const Index4D &x, MType type){
	energy_t min_energy = INF;

	Matrix4D &PfromXprime = PfromXprime_by_mtype(type);
	for(cand_pos_t d=x.i(); d<x.j(); ++d){
		energy_t tmp = PfromXprime.get(x.i(),d,x.k(),x.l()) + calc_WP(d+1,x.j());
		min_energy = std::min(min_energy,tmp);
	}

	if (min_energy < INF/2){
		Matrix4D &PfromX = PfromX_by_mtype(type);
		PfromX.set(x.i(),x.j(),x.k(),x.l(),min_energy);
	}
}

void pseudo_loop::compute_PfromR(const Index4D &x, MType type){
	energy_t min_energy = INF;

	Matrix4D &PfromXprime = PfromXprime_by_mtype(type);
	for(cand_pos_t d=x.k()+1; d<=x.l(); ++d){
		energy_t tmp = calc_WP(x.k(),d-1) + PfromXprime.get(x.i(),x.j(),d,x.l());
		min_energy = std::min(min_energy,tmp);
	}

	if (min_energy < INF/2){
		Matrix4D &PfromX = PfromX_by_mtype(type);
		PfromX.set(x.i(),x.j(),x.k(),x.l(),min_energy);
	}
}

void pseudo_loop::compute_PfromO(const Index4D &x, MType type){
	energy_t min_energy = INF;

	Matrix4D &PfromXprime = PfromXprime_by_mtype(type);
	for(cand_pos_t d=x.i()+1; d<=x.j(); ++d){
		energy_t tmp = calc_WP(x.i(),d-1) + PfromXprime.get(d,x.j(),x.k(),x.l());
		min_energy = std::min(min_energy,tmp);
	}

	if (min_energy < INF/2){
		Matrix4D &PfromX = PfromX_by_mtype(type);
		PfromX.set(x.i(),x.j(),x.k(),x.l(),min_energy);
	}
}
/**
 * 
 * 
 * 
 */
void pseudo_loop::compute_PfromLprime(const Index4D &x, MType type){
	energy_t min_energy = INF;

	for(cand_pos_t d=x.i(); d<x.j(); ++d){
		energy_t tmp = calc_PfromXdoubleprime(x.i(),d,x.k(),x.l(),type) + calc_WP(d+1,x.j());
		min_energy = std::min(min_energy,tmp);
	}
	
	if (min_energy < INF/2){
		Matrix4D &PfromXprime = PfromXprime_by_mtype(type);
		PfromXprime.set(x.i(),x.j(),x.k(),x.l(),min_energy);
	}
}

void pseudo_loop::compute_PfromMprime(const Index4D &x, MType type){
	energy_t min_energy = INF;

	for(cand_pos_t d=x.k()+1; d<=x.l(); ++d){
		energy_t tmp= calc_WP(x.k(),d-1) + calc_PfromXdoubleprime(x.i(),x.j(),d,x.l(),type);
		min_energy = std::min(min_energy,tmp);
	}
	
	if (min_energy < INF/2){
		Matrix4D &PfromXprime = PfromXprime_by_mtype(type);
		PfromXprime.set(x.i(),x.j(),x.k(),x.l(),min_energy);
	}
}

void pseudo_loop::compute_PfromRprime(const Index4D &x, MType type){
	energy_t min_energy = INF;

	for(cand_pos_t d=x.k(); d<x.l(); ++d){
		energy_t tmp = calc_PfromXdoubleprime(x.i(),x.j(),x.k(),d,type) + calc_WP(d+1,x.l());
		min_energy = std::min(min_energy,tmp);
	}
	
	if (min_energy < INF/2){
		Matrix4D &PfromXprime = PfromXprime_by_mtype(type);
		PfromXprime.set(x.i(),x.j(),x.k(),x.l(),min_energy);
	}
}

void pseudo_loop::compute_PfromOprime(const Index4D &x, MType type){
	energy_t min_energy = INF;

	for(cand_pos_t d=x.k(); d<x.l(); ++d){
		energy_t tmp = calc_PfromXdoubleprime(x.i(),x.j(),x.k(),d,type) + calc_WP(d+1,x.j());
		min_energy = std::min(min_energy,tmp);
	}
	
	if (min_energy < INF/2){
		Matrix4D &PfromXprime = PfromXprime_by_mtype(type);
		PfromXprime.set(x.i(),x.j(),x.k(),x.l(),min_energy);
	}
}

void pseudo_loop::compute_PLmloop00(const Index4D &x, MType type){
	energy_t min_energy = INF,tmp=INF;
	const cand_pos_t i = x.i(), j = x.j(), k = x.k(), l = x.l();

	Matrix4D &PXmloop10 = PXmloop10_by_mtype(type);
	Matrix4D &PXmloop01 = PXmloop01_by_mtype(type);
    for(cand_pos_t d = i; d<=j; ++d){
        if (d>i){
            tmp = calc_WB(i,d-1) + PXmloop10.get(d,j,k,l);
            min_energy = std::min(min_energy,tmp);
        }
        if(d<j){
            tmp = PXmloop01.get(i,d,k,l) + calc_WB(d+1,j);
            min_energy = std::min(min_energy,tmp);
        }
    }
	if (min_energy < INF/2){
		Matrix4D &PXmloop00 = PXmloop00_by_mtype(type);
		PXmloop00.set(i,j,k,l,min_energy);
	}
}
void pseudo_loop::compute_PMmloop00(const Index4D &x, MType type){
	energy_t min_energy = INF,tmp=INF;
	const cand_pos_t i = x.i(), j = x.j(), k = x.k(), l = x.l();

	Matrix4D &PXmloop10 = PXmloop10_by_mtype(type);
	Matrix4D &PXmloop01 = PXmloop01_by_mtype(type);
    for(cand_pos_t d=i; d<j; ++d){
        tmp=PXmloop10.get(i,d,k,l) + calc_WB(d+1,j);
        min_energy = std::min(min_energy,tmp);
    }
    for(cand_pos_t d=k+1; d<=l; ++d){
        tmp=PXmloop01.get(i,j,d,l) + calc_WB(k,d-1);
        min_energy = std::min(min_energy,tmp);
    }
	if (min_energy < INF/2){
		Matrix4D &PXmloop00 = PXmloop00_by_mtype(type);
		PXmloop00.set(i,j,k,l,min_energy);
	}
}
void pseudo_loop::compute_PRmloop00(const Index4D &x, MType type){
	energy_t min_energy = INF,tmp=INF;
	const cand_pos_t i = x.i(), j = x.j(), k = x.k(), l = x.l();

	Matrix4D &PXmloop10 = PXmloop10_by_mtype(type);
	Matrix4D &PXmloop01 = PXmloop01_by_mtype(type);
	for(cand_pos_t d=k; d<=l; ++d){
        if(d>k){
            tmp=calc_WB(k,d-1) + PXmloop10.get(i,j,d,l);
            min_energy = std::min(min_energy,tmp);
        }
        if (d<l){
            tmp = PXmloop01.get(i,j,k,d)+calc_WB(d+1,l);
            min_energy = std::min(min_energy,tmp);
        }
    }
	if (min_energy < INF/2){
		Matrix4D &PXmloop00 = PXmloop00_by_mtype(type);
		PXmloop00.set(i,j,k,l,min_energy);
	}
}
void pseudo_loop::compute_POmloop00(const Index4D &x, MType type){
	energy_t min_energy = INF,tmp=INF;
	const cand_pos_t i = x.i(), j = x.j(), k = x.k(), l = x.l();

	Matrix4D &PXmloop10 = PXmloop10_by_mtype(type);
	Matrix4D &PXmloop01 = PXmloop01_by_mtype(type);
    for(cand_pos_t d=i+1; d<=j; ++d){
        tmp = calc_WB(i,d-1)+PXmloop10.get(d,j,k,l);
        min_energy = std::min(min_energy,tmp);
    }
    for(cand_pos_t d=k; d<l; ++d){
        tmp = PXmloop01.get(i,j,k,d)+calc_WB(d+1,l);
        min_energy = std::min(min_energy,tmp);
    }
	if (min_energy < INF/2){
		Matrix4D &PXmloop00 = PXmloop00_by_mtype(type);
		PXmloop00.set(i,j,k,l,min_energy);
	}
}
/**
 * 
 * 
 * 
 */
void pseudo_loop::compute_PLmloop10(const Index4D &x, MType type){
	energy_t min_energy = INF,tmp=INF;
	const cand_pos_t i = x.i(), j = x.j(), k = x.k(), l = x.l();

	Matrix4D &PX = PX_by_mtype(type);
	for(cand_pos_t d = i+1; d <= j; ++d){
        tmp = PX.get(i,d,k,l) + calc_WB(d+1,j);
        min_energy = std::min(min_energy,tmp);
    }
	if (min_energy < INF/2){
		Matrix4D &PXmloop10 = PXmloop10_by_mtype(type);
		PXmloop10.set(i,j,k,l,min_energy);
	}

}
void pseudo_loop::compute_PMmloop10(const Index4D &x, MType type){
	energy_t min_energy = INF,tmp=INF;
	const cand_pos_t i = x.i(), j = x.j(), k = x.k(), l = x.l();

	Matrix4D &PX = PX_by_mtype(type);
    for(cand_pos_t d = i+1; d <=j; ++d){
        tmp = PX.get(i,d,k,d) + calc_WB(d+1,j);
        min_energy = std::min(min_energy,tmp);
    }

	if (min_energy < INF/2){
		Matrix4D &PXmloop10 = PXmloop10_by_mtype(type);
		PXmloop10.set(i,j,k,l,min_energy);
	}
}
void pseudo_loop::compute_PRmloop10(const Index4D &x, MType type){
	energy_t min_energy = INF,tmp=INF;
	const cand_pos_t i = x.i(), j = x.j(), k = x.k(), l = x.l();

	Matrix4D &PX = PX_by_mtype(type);
    for(cand_pos_t d = k+1; d <= l; ++d){
        tmp = PX.get(i,j,d,l) + calc_WB(d+1,l);
        min_energy = std::min(min_energy,tmp);
    }
	if (min_energy < INF/2){
		Matrix4D &PXmloop10 = PXmloop10_by_mtype(type);
		PXmloop10.set(i,j,k,l,min_energy);
	}

}
void pseudo_loop::compute_POmloop10(const Index4D &x, MType type){
	energy_t min_energy = INF,tmp=INF;
	const cand_pos_t i = x.i(), j = x.j(), k = x.k(), l = x.l();

	Matrix4D &PX = PX_by_mtype(type);
	for(cand_pos_t d=i+1; d<=j;++d){
        tmp=PX.get(d,j,k,l) + calc_WB(d+1,l);
        min_energy = std::min(min_energy,tmp);
    }
	if (min_energy < INF/2){
		Matrix4D &PXmloop10 = PXmloop10_by_mtype(type);
		PXmloop10.set(i,j,k,l,min_energy);
	}
}
 /**
  * 
  * 
  * 
  */
void pseudo_loop::compute_PLmloop01(const Index4D &x, MType type){
	energy_t min_energy = INF,tmp=INF;
	const cand_pos_t i = x.i(), j = x.j(), k = x.k(), l = x.l();

	Matrix4D &PX = PX_by_mtype(type);
	for(cand_pos_t d = i; d < j; ++d){
        tmp = cp_penalty*(d-i) + PX.get(d,j,k,l);
        min_energy = std::min(min_energy,tmp);
    }
	if (min_energy < INF/2){
		Matrix4D &PXmloop01 = PXmloop01_by_mtype(type);
		PXmloop01.set(i,j,k,l,min_energy);
	}

}
void pseudo_loop::compute_PMmloop01(const Index4D &x, MType type){
	energy_t min_energy = INF,tmp=INF;
	const cand_pos_t i = x.i(), j = x.j(), k = x.k(), l = x.l();

	Matrix4D &PX = PX_by_mtype(type);
	for(cand_pos_t d = k; d < l; ++d){
        tmp = cp_penalty*(d-k) + PX.get(i,j,d,l);
        min_energy = std::min(min_energy,tmp);
    }

	if (min_energy < INF/2){
		Matrix4D &PXmloop01 = PXmloop01_by_mtype(type);
		PXmloop01.set(i,j,k,l,min_energy);
	}
}
void pseudo_loop::compute_PRmloop01(const Index4D &x, MType type){
	energy_t min_energy = INF,tmp=INF;
	const cand_pos_t i = x.i(), j = x.j(), k = x.k(), l = x.l();

	Matrix4D &PX = PX_by_mtype(type);
    for(cand_pos_t d = k; d < l; d++){
        tmp = cp_penalty*(d-k) + PX.get(i,j,d,l);
        min_energy = std::min(min_energy,tmp);
    }
	if (min_energy < INF/2){
		Matrix4D &PXmloop01 = PXmloop01_by_mtype(type);
		PXmloop01.set(i,j,k,l,min_energy);
	}
}
void pseudo_loop::compute_POmloop01(const Index4D &x, MType type){
	energy_t min_energy = INF;
	const cand_pos_t i = x.i(), j = x.j(), k = x.k(), l = x.l();

	Matrix4D &PX = PX_by_mtype(type);
	for(cand_pos_t d = k; d < l; ++d){
        energy_t tmp = cp_penalty*(d-i) + PX.get(d,j,k,l);
        min_energy = std::min(min_energy,tmp);
    }

	if (min_energy < INF/2){
		Matrix4D &PXmloop01 = PXmloop01_by_mtype(type);
		PXmloop01.set(i,j,k,l,min_energy);
	}
}
/**
 * 
 * 
 * 
 */
energy_t pseudo_loop::calc_PfromLdoubleprime(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l){
	if (!(i <= j && j < k-1 && k <= l)){
		return INF;
	}
	assert(!(i<=0 || l> n));

	energy_t b1 = PR.get(i,j,k,l) + gamma2(l,k) + PB_penalty;
	energy_t b2 = PM.get(i,j,k,l) + gamma2(j,k) + PB_penalty;
	energy_t b3 = calc_PO(i,j,k,l) + gamma2(l,i) + PB_penalty;

	return std::min({b1,b2,b3});
}

energy_t pseudo_loop::calc_PfromRdoubleprime(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l){
	if (!(i <= j && j < k-1 && k <= l)){
		return INF;
	}
	assert(!(i<=0 || l> n));

	energy_t b1 = PM.get(i,j,k,l) + gamma2(j,k) + PB_penalty;
	energy_t b2 = calc_PO(i,j,k,l) + gamma2(l,i) + PB_penalty;

	return std::min({b1,b2});
}

energy_t pseudo_loop::calc_PfromMdoubleprime(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l){
	if (!(i <= j && j < k-1 && k <= l)){
		return INF;
	}
	assert(!(i<=0 || l> n));

	energy_t b1 = PL.get(i,j,k,l) + gamma2(j,i) + PB_penalty;
	energy_t b2 = PR.get(i,j,k,l) + gamma2(l,k) + PB_penalty;
	energy_t b3 = POm.get(i,j,k,l) + gamma2(l,i) + PB_penalty;

	return std::min({b1,b2,b3});
}

energy_t pseudo_loop::calc_PfromOdoubleprime(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l){
	if (!(i <= j && j < k-1 && k <= l)){
		return INF;
	}
	assert(!(i<=0 || l> n));

	energy_t b1 = PL.get(i,j,k,l) + gamma2(j,i) + PB_penalty;
	energy_t b2 = PR.get(i,j,k,l) + gamma2(l,k) + PB_penalty;

	return std::min({b1,b2});
}

energy_t pseudo_loop::calc_PfromLreOdoubleprime(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l){
	if (!(i <= j && j < k-1 && k <= l)){
		return INF;
	}
	assert(!(i<=0 || l> n));

	energy_t b1 = PMreO.get(i,j,k,l) + gamma2(j,k) + PB_penalty;
	energy_t b2 = calc_PO(i,j,k,l) + gamma2(l,i) + PB_penalty;

	return std::min({b1,b2});
}

energy_t pseudo_loop::calc_PfromLreRdoubleprime(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l){
	if (!(i <= j && j < k-1 && k <= l)){
		return INF;
	}
	assert(!(i<=0 || l> n));

	energy_t b1 = PMreR.get(i,j,k,l) + gamma2(j,k) + PB_penalty;
	energy_t b2 = PR.get(i,j,k,l) + gamma2(l,k) + PB_penalty;

	return std::min({b1,b2});
}

energy_t pseudo_loop::calc_PfromMreOdoubleprime(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l){
	if (!(i <= j && j < k-1 && k <= l)){
		return INF;
	}
	assert(!(i<=0 || l> n));

	energy_t b1 = PLreO.get(i,j,k,l) + gamma2(j,i) + PB_penalty;
	energy_t b2 = calc_PO(i,j,k,l) + gamma2(l,i) + PB_penalty;

	return std::min({b1,b2});
}

energy_t pseudo_loop::calc_PfromMreRdoubleprime(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l){
	if (!(i <= j && j < k-1 && k <= l)){
		return INF;
	}
	assert(!(i<=0 || l> n));

	energy_t b1 = PLreR.get(i,j,k,l) + gamma2(j,i) + PB_penalty;
	energy_t b2 = PR.get(i,j,k,l) + gamma2(l,k) + PB_penalty;

	return std::min({b1,b2});
}

energy_t pseudo_loop::calc_PfromLMreRdoubleprime(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l){
	if (!(i <= j && j < k-1 && k <= l)){
		return INF;
	}
	assert(!(i<=0 || l> n));

	energy_t b1 = PMreR.get(i,j,k,l) + gamma2(j,k) + PB_penalty;

	return b1;
}

energy_t pseudo_loop::calc_PfromLMorOdoubleprime(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l){
	if (!(i <= j && j < k-1 && k <= l)){
		return INF;
	}
	assert(!(i<=0 || l> n));

	energy_t b1 = PM.get(i,j,k,l) + gamma2(j,k) + PB_penalty;
	energy_t b2 = calc_PO(i,j,k,l) + gamma2(l,i) + PB_penalty;

	return std::min(b1,b2);
}

energy_t pseudo_loop::calc_WB(cand_pos_t i, cand_pos_t l){
	if (i<=0 || l<=0 || i>n || l>n){
		return INF;
	}
	if (i>l) return 0;
	return (std::min(cp_penalty*(l-i+1),WBP.get(i,l)));
}

energy_t pseudo_loop::calc_WP(cand_pos_t i, cand_pos_t l){
	if (i<=0 || l<=0 || i>n || l>n){
		return INF;
	}
	if (i>l) return 0; // needed as this will happen 
	return (std::min(PUP_penalty*(l-i+1),WPP.get(i,l)));
}


////////////////////////// Traceback ////////////////////////////////

void pseudo_loop::Trace_WB(cand_pos_t i, cand_pos_t l, energy_t e){
	if (debug) printf("WB at %d and %d with %d\n", i, l, e);
	if (i>l) return;
	if(e==cp_penalty*(l-i+1)) return;
	if(e==WBP.get(i,l)){
		Trace_WBP(i,l,WBP.get(i,l));
	}
	__builtin_unreachable();
}
void pseudo_loop::Trace_WBP(cand_pos_t i, cand_pos_t l, energy_t e){
	if (debug) printf("WBP at %d and %d with %d\n", i, l, e);
	energy_t tmp = WBP.get(i,l-1)+cp_penalty;
	if(e==tmp){
		Trace_WBP(i,l-1,WBP.get(i,l-1));
		return;
	}
	for(cand_pos_t d=i; d< l; ++d){
		tmp = calc_WB(i,d-1) + V->get_energy(d,l) + bp_penalty + PPS_penalty;
		if(e==tmp){
			Trace_WB(i,d-1,calc_WB(i,d-1));
			W->Trace_V(d,l,V->get_energy(d,l));
			return;
		}
		tmp = calc_WB(i,d-1) + P.get(d,l) + PSM_penalty + PPS_penalty;
		if(e==tmp){
			Trace_WB(i,d-1,calc_WB(i,d-1));
			Trace_P(d,l,P.get(d,l));
			return;
		}
	}
	__builtin_unreachable();
}
void pseudo_loop::Trace_WP(cand_pos_t i, cand_pos_t l, energy_t e){
	if (debug) printf("WP at %d and %d with %d\n", i, l, e);
	if (i>l) return;
	if(e==PUP_penalty*(l-i+1)) return;
	if(e==WPP.get(i,l)){
		Trace_WBP(i,l,WPP.get(i,l));
	}
	__builtin_unreachable();
}
void pseudo_loop::Trace_WPP(cand_pos_t i, cand_pos_t l, energy_t e){
	if (debug) printf("WPP at %d and %d with %d\n", i, l, e);
	energy_t tmp = WPP.get(i,l-1)+PUP_penalty;
	if(e==tmp){
		Trace_WPP(i,l-1,WPP.get(i,l-1));
		return;
	}
	for(cand_pos_t d=i; d<l; ++d){
		tmp = calc_WP(i,d-1) + V->get_energy(d,l) + gamma2(l,d) + PPS_penalty;
		if(e==tmp){
			Trace_WP(i,d-1,calc_WP(i,d-1));
			W->Trace_V(d,l,V->get_energy(d,l));
			return;
		}
		tmp = calc_WP(i,d-1) + P.get(d,l) + PSP_penalty + PPS_penalty;
		if(e==tmp){
			Trace_WP(i,d-1,calc_WP(i,d-1));
			Trace_P(d,l,P.get(d,l));
			return;
		}
	}
	__builtin_unreachable();
}

void pseudo_loop::Trace_P(cand_pos_t i, cand_pos_t l, energy_t e){
	if (debug) printf("P at %d and %d with %d\n", i, l, e);
	energy_t tmp = INF;
	for(cand_pos_t j=i; j<l; ++j){
		for (cand_pos_t d=j+1; d<l; ++d){
			for (cand_pos_t k=d+1; k<l; ++k){
				// b1 = PK1Om.get(i,j,d+1,k) + PK2Om.get(j+1,d,k+1,l);
				// b2 = PK1Om.get(i,j,d+1,k) + PK2Os.get(j+1,d,k+1,l);
				// b3 = PK1Om.get(i,j,d+1,k) + PK2LreO.get(j+1,d,k+1,l);
				// b4 = PK1Om.get(i,j,d+1,k) + PK2MreO.get(j+1,d,k+1,l);
				// b5 = PK1Om.get(i,j,d+1,k) + PK2LreR.get(j+1,d,k+1,l);
				// b6 = PK1Om.get(i,j,d+1,k) + PK2MreR.get(j+1,d,k+1,l);
				// b7 = PK1Om.get(i,j,d+1,k) + PK2R.get(j+1,d,k+1,l);

				// b8 = PK1Os.get(i,j,d+1,k) + PK2Om.get(j+1,d,k+1,l);
				tmp = PK1Os.get(i,j,d+1,k) + PK2Os.get(j+1,d,k+1,l);
				if(e==tmp){
					Trace_PX1(i,j,d+1,k,MType::Os,PK1Os.get(i,j,d+1,k));
					Trace_PX2(j+1,d,k+1,l,MType::Os,PK2Os.get(j+1,d,k+1,l));
					return;
				}
				// b10 = PK1Os.get(i,j,d+1,k) + PK2MreO.get(j+1,d,k+1,l);
				// b11 = PK1Os.get(i,j,d+1,k) + PK2MreR.get(j+1,d,k+1,l);
				// b12 = PK1Os.get(i,j,d+1,k) + PK2R.get(j+1,d,k+1,l);

				// b13 = PK1LreO.get(i,j,d+1,k) + PK2Om.get(j+1,d,k+1,l);
				// b14 = PK1LreO.get(i,j,d+1,k) + PK2LreO.get(j+1,d,k+1,l);
				// b15 = PK1LreO.get(i,j,d+1,k) + PK2MreO.get(j+1,d,k+1,l);

				// b16 = PK1LMreR.get(i,j,d+1,k) + PK2Om.get(j+1,d,k+1,l);
				// b17 = PK1LMreR.get(i,j,d+1,k) + PK2LreO.get(j+1,d,k+1,l);
				// b18 = PK1LMreR.get(i,j,d+1,k) + PK2MreO.get(j+1,d,k+1,l);

				// b19 = PK1LMorO.get(i,j,d+1,k) + PK2LreR.get(j+1,d,k+1,l);
				// b20 = PK1LMorO.get(i,j,d+1,k) + PK2MreO.get(j+1,d,k+1,l);
				// b21 = PK1LMorO.get(i,j,d+1,k) + PK2MreO.get(j+1,d,k+1,l);
			}
		}
	}
	__builtin_unreachable();
}
void pseudo_loop::Trace_PX1(cand_pos_t i,cand_pos_t j,cand_pos_t k, cand_pos_t l, MType type, energy_t e){
	if (debug) std::cout << "PX1 at " << i << " and " << j << " and " << k << " and " << l << " with type: " << type << " and en: " << e << std::endl;
	Matrix4D &PX = PX_by_mtype(type);
	energy_t tmp = INF;
	for(cand_pos_t d = i+1;d<=j;++d){
		tmp = PX.get(i,d,k,l) + calc_WP(d+1,j) + PB_penalty;
		if(e==tmp){
			Trace_PX(i,d,k,l,type,PX.get(i,d,k,l));
			Trace_WP(d+1,j,calc_WP(d+1,j));
			return;
		}
	}
	__builtin_unreachable();
}
void pseudo_loop::Trace_PX2(cand_pos_t i,cand_pos_t j,cand_pos_t k, cand_pos_t l, MType type, energy_t e){
	if (debug) std::cout << "PX2 at " << i << " and " << j << " and " << k << " and " << l << " with type: " << type << " and en: " << e << std::endl;
	Matrix4D &PK1X = PK1X_by_mtype(type);
	energy_t tmp = INF;
	for(cand_pos_t d = k;d<l;++d){
		tmp = PK1X.get(i,j,d,l) + calc_WP(k,d-1);
		if(e==tmp){
			Trace_PX1(i,j,d,l,type,PK1X.get(i,j,d,l));
			Trace_WP(k,d-1,calc_WP(k,d-1));
			return;
		}
	}
	__builtin_unreachable();
}
void pseudo_loop::Trace_PX(cand_pos_t i,cand_pos_t j,cand_pos_t k, cand_pos_t l, MType type, energy_t e){
	if (debug) std::cout << "PX at " << i << " and " << j << " and " << k << " and " << l << " with type: " << type << " and en: " << e << std::endl;
	const Index4D x(i,j,k,l);
	const int ptype_closing = can_pair(x.lend(type),x.rend(type));//pair[S_[x.lend(type)]][S_[x.rend(type)]];
	f[x.lend(type)].pair = x.rend(type);
	f[x.rend(type)].pair = x.lend(type);
	// std::cout << x.lend(type) << " " << x.rend(type) << " " << S_[x.lend(type)] << "\t" << S_[x.rend(type)] << pair[S_[x.lend(type)]][S_[x.rend(type)]] << "\t" << ptype_closing << std::endl;
	// exit(0);
	if (ptype_closing>0){
		if(type == MType::Os){
			if(x.i()==x.j() && x.k()==x.l()){
				if(e==gamma2(x.l(),x.i())) return;
			}
		} else if (x.difference(type)>TURN){
			Index4D xp(x);
			// xp.shrink(type);
			// Matrix4D &PfromX = PfromX_by_mtype(type);
			// b3 = PfromX.get(xp) + penalty(xp, gamma2, type);
		}
		// energy_t tmp = calc_PXmloop(x,type);
		// if(e==tmp){
		// 	Trace_PXmloop(i,j,k,l,type,tmp);
		// 	return;
		// }
		energy_t tmp = calc_PXiloop(x, type);
		if(e==tmp){
			Trace_PXiloop(x,type,tmp);
			return;
		}
	}
	__builtin_unreachable();
}

void pseudo_loop::Trace_PLiloop(const Index4D &x, MType type, energy_t e){
	if (debug) std::cout << "PLiloop at " << x.i() << " and " << x.j() << " and " << x.k() << " and " << x.l() << " with type: " << type << " and en: " << e << std::endl;
	assert(!impossible_case(x));
	const cand_pos_t i = x.i(), j = x.j(), k = x.k(), l = x.l();
	Matrix4D &PX = PX_by_mtype(type);
	energy_t tmp = INF;
	if (i+TURN+2<j) { 
		tmp = PX.get(i+1,j-1,k,l) + get_e_stP(i,j);
		if(e==tmp){
			Trace_PX(i+1,j-1,k,l,type,PX.get(i+1,j-1,k,l));
			return;
		}
	}
	cand_pos_t max_d = std::min(j,i+MAXLOOP);
	for(cand_pos_t d= i+1; d<max_d; ++d){
		cand_pos_t min_dp = std::max(d+TURN,j-MAXLOOP);
		for(cand_pos_t dp = j-1; dp > min_dp; --dp){
			if (!can_pair(d,dp)) continue;
			tmp = get_e_intP(i,d,dp,j) + PX.get(d,dp,k,l);
			if(e==tmp){
				Trace_PX(d,dp,k,l,type,PX.get(d,dp,k,l));
				return;
			}
		}
	}
}
void pseudo_loop::Trace_PMiloop(const Index4D &x, MType type, energy_t e){
	if (debug) std::cout << "PMiloop at " << x.i() << " and " << x.j() << " and " << x.k() << " and " << x.l() << " with type: " << type << " and en: " << e << std::endl;
	assert(!impossible_case(x));
	const cand_pos_t i = x.i(), j = x.j(), k = x.k(), l = x.l();
	Matrix4D &PX = PX_by_mtype(type);
	energy_t tmp = INF;
	if (i<j && k<l) {
		tmp = PX.get(i,j-1,k+1,l) + get_e_stP(j-1,k+1);
		if(e==tmp){
			Trace_PX(i,j-1,k+1,l,type,PX.get(i,j-1,k+1,l));
			return;
		}
	}
	cand_pos_t max_d = std::max(i,j-MAXLOOP);
	for(cand_pos_t d= j-1; d>max_d; --d){
		cand_pos_t min_dp = std::min(l,k+MAXLOOP); // could switch these here so that we are increasing in the first for like all the others
		for (cand_pos_t dp=k+1; dp <min_dp; ++dp) {
			if (!can_pair(d,dp)) continue;
			tmp = get_e_intP(d,j,k,dp) + PX.get(i,d,dp,l);
			if(e==tmp){
				Trace_PX(i,d,dp,l,type,PX.get(i,d,dp,l));
				return;
			}
		}
	}

}
void pseudo_loop::Trace_PRiloop(const Index4D &x, MType type, energy_t e){
	if (debug) std::cout << "PRiloop at " << x.i() << " and " << x.j() << " and " << x.k() << " and " << x.l() << " with type: " << type << " and en: " << e << std::endl;
	assert(!impossible_case(x));
	const cand_pos_t i = x.i(), j = x.j(), k = x.k(), l = x.l();
	Matrix4D &PX = PX_by_mtype(type);
	energy_t tmp = INF;
	if (k+TURN+2<l) { 
		tmp = PX.get(i,j,k+1,l-1) + get_e_stP(k,l);
		if(e==tmp){
			Trace_PX(i,j,k+1,l-1,type,PX.get(i,j,k+1,l-1));
			return;
		}
	}
	cand_pos_t max_d = std::min(l,k+MAXLOOP);
	for(cand_pos_t d= k+1; d<max_d; ++d){
		cand_pos_t min_dp = std::max(d+TURN,l-MAXLOOP);
		for(cand_pos_t dp=l-1; dp > min_dp; --dp){
			if (!can_pair(d,dp)) continue;
			tmp = get_e_intP(k,d,dp,l) + PX.get(i,j,d,dp);
			if(e==tmp){
				Trace_PX(i,j,d,dp,type,PX.get(i,j,d,dp));
				return;
			}
		}
	}

}
void pseudo_loop::Trace_POiloop(const Index4D &x, MType type, energy_t e){
	if (debug) std::cout << "POiloop at " << x.i() << " and " << x.j() << " and " << x.k() << " and " << x.l() << " with type: " << type << " and en: " << e << std::endl;
	assert(!impossible_case(x));
	const cand_pos_t i = x.i(), j = x.j(), k = x.k(), l = x.l();
	Matrix4D &PX = PX_by_mtype(type);
	energy_t tmp = INF;
	if (i<j && k<l ) { 
		tmp = PX.get(i+1,j,k,l-1) + get_e_stP(i,l);
		if(e==tmp){
			Trace_PX(i+1,j,k,l-1,type,PX.get(i+1,j,k,l-1));
			return;
		}
	}
	cand_pos_t max_d = std::min(j,i+MAXLOOP);
	for(cand_pos_t d= i+1; d<max_d; ++d){
		cand_pos_t min_dp = std::max(l-MAXLOOP,k);
		for (cand_pos_t dp=l-1; dp >min_dp; --dp) {
			if (!can_pair(d,dp)) continue;
			tmp = get_e_intP(i,d,dp,l) + PX.get(d,j,k,dp);
			if(e==tmp){
				Trace_PX(d,j,k,dp,type,PX.get(d,j,k,dp));
				return;
			}
		}
	}

}