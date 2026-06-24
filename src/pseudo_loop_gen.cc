#include "pseudo_loop.hh"
#include "h_externs.hh"
#include <stdio.h>
#include <string>
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <algorithm>
#include <cassert>

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
			min_energy = std::min(min_energy,get_e_intP(i,d,dp,l) + PX.get(d,j,dp,k));
		}
	}
	return min_energy;
}

void pseudo_loop::compute_PfromL(const Index4D &x, MType type){
	energy_t min_energy = INF;

	Matrix4D &PfromXprime = PfromXprime_by_mtype(type);
	for(cand_pos_t d=x.i()+1; d<=x.j(); ++d){
		energy_t tmp = get_WP(x.i(),d-1) + PfromXprime.get(d,x.j(),x.k(),x.l());
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
		energy_t tmp = PfromXprime.get(x.i(),d,x.k(),x.l())+get_WP(d+1,x.j());
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
		energy_t tmp = get_WP(x.k(),d-1) + PfromXprime.get(x.i(),x.j(),d,x.l());
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
		energy_t tmp = get_WP(x.i(),d-1) + PfromXprime.get(d,x.j(),x.k(),x.l());
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

// These are not fully generic as I have to do the getters for  PfromXdoubleprime
// void pseudo_loop::compute_PfromLprime(const Index4D &x, MType type){
// 	energy_t min_energy = INF;

// 	for(cand_pos_t d=x.i(); d<x.j(); ++d){
// 		energy_t tmp=get_PfromLdoubleprime(x.i(),d,x.k(),x.l())+get_WP(d+1,x.j());
// 		min_energy = std::min(min_energy,tmp);
// 	}
	
// 	if (min_energy < INF/2){
// 		Matrix4D &PfromXprime = PfromXprime_by_mtype(type);
// 		PfromXprime.set(x.i(),x.j(),x.k(),x.l(),min_energy);
// 	}
// }

// void pseudo_loop::compute_PfromMprime(const Index4D &x, MType type){
// 	energy_t min_energy = INF;

// 	for(cand_pos_t d=x.k()+1; d<=x.l(); ++d){
// 		energy_t tmp= get_WP(x.k(),d-1) + get_PfromMdoubleprime(x.i(),x.j(),d,x.l());
// 		min_energy = std::min(min_energy,tmp);
// 	}
	
// 	if (min_energy < INF/2){
// 		Matrix4D &PfromXprime = PfromXprime_by_mtype(type);
// 		PfromXprime.set(x.i(),x.j(),x.k(),x.l(),min_energy);
// 	}
// }

// void pseudo_loop::compute_PfromRprime(const Index4D &x, MType type){
// 	energy_t min_energy = INF;

// 	for(cand_pos_t d=x.k(); d<x.l(); ++d){
// 		energy_t tmp=get_PfromRdoubleprime(x.i(),x.j(),x.k(),d)+get_WP(d+1,x.l());
// 		min_energy = std::min(min_energy,tmp);
// 	}
	
// 	if (min_energy < INF/2){
// 		Matrix4D &PfromXprime = PfromXprime_by_mtype(type);
// 		PfromXprime.set(x.i(),x.j(),x.k(),x.l(),min_energy);
// 	}
// }

// void pseudo_loop::compute_PfromOprime(const Index4D &x, MType type){
// 	energy_t min_energy = INF;

// 	for(cand_pos_t d=x.k(); d<x.l(); ++d){
// 		energy_t tmp=get_PfromOdoubleprime(x.i(),x.j(),x.k(),d)+get_WP(d+1,x.j());
// 		min_energy = std::min(min_energy,tmp);
// 	}
	
// 	if (min_energy < INF/2){
// 		Matrix4D &PfromXprime = PfromXprime_by_mtype(type);
// 		PfromXprime.set(x.i(),x.j(),x.k(),x.l(),min_energy);
// 	}
// }

void pseudo_loop::compute_PLmloop00(const Index4D &x, MType type){
	energy_t min_energy = INF,tmp=INF;
	const cand_pos_t i = x.i(), j = x.j(), k = x.k(), l = x.l();

	Matrix4D &PXmloop10 = PXmloop10_by_mtype(type);
	Matrix4D &PXmloop01 = PXmloop01_by_mtype(type);
    for(cand_pos_t d = i; d<=j; ++d){
        if (d>i){
            tmp=get_WB(i,d-1)+PXmloop10.get(d,j,k,l);
            min_energy = std::min(min_energy,tmp);
        }
        if(d<j){
            tmp=PXmloop01.get(i,d,k,l)+get_WB(d+1,j);
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
        tmp=PXmloop10.get(i,d,k,l)+get_WB(d+1,j);
        min_energy = std::min(min_energy,tmp);
    }
    for(cand_pos_t d=k+1; d<=l; ++d){
        tmp=PXmloop01.get(i,j,d,l)+get_WB(k,d-1);
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
            tmp=get_WB(k,d-1)+PXmloop10.get(i,j,d,l);
            min_energy = std::min(min_energy,tmp);
        }
        if (d<l){
            tmp = PXmloop01.get(i,j,k,d)+get_WB(d+1,l);
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
        tmp = get_WB(i,d-1)+PXmloop10.get(d,j,k,l);
        min_energy = std::min(min_energy,tmp);
    }
    for(cand_pos_t d=k; d<l; ++d){
        tmp = PXmloop01.get(i,j,k,d)+get_WB(d+1,l);
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
        tmp = PX.get(i,d,k,l) + get_WB(d+1,j);
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
        tmp = PX.get(i,d,k,d) + get_WB(d+1,j);
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
        tmp = PX.get(i,j,d,l) + get_WB(d+1,l);
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
        tmp=PX.get(d,j,k,l) + get_WB(d+1,l);
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
// energy_t pseudo_loop::calc_PLiloop(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
// 	if (!(i <= j && j < k-1 && k <= l)){
// 		return INF;
// 	}
// 	assert(!(i<=0 || l> n));

// 	if(!can_pair(i,j)) return INF;

// 	energy_t min_energy = INF;
// 	if (i+TURN+2<j) { 
// 		min_energy = PL.get(i+1,j-1,k,l) + get_e_stP(i,j);
// 	}
// 	cand_pos_t max_d = std::min(j,i+MAXLOOP);
// 	for(cand_pos_t d= i+1; d<max_d; ++d){
// 		cand_pos_t min_dp = std::max(d+TURN,j-MAXLOOP);
// 		for(cand_pos_t dp = j-1; dp > min_dp; --dp){
// 			if (!can_pair(d,dp)) continue;
// 			min_energy = std::min(min_energy,get_e_intP(i,d,dp,j) + PL.get(d,dp,k,l));
// 		}
// 	}
// 	return min_energy;
// }

// energy_t pseudo_loop::calc_PRiloop(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
// 	if (!(i <= j && j < k-1 && k <= l)){
// 		return INF;
// 	}
// 	assert(!(i<=0 || l> n));

// 	if(!can_pair(k,l)) return INF;

// 	energy_t min_energy = INF;
// 	if (k+TURN+2<l) { 
// 		min_energy = PR.get(i,j,k+1,l-1) + get_e_stP(k,l);
// 	}
// 	cand_pos_t max_d = std::min(l,k+MAXLOOP);
// 	for(cand_pos_t d= k+1; d<max_d; ++d){
// 		cand_pos_t min_dp = std::max(d+TURN,l-MAXLOOP);
// 		for(cand_pos_t dp=l-1; dp > min_dp; --dp){
// 			if (!can_pair(d,dp)) continue;
// 			min_energy = std::min(min_energy,get_e_intP(k,d,dp,l) + PR.get(i,j,d,dp));
// 		}
// 	}
// 	return min_energy;
// }

// energy_t pseudo_loop::calc_PMiloop(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
// 	if (!(i <= j && j < k-1 && k <= l)){
// 		return INF;
// 	}
// 	assert(!(i<=0 || l> n));

// 	if(!can_pair(j,k)) return INF;

// 	energy_t min_energy = INF;
// 	if ( i<j && k<l ) {
// 		min_energy = PM.get(i,j-1,k+1,l) + get_e_stP(j-1,k+1);
// 	}
// 	cand_pos_t max_d = std::max(i,j-MAXLOOP);
// 	for(cand_pos_t d= j-1; d>max_d; --d){
// 		cand_pos_t min_dp = std::min(l,k+MAXLOOP);
// 		for (cand_pos_t dp=k+1; dp <min_dp; ++dp) {
// 			if (!can_pair(d,dp)) continue;
// 			min_energy = std::min(min_energy,get_e_intP(d,j,k,dp) + PM.get(i,d,dp,l));
// 		}
// 	}
// 	return min_energy;
// }

// energy_t pseudo_loop::calc_POmiloop(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
// 	if (!(i <= j && j < k-1 && k <= l)){
// 		return INF;
// 	}
// 	assert(!(i<=0 || l> n));

// 	if(!can_pair(i,l)) return INF;

// 	energy_t min_energy = INF;
// 	if (i<j && k<l ) { 
// 		min_energy = POm.get(i+1,j,k,l-1) + get_e_stP(i,l);
// 	}
// 	cand_pos_t max_d = std::min(j,i+MAXLOOP);
// 	for(cand_pos_t d= i+1; d<max_d; ++d){
// 		cand_pos_t min_dp = std::max(l-MAXLOOP,k);
// 		for (cand_pos_t dp=l-1; dp >min_dp; --dp) {
// 			if (!can_pair(d,dp)) continue;
// 			min_energy = std::min(min_energy,get_e_intP(i,d,dp,l) + POm.get(d,j,dp,k));
// 		}
// 	}
// 	return min_energy;
// }

// energy_t pseudo_loop::calc_POsiloop(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
// 	if (!(i <= j && j < k-1 && k <= l)){
// 		return INF;
// 	}
// 	assert(!(i<=0 || l> n));

// 	if(!can_pair(i,l)) return INF;

// 	energy_t min_energy = INF;
// 	if (i<j && k<l ) { 
// 		min_energy = POs.get(i+1,j,k,l-1) + get_e_stP(i,l);
// 	}
// 	cand_pos_t max_d = std::min(j,i+MAXLOOP);
// 	for(cand_pos_t d= i+1; d<max_d; ++d){
// 		cand_pos_t min_dp = std::max(l-MAXLOOP,k);
// 		for (cand_pos_t dp=l-1; dp >min_dp; --dp) {
// 			if (!can_pair(d,dp)) continue;
// 			min_energy = std::min(min_energy,get_e_intP(i,d,dp,l) + POs.get(d,j,dp,k));
// 		}
// 	}
// 	return min_energy;
// }

// energy_t pseudo_loop::calc_PLreOiloop(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
// 	if (!(i <= j && j < k-1 && k <= l)){
// 		return INF;
// 	}
// 	assert(!(i<=0 || l> n));

// 	if(!can_pair(i,j)) return INF;

// 	energy_t min_energy = INF;
// 	if (i+TURN+2<j) { 
// 		min_energy = PLreO.get(i+1,j-1,k,l) + get_e_stP(i,j);
// 	}
// 	cand_pos_t max_d = std::min(j,i+MAXLOOP);
// 	for(cand_pos_t d= i+1; d<max_d; ++d){
// 		cand_pos_t min_dp = std::max(d+TURN,j-MAXLOOP);
// 		for(cand_pos_t dp = j-1; dp > min_dp; --dp){
// 			if (!can_pair(d,dp)) continue;
// 			min_energy = std::min(min_energy,get_e_intP(i,d,dp,j) + PLreO.get(d,dp,k,l));
// 		}
// 	}
// 	return min_energy;
// }

// energy_t pseudo_loop::calc_PLreRiloop(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
// 	if (!(i <= j && j < k-1 && k <= l)){
// 		return INF;
// 	}
// 	assert(!(i<=0 || l> n));

// 	if(!can_pair(i,j)) return INF;

// 	energy_t min_energy = INF;
// 	if (i+TURN+2<j) { 
// 		min_energy = PLreR.get(i+1,j-1,k,l) + get_e_stP(i,j);
// 	}
// 	cand_pos_t max_d = std::min(j,i+MAXLOOP);
// 	for(cand_pos_t d= i+1; d<max_d; ++d){
// 		cand_pos_t min_dp = std::max(d+TURN,j-MAXLOOP);
// 		for(cand_pos_t dp = j-1; dp > min_dp; --dp){
// 			if (!can_pair(d,dp)) continue;
// 			min_energy = std::min(min_energy,get_e_intP(i,d,dp,j) + PLreR.get(d,dp,k,l));
// 		}
// 	}
// 	return min_energy;
// }

// energy_t pseudo_loop::calc_PMreOiloop(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
// 	if (!(i <= j && j < k-1 && k <= l)){
// 		return INF;
// 	}
// 	assert(!(i<=0 || l> n));

// 	if(!can_pair(j,k)) return INF;

// 	energy_t min_energy = INF;
// 	if ( i<j && k<l ) {
// 		min_energy = PMreO.get(i,j-1,k+1,l) + get_e_stP(j-1,k+1);
// 	}
// 	cand_pos_t max_d = std::max(i,j-MAXLOOP);
// 	for(cand_pos_t d= j-1; d>max_d; --d){
// 		cand_pos_t min_dp = std::min(l,k+MAXLOOP);
// 		for (cand_pos_t dp=k+1; dp <min_dp; ++dp) {
// 			if (!can_pair(d,dp)) continue;
// 			min_energy = std::min(min_energy,get_e_intP(d,j,k,dp) + PMreO.get(i,d,dp,l));
// 		}
// 	}
// 	return min_energy;
// }

// energy_t pseudo_loop::calc_PMreRiloop(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
// 	if (!(i <= j && j < k-1 && k <= l)){
// 		return INF;
// 	}
// 	assert(!(i<=0 || l> n));

// 	if(!can_pair(j,k)) return INF;

// 	energy_t min_energy = INF;
// 	if ( i<j && k<l ) {
// 		min_energy = PMreR.get(i,j-1,k+1,l) + get_e_stP(j-1,k+1);
// 	}
// 	cand_pos_t max_d = std::max(i,j-MAXLOOP);
// 	for(cand_pos_t d= j-1; d>max_d; --d){
// 		cand_pos_t min_dp = std::min(l,k+MAXLOOP);
// 		for (cand_pos_t dp=k+1; dp <min_dp; ++dp) {
// 			if (!can_pair(d,dp)) continue;
// 			min_energy = std::min(min_energy,get_e_intP(d,j,k,dp) + PMreR.get(i,d,dp,l));
// 		}
// 	}
// 	return min_energy;
// }

// energy_t pseudo_loop::calc_PLMreRiloop(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
// 	if (!(i <= j && j < k-1 && k <= l)){
// 		return INF;
// 	}
// 	assert(!(i<=0 || l> n));

// 	if(!can_pair(i,j)) return INF;

// 	energy_t min_energy = INF;
// 	if (i+TURN+2<j) { 
// 		min_energy = PLMreR.get(i+1,j-1,k,l) + get_e_stP(i,j);
// 	}
// 	cand_pos_t max_d = std::min(j,i+MAXLOOP);
// 	for(cand_pos_t d= i+1; d<max_d; ++d){
// 		cand_pos_t min_dp = std::max(d+TURN,j-MAXLOOP);
// 		for(cand_pos_t dp = j-1; dp > min_dp; --dp){
// 			if (!can_pair(d,dp)) continue;
// 			min_energy = std::min(min_energy,get_e_intP(i,d,dp,j) + PLMreR.get(d,dp,k,l));
// 		}
// 	}
// 	return min_energy;
// }

// energy_t pseudo_loop::calc_PLMorOiloop(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
// 	if (!(i <= j && j < k-1 && k <= l)){
// 		return INF;
// 	}
// 	assert(!(i<=0 || l> n));

// 	if(!can_pair(i,j)) return INF;

// 	energy_t min_energy = INF;
// 	if (i+TURN+2<j) { 
// 		min_energy = PLMorO.get(i+1,j-1,k,l) + get_e_stP(i,j);
// 	}
// 	cand_pos_t max_d = std::min(j,i+MAXLOOP);
// 	for(cand_pos_t d= i+1; d<max_d; ++d){
// 		cand_pos_t min_dp = std::max(d+TURN,j-MAXLOOP);
// 		for(cand_pos_t dp = j-1; dp > min_dp; --dp){
// 			if (!can_pair(d,dp)) continue;
// 			min_energy = std::min(min_energy,get_e_intP(i,d,dp,j) + PLMorO.get(d,dp,k,l));
// 		}
// 	}
// 	return min_energy;
// }
/**
 * 
 * 
 * 
 */

// void pseudo_loop::compute_PfromL(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
// 	energy_t min_energy = INF;

// 	for(cand_pos_t d=i+1; d<=j; ++d){
// 		energy_t tmp = get_WP(i,d-1) + PfromLprime.get(d,j,k,l);
// 		min_energy = std::min(min_energy,tmp);
// 	}

// 	if (min_energy < INF/2){
// 		PfromL.set(i,j,k,l,min_energy);
// 	}
// }

// void pseudo_loop::compute_PfromR(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
// 	energy_t min_energy = INF;

// 	for(cand_pos_t d=k+1; d<=l; ++d){
// 		energy_t tmp = get_WP(k,d-1) + PfromRprime.get(i,j,d,l);
// 		min_energy = std::min(min_energy,tmp);
// 	}

// 	if (min_energy < INF/2){
// 		PfromR.set(i,j,k,l,min_energy);
// 	}
// }

// void pseudo_loop::compute_PfromM(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
// 	energy_t min_energy = INF;

// 	for(cand_pos_t d=i; d<j; ++d){
// 		energy_t tmp = PfromMprime.get(i,d,k,l)+get_WP(d+1,j);
// 		min_energy = std::min(min_energy,tmp);
// 	}

// 	if (min_energy < INF/2){
// 		PfromM.set(i,j,k,l,min_energy);
// 	}
// }

// void pseudo_loop::compute_PfromO(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
// 	energy_t min_energy = INF;

// 	for(cand_pos_t d=i+1; d<=j; ++d){
// 		energy_t tmp = get_WP(i,d-1) + PfromOprime.get(d,j,k,l);
// 		min_energy = std::min(min_energy,tmp);
// 	}

// 	if (min_energy < INF/2){
// 		PfromO.set(i,j,k,l,min_energy);
// 	}
// }

// void pseudo_loop::compute_PfromLreO(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
// 	energy_t min_energy = INF;

// 	for(cand_pos_t d=i+1; d<=j; ++d){
// 		energy_t tmp = get_WP(i,d-1) + PfromLreOprime.get(d,j,k,l);
// 		min_energy = std::min(min_energy,tmp);
// 	}

// 	if (min_energy < INF/2){
// 		PfromLreO.set(i,j,k,l,min_energy);
// 	}
// }

// void pseudo_loop::compute_PfromLreR(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
// 	energy_t min_energy = INF;

// 	for(cand_pos_t d=i+1; d<=j; ++d){
// 		energy_t tmp = get_WP(i,d-1) + PfromLreRprime.get(d,j,k,l);
// 		min_energy = std::min(min_energy,tmp);
// 	}

// 	if (min_energy < INF/2){
// 		PfromLreR.set(i,j,k,l,min_energy);
// 	}
// }

// void pseudo_loop::compute_PfromMreO(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
// 	energy_t min_energy = INF;

// 	for(cand_pos_t d=i+1; d<j; ++d){
// 		energy_t tmp = PfromMreOprime.get(i,d,k,l)+get_WP(d+1,j);
// 		min_energy = std::min(min_energy,tmp);
// 	}

// 	if (min_energy < INF/2){
// 		PfromMreO.set(i,j,k,l,min_energy);
// 	}
// }

// void pseudo_loop::compute_PfromMreR(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
// 	energy_t min_energy = INF;

// 	for(cand_pos_t d=i+1; d<j; ++d){
// 		energy_t tmp = PfromMreRprime.get(i,d,k,l)+get_WP(d+1,j);
// 		min_energy = std::min(min_energy,tmp);
// 	}

// 	if (min_energy < INF/2){
// 		PfromMreR.set(i,j,k,l,min_energy);
// 	}
// }

// void pseudo_loop::compute_PfromLMreR(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
// 	energy_t min_energy = INF;

// 	for(cand_pos_t d=i+1; d<=j; ++d){
// 		energy_t tmp = get_WP(i,d-1) + PfromLMreRprime.get(d,j,k,l);
// 		min_energy = std::min(min_energy,tmp);
// 	}

// 	if (min_energy < INF/2){
// 		PfromLMreR.set(i,j,k,l,min_energy);
// 	}
// }

// void pseudo_loop::compute_PfromLMorO(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
// 	energy_t min_energy = INF;

// 	for(cand_pos_t d=i+1; d<=j; ++d){
// 		energy_t tmp = get_WP(i,d-1) + PfromLMorOprime.get(d,j,k,l);
// 		min_energy = std::min(min_energy,tmp);
// 	}

// 	if (min_energy < INF/2){
// 		PfromLMorO.set(i,j,k,l,min_energy);
// 	}
// }

/**
 * 
 * 
 * 
 * 
 */
void pseudo_loop::compute_PfromLprime(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	energy_t min_energy = INF;

	for(cand_pos_t d=i; d<j; ++d){
		energy_t tmp=get_PfromLdoubleprime(i,d,k,l)+get_WP(d+1,j);
		min_energy = std::min(min_energy,tmp);
	}
	
	if (min_energy < INF/2){
		PfromLprime.set(i,j,k,l,min_energy);
	}
}

void pseudo_loop::compute_PfromRprime(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	energy_t min_energy = INF;

	for(cand_pos_t d=k; d<l; ++d){
		energy_t tmp=get_PfromRdoubleprime(i,j,k,d)+get_WP(d+1,l);
		min_energy = std::min(min_energy,tmp);
	}
	
	if (min_energy < INF/2){
		PfromRprime.set(i,j,k,l,min_energy);
	}
}

void pseudo_loop::compute_PfromMprime(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	energy_t min_energy = INF;

	for(cand_pos_t d=k+1; d<=l; ++d){
		energy_t tmp= get_WP(k,d-1) + get_PfromMdoubleprime(i,j,d,l);
		min_energy = std::min(min_energy,tmp);
	}
	
	if (min_energy < INF/2){
		PfromMprime.set(i,j,k,l,min_energy);
	}
}

void pseudo_loop::compute_PfromOprime(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	energy_t min_energy = INF;

	for(cand_pos_t d=k; d<l; ++d){
		energy_t tmp=get_PfromOdoubleprime(i,j,k,d)+get_WP(d+1,j);
		min_energy = std::min(min_energy,tmp);
	}
	
	if (min_energy < INF/2){
		PfromOprime.set(i,j,k,l,min_energy);
	}
}

void pseudo_loop::compute_PfromLreOprime(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	energy_t min_energy = INF;

	for(cand_pos_t d=i; d<j; ++d){
		energy_t tmp=get_PfromLreOdoubleprime(i,d,k,l)+get_WP(d+1,j);
		min_energy = std::min(min_energy,tmp);
	}
	
	if (min_energy < INF/2){
		PfromLreOprime.set(i,j,k,l,min_energy);
	}
}

void pseudo_loop::compute_PfromLreRprime(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	energy_t min_energy = INF;

	for(cand_pos_t d=i; d<j; ++d){
		energy_t tmp=get_PfromLreRdoubleprime(i,d,k,l)+get_WP(d+1,j);
		min_energy = std::min(min_energy,tmp);
	}
	
	if (min_energy < INF/2){
		PfromLreRprime.set(i,j,k,l,min_energy);
	}
}

void pseudo_loop::compute_PfromMreOprime(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	energy_t min_energy = INF;

	for(cand_pos_t d=k+1; d<l; ++d){
		energy_t tmp=get_PfromMreOdoubleprime(i,j,d,l)+get_WP(k,d-1);
		min_energy = std::min(min_energy,tmp);
	}
	
	if (min_energy < INF/2){
		PfromMreOprime.set(i,j,k,l,min_energy);
	}
}

void pseudo_loop::compute_PfromMreRprime(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	energy_t min_energy = INF;

	for(cand_pos_t d=k+1; d<l; ++d){
		energy_t tmp=get_PfromMreRdoubleprime(i,j,d,l)+get_WP(k,d-1);
		min_energy = std::min(min_energy,tmp);
	}
	
	if (min_energy < INF/2){
		PfromMreRprime.set(i,j,k,l,min_energy);
	}
}

void pseudo_loop::compute_PfromLMreRprime(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	energy_t min_energy = INF;

	for(cand_pos_t d=i; d<j; ++d){
		energy_t tmp=get_PfromLMreRdoubleprime(i,d,k,l)+get_WP(d+1,j);
		min_energy = std::min(min_energy,tmp);
	}
	
	if (min_energy < INF/2){
		PfromLMreRprime.set(i,j,k,l,min_energy);
	}
}

void pseudo_loop::compute_PfromLMorOprime(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	energy_t min_energy = INF;

	for(cand_pos_t d=i; d<j; ++d){
		energy_t tmp=get_PfromLMorOdoubleprime(i,d,k,l)+get_WP(d+1,j);
		min_energy = std::min(min_energy,tmp);
	}
	
	if (min_energy < INF/2){
		PfromLMorOprime.set(i,j,k,l,min_energy);
	}
}
/**
 * 
 * 
 */
// void pseudo_loop::compute_PLmloop00(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
// 	energy_t min_energy = INF,tmp=INF;

//     for(cand_pos_t d = i; d<=j; ++d){
//         if (d>i){
//             tmp=get_WB(i,d-1)+PLmloop10.get(d,j,k,l);
//             min_energy = std::min(min_energy,tmp);
//         }
//         if(d<j){
//             tmp=PLmloop01.get(i,d,k,l)+get_WB(d+1,j);
//             min_energy = std::min(min_energy,tmp);
//         }

//     }
// 	if (min_energy < INF/2){
// 		PLmloop00.set(i,j,k,l,min_energy);
// 	}
// }
// void pseudo_loop::compute_PMmloop00(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
// 	energy_t min_energy = INF,tmp=INF;

//     for(cand_pos_t d=i; d<j; ++d){
//         tmp=PMmloop10.get(i,d,k,l)+get_WB(d+1,j);
//         min_energy = std::min(min_energy,tmp);
//     }
//     for(cand_pos_t d=k+1; d<=l; ++d){
//         tmp=PMmloop01.get(i,j,d,l)+get_WB(k,d-1);
//         min_energy = std::min(min_energy,tmp);
//     }

// 	if (min_energy < INF/2){
// 		PMmloop00.set(i,j,k,l,min_energy);
// 	}
// }
// void pseudo_loop::compute_PRmloop00(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
// 	energy_t min_energy = INF,tmp=INF;

// 	for(cand_pos_t d=k; d<=l; ++d){
//         if(d>k){
//             tmp=get_WB(k,d-1)+PRmloop10.get(i,j,d,l);
//             min_energy = std::min(min_energy,tmp);
//         }
//         if (d<l){
//             tmp = PRmloop01.get(i,j,k,d)+get_WB(d+1,l);
//             min_energy = std::min(min_energy,tmp);
//         }
//     }
// 	if (min_energy < INF/2){
// 		PRmloop00.set(i,j,k,l,min_energy);
// 	}
// }
// void pseudo_loop::compute_POmmloop00(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
// 	energy_t min_energy = INF,tmp=INF;

//     for(cand_pos_t d=i+1; d<=j; ++d){
//         tmp = get_WB(i,d-1)+POmmloop10.get(d,j,k,l);
//         min_energy = std::min(min_energy,tmp);
//     }
//     for(cand_pos_t d=k; d<l; ++d){
//         tmp = POmmloop01.get(i,j,k,d)+get_WB(d+1,l);
//         min_energy = std::min(min_energy,tmp);
//     }

// 	if (min_energy < INF/2){
// 		POmmloop00.set(i,j,k,l,min_energy);
// 	}
// }
// void pseudo_loop::compute_POsmloop00(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
// 	energy_t min_energy = INF,tmp=INF;

//     for(cand_pos_t d=i+1; d<=j; ++d){
//         tmp = get_WB(i,d-1)+POsmloop10.get(d,j,k,l);
//         min_energy = std::min(min_energy,tmp);
//     }
//     for(cand_pos_t d=k; d<l; ++d){
//         tmp = POsmloop01.get(i,j,k,d)+get_WB(d+1,l);
//         min_energy = std::min(min_energy,tmp);
//     }

// 	if (min_energy < INF/2){
// 		POsmloop00.set(i,j,k,l,min_energy);
// 	}
// }
// void pseudo_loop::compute_PLreOmloop00(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
// 	energy_t min_energy = INF,tmp=INF;

//     for(cand_pos_t d = i; d<=j; ++d){
//         if (d>i){
//             tmp=get_WB(i,d-1)+PLreOmloop10.get(d,j,k,l);
//             min_energy = std::min(min_energy,tmp);
//         }
//         if(d<j){
//             tmp=PLreOmloop01.get(i,d,k,l)+get_WB(d+1,j);
//             min_energy = std::min(min_energy,tmp);
//         }

//     }
// 	if (min_energy < INF/2){
// 		PLreOmloop00.set(i,j,k,l,min_energy);
// 	}
// }
// void pseudo_loop::compute_PLreRmloop00(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
// 	energy_t min_energy = INF,tmp=INF;

//     for(cand_pos_t d = i; d<=j; ++d){
//         if (d>i){
//             tmp=get_WB(i,d-1)+PLreRmloop10.get(d,j,k,l);
//             min_energy = std::min(min_energy,tmp);
//         }
//         if(d<j){
//             tmp=PLreRmloop01.get(i,d,k,l)+get_WB(d+1,j);
//             min_energy = std::min(min_energy,tmp);
//         }

//     }
// 	if (min_energy < INF/2){
// 		PLreRmloop00.set(i,j,k,l,min_energy);
// 	}
// }
// void pseudo_loop::compute_PMreOmloop00(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
// 	energy_t min_energy = INF,tmp=INF;

//     for(cand_pos_t d=i; d<j; ++d){
//         tmp=PMreOmloop10.get(i,d,k,l)+get_WB(d+1,j);
//         min_energy = std::min(min_energy,tmp);
//     }
//     for(cand_pos_t d=k+1; d<=l; ++d){
//         tmp=PMreOmloop01.get(i,j,d,l)+get_WB(k,d-1);
//         min_energy = std::min(min_energy,tmp);
//     }

// 	if (min_energy < INF/2){
// 		PMreOmloop00.set(i,j,k,l,min_energy);
// 	}
// }
// void pseudo_loop::compute_PMreRmloop00(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
// 	energy_t min_energy = INF,tmp=INF;

//     for(cand_pos_t d=i; d<j; ++d){
//         tmp=PMreRmloop10.get(i,d,k,l)+get_WB(d+1,j);
//         min_energy = std::min(min_energy,tmp);
//     }
//     for(cand_pos_t d=k+1; d<=l; ++d){
//         tmp=PMreRmloop01.get(i,j,d,l)+get_WB(k,d-1);
//         min_energy = std::min(min_energy,tmp);
//     }

// 	if (min_energy < INF/2){
// 		PMreRmloop00.set(i,j,k,l,min_energy);
// 	}
// }
// void pseudo_loop::compute_PLMreRmloop00(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
// 	energy_t min_energy = INF,tmp=INF;

//     for(cand_pos_t d = i; d<=j; ++d){
//         if (d>i){
//             tmp=get_WB(i,d-1)+PLMreRmloop10.get(d,j,k,l);
//             min_energy = std::min(min_energy,tmp);
//         }
//         if(d<j){
//             tmp=PLMreRmloop01.get(i,d,k,l)+get_WB(d+1,j);
//             min_energy = std::min(min_energy,tmp);
//         }

//     }
// 	if (min_energy < INF/2){
// 		PLMreRmloop00.set(i,j,k,l,min_energy);
// 	}
// }
// void pseudo_loop::compute_PLMorOmloop00(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
// 	energy_t min_energy = INF,tmp=INF;

//     for(cand_pos_t d = i; d<=j; ++d){
//         if (d>i){
//             tmp=get_WB(i,d-1)+PLMorOmloop10.get(d,j,k,l);
//             min_energy = std::min(min_energy,tmp);
//         }
//         if(d<j){
//             tmp=PLMorOmloop01.get(i,d,k,l)+get_WB(d+1,j);
//             min_energy = std::min(min_energy,tmp);
//         }

//     }
// 	if (min_energy < INF/2){
// 		PLMorOmloop00.set(i,j,k,l,min_energy);
// 	}
// }
/**
 * 
 * 
 * 
 */
// void pseudo_loop::compute_PLmloop10(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
// 	energy_t min_energy = INF,tmp=INF;

// 	for(cand_pos_t d = i+1; d <= j; ++d){
//         tmp = PL.get(i,d,k,l) + get_WB(d+1,j);
//         min_energy = std::min(min_energy,tmp);
//     }
// 	if (min_energy < INF/2){
// 		PLmloop10.set(i,j,k,l,min_energy);
// 	}

// }
// void pseudo_loop::compute_PMmloop10(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
// 	energy_t min_energy = INF,tmp=INF;

//     for(cand_pos_t d = i+1; d <=j; ++d){
//         tmp = PM.get(i,d,k,d) + get_WB(d+1,j);
//         min_energy = std::min(min_energy,tmp);
//     }

// 	if (min_energy < INF/2){
// 		PMmloop10.set(i,j,k,l,min_energy);
// 	}
// }
// void pseudo_loop::compute_PRmloop10(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
// 	energy_t min_energy = INF,tmp=INF;

//     for(cand_pos_t d = k+1; d <= l; ++d){
//         tmp = PR.get(i,j,d,l) + get_WB(d+1,l);
//         min_energy = std::min(min_energy,tmp);
//     }
// 	if (min_energy < INF/2){
// 		PRmloop10.set(i,j,k,l,min_energy);
// 	}

// }
// void pseudo_loop::compute_POmmloop10(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
// 	energy_t min_energy = INF,tmp=INF;

// 	for(cand_pos_t d=i+1; d<=j;++d){
//         tmp=POm.get(d,j,k,l) + get_WB(d+1,l);
//         min_energy = std::min(min_energy,tmp);
//     }
// 	if (min_energy < INF/2){
// 		POmmloop10.set(i,j,k,l,min_energy);
// 	}
// }
// void pseudo_loop::compute_POsmloop10(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
// 	energy_t min_energy = INF,tmp=INF;

// 	for(cand_pos_t d=i+1; d<=j;++d){
//         tmp= POs.get(d,j,k,l) + get_WB(d+1,l);
//         min_energy = std::min(min_energy,tmp);
//     }
// 	if (min_energy < INF/2){
// 		POsmloop10.set(i,j,k,l,min_energy);
// 	}
// }
// void pseudo_loop::compute_PLreOmloop10(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
// 	energy_t min_energy = INF,tmp=INF;

// 	for(cand_pos_t d = i+1; d <= j; ++d){
//         tmp = PLreO.get(i,d,k,l) + get_WB(d+1,j);
//         min_energy = std::min(min_energy,tmp);
//     }
// 	if (min_energy < INF/2){
// 		PLreOmloop10.set(i,j,k,l,min_energy);
// 	}

// }
// void pseudo_loop::compute_PLreRmloop10(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
// 	energy_t min_energy = INF,tmp=INF;

// 	for(cand_pos_t d = i+1; d <= j; ++d){
//         tmp = PLreR.get(i,d,k,l) + get_WB(d+1,j);
//         min_energy = std::min(min_energy,tmp);
//     }
// 	if (min_energy < INF/2){
// 		PLreRmloop10.set(i,j,k,l,min_energy);
// 	}

// }
// void pseudo_loop::compute_PMreOmloop10(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
// 	energy_t min_energy = INF,tmp=INF;

//     for(cand_pos_t d = i+1; d <=j; ++d){
//         tmp = PMreO.get(i,d,k,l) + get_WB(d+1,j);
//         min_energy = std::min(min_energy,tmp);
//     }
// 	if (min_energy < INF/2){
// 		PMreOmloop10.set(i,j,k,l,min_energy);
// 	}
// }
// void pseudo_loop::compute_PMreRmloop10(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
// 	energy_t min_energy = INF,tmp=INF;

//     for(cand_pos_t d = i+1; d <=j; ++d){
//         tmp = PMreR.get(i,d,k,l) + get_WB(d+1,j);
//         min_energy = std::min(min_energy,tmp);
//     }
// 	if (min_energy < INF/2){
// 		PMreRmloop10.set(i,j,k,l,min_energy);
// 	}
// }
// void pseudo_loop::compute_PLMreRmloop10(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
// 	energy_t min_energy = INF,tmp=INF;

// 	for(cand_pos_t d = i+1; d <= j; ++d){
//         tmp = PLMreR.get(i,d,k,l) + get_WB(d+1,j);
//         min_energy = std::min(min_energy,tmp);
//     }
// 	if (min_energy < INF/2){
// 		PLMreRmloop10.set(i,j,k,l,min_energy);
// 	}

// }
// void pseudo_loop::compute_PLMorOmloop10(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
// 	energy_t min_energy = INF,tmp=INF;

// 	for(cand_pos_t d = i+1; d <= j; ++d){
//         tmp = PLMorO.get(i,d,k,l) + get_WB(d+1,j);
//         min_energy = std::min(min_energy,tmp);
//     }
// 	if (min_energy < INF/2){
// 		PLMorOmloop10.set(i,j,k,l,min_energy);
// 	}

// }
/**
 * 
 * 
 * 
 */
// void pseudo_loop::compute_PLmloop01(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
// 	energy_t min_energy = INF,tmp=INF;

// 	for(cand_pos_t d = i; d < j; ++d){
//         tmp = cp_penalty*(d-i) + PL.get(d,j,k,l);
//         min_energy = std::min(min_energy,tmp);
//     }
// 	if (min_energy < INF/2){
// 		PLmloop01.set(i,j,k,l,min_energy);
// 	}

// }
// void pseudo_loop::compute_PMmloop01(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
// 	energy_t min_energy = INF,tmp=INF;

// 	for(cand_pos_t d = k; d < l; ++d){
//         tmp = cp_penalty*(d-k) + PM.get(i,j,k,d);
//         min_energy = std::min(min_energy,tmp);
//     }

// 	if (min_energy < INF/2){
// 		PMmloop01.set(i,j,k,l,min_energy);
// 	}
// }
// void pseudo_loop::compute_PRmloop01(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
// 	energy_t min_energy = INF,tmp=INF;

//     for(cand_pos_t d = k; d < l; d++){
//         tmp = cp_penalty*(d-k) + PR.get(i,j,d,l);
//         min_energy = std::min(min_energy,tmp);
//     }
// 	if (min_energy < INF/2){
// 		PRmloop01.set(i,j,k,l,min_energy);
// 	}

// }
// void pseudo_loop::compute_POmmloop01(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
// 	energy_t min_energy = INF;

// 	for(cand_pos_t d = k; d < l; ++d){
//         energy_t tmp = cp_penalty*(d-i) + POm.get(d,j,k,l);
//         min_energy = std::min(min_energy,tmp);
//     }

// 	if (min_energy < INF/2){
// 		POmmloop01.set(i,j,k,l,min_energy);
// 	}

// }
// void pseudo_loop::compute_POsmloop01(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
// 	energy_t min_energy = INF;

// 	for(cand_pos_t d = k; d < l; ++d){
//         energy_t tmp = cp_penalty*(d-i) + POs.get(d,j,k,l);
//         min_energy = std::min(min_energy,tmp);
//     }

// 	if (min_energy < INF/2){
// 		POsmloop01.set(i,j,k,l,min_energy);
// 	}

// }
// void pseudo_loop::compute_PLreOmloop01(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
// 	energy_t min_energy = INF,tmp=INF;

// 	for(cand_pos_t d = i; d < j; ++d){
//         tmp = cp_penalty*(d-i) + PLreO.get(d,j,k,l);
//         min_energy = std::min(min_energy,tmp);
//     }
// 	if (min_energy < INF/2){
// 		PLreOmloop01.set(i,j,k,l,min_energy);
// 	}

// }
// void pseudo_loop::compute_PLreRmloop01(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
// 	energy_t min_energy = INF,tmp=INF;

// 	for(cand_pos_t d = i; d < j; ++d){
//         tmp = cp_penalty*(d-i) + PLreR.get(d,j,k,l);
//         min_energy = std::min(min_energy,tmp);
//     }
// 	if (min_energy < INF/2){
// 		PLreRmloop01.set(i,j,k,l,min_energy);
// 	}

// }
// void pseudo_loop::compute_PMreOmloop01(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
// 	energy_t min_energy = INF,tmp=INF;

// 	for(cand_pos_t d = k; d < l; ++d){
//         tmp = cp_penalty*(d-k) + PMreO.get(i,j,k,d);
//         min_energy = std::min(min_energy,tmp);
//     }

// 	if (min_energy < INF/2){
// 		PMreOmloop01.set(i,j,k,l,min_energy);
// 	}
// }
// void pseudo_loop::compute_PMreRmloop01(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
// 	energy_t min_energy = INF,tmp=INF;

// 	for(cand_pos_t d = k; d < l; ++d){
//         tmp = cp_penalty*(d-k) + PMreR.get(i,j,k,d);
//         min_energy = std::min(min_energy,tmp);
//     }

// 	if (min_energy < INF/2){
// 		PMreRmloop01.set(i,j,k,l,min_energy);
// 	}
// }
// void pseudo_loop::compute_PLMreRmloop01(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
// 	energy_t min_energy = INF,tmp=INF;

// 	for(cand_pos_t d = i; d < j; ++d){
//         tmp = cp_penalty*(d-i) + PLMreR.get(d,j,k,l);
//         min_energy = std::min(min_energy,tmp);
//     }
// 	if (min_energy < INF/2){
// 		PLMreRmloop01.set(i,j,k,l,min_energy);
// 	}

// }
// void pseudo_loop::compute_PLMorOmloop01(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
// 	energy_t min_energy = INF,tmp=INF;

// 	for(cand_pos_t d = i; d < j; ++d){
//         tmp = cp_penalty*(d-i) + PLMorO.get(d,j,k,l);
//         min_energy = std::min(min_energy,tmp);
//     }
// 	if (min_energy < INF/2){
// 		PLMorOmloop01.set(i,j,k,l,min_energy);
// 	}

// }