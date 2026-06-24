#include "pseudo_loop.hh"
#include "h_externs.hh"
#include <stdio.h>
#include <string>
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <algorithm>
#include <cassert>

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

	TriangleMatrix::new_index(index,n+1);
	Matrix4D::construct_index(index3D,n);
	init_can_pair();

    WBP.init(n+1,index);
	WPP.init(n+1,index);
	P.init(n+1,index);
	
	// 4D matrix initialization
	// PK.init(n,index3D);
	// PO.init(n,index3D);

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

			compute_PfromLprime(i,j,k,l);
			compute_PfromRprime(i,j,k,l);
			compute_PfromMprime(i,j,k,l);
			compute_PfromOprime(i,j,k,l);
			compute_PfromLreOprime(i,j,k,l);
			compute_PfromLreRprime(i,j,k,l);
			compute_PfromMreOprime(i,j,k,l);
			compute_PfromMreRprime(i,j,k,l);
			compute_PfromLMreRprime(i,j,k,l);
			compute_PfromLMorOprime(i,j,k,l);

			compute_PK1X(i,j,k,l,MType::Om);
			compute_PK1X(i,j,k,l,MType::Os);
			compute_PK1X(i,j,k,l,MType::LreO);
			compute_PK1X(i,j,k,l,MType::LreR);
			compute_PK1X(i,j,k,l,MType::MreO);
			compute_PK1X(i,j,k,l,MType::MreR);
			compute_PK1X(i,j,k,l,MType::R);
			compute_PK1X(i,j,k,l,MType::LMreR);
			compute_PK1X(i,j,k,l,MType::LMorO);

			compute_PK2X(i,j,k,l,MType::Om);
			compute_PK2X(i,j,k,l,MType::Os);
			compute_PK2X(i,j,k,l,MType::LreO);
			compute_PK2X(i,j,k,l,MType::LreR);
			compute_PK2X(i,j,k,l,MType::MreO);
			compute_PK2X(i,j,k,l,MType::MreR);
			compute_PK2X(i,j,k,l,MType::R);
			if(i==1 && j==4 && k==14 && l==17) std::cout << POs.get(i,j,k,l) << std::endl;
			if(i==5 && j==12 && k==18 && l==25) std::cout << POs.get(i,j,k,l) << std::endl;
		}
	}

}

void pseudo_loop::compute_WBP(cand_pos_t i, cand_pos_t l){
	energy_t min_energy= INF, b1 = INF, b2=INF, b3 = INF,tmp =INF;

	for(cand_pos_t d=i; d< l; ++d){
		tmp = get_WB(i,d-1) + V->get_energy(d,l) + bp_penalty + PPS_penalty;
		b1 = std::min(b1,tmp);
		tmp = get_WB(i,d-1) + P.get(d,l) + PSM_penalty + PPS_penalty;
		b2 = std::min(b2,tmp);
	}
	b3 = WBP.get(i,l-1)+cp_penalty;
	min_energy = std::min({b1,b2,b3});
	if (min_energy < INF/2){
		WBP.set(i,l) = min_energy;
	}
}

void pseudo_loop::compute_WPP(cand_pos_t i, cand_pos_t l){
	energy_t min_energy = INF, b1 = INF, b2=INF, b3 =INF, tmp = INF;

	for(cand_pos_t d=i; d<l; ++d){
		tmp = get_WP(i,d-1) + V->get_energy(d,l) + gamma2(l,d) + PPS_penalty;
		b1 = std::min(b1,tmp);
		tmp = get_WP(i,d-1) + P.get(d,l) + PSP_penalty + PPS_penalty;
		b2 = std::min(b2,tmp);
	}
	b3 = WPP.get(i,l-1)+PUP_penalty;
	min_energy = std::min({b1,b2,b3});
	if (min_energy < INF/2){
		WPP.set(i,l) = min_energy;
	}
}

void pseudo_loop::compute_P(cand_pos_t i, cand_pos_t l){
	energy_t min_energy = INF;
	energy_t b1=INF,b2=INF,b3=INF,b4=INF,b5=INF,b6=INF,b7=INF,b8=INF,b9=INF,b10=INF,b11=INF,b12=INF,b13=INF,b14=INF,b15=INF,b16=INF,b17=INF,b18=INF,b19=INF,b20=INF,b21=INF;
	for(cand_pos_t j=i; j<l; ++j){
		for (cand_pos_t d=j+1; d<l; ++d){
			for (cand_pos_t k=d+1; k<l; ++k){
				b1 = PK1Om.get(i,j,d+1,k) + PK2Om.get(j+1,d,k+1,l);
				b2 = PK1Om.get(i,j,d+1,k) + PK2Os.get(j+1,d,k+1,l);
				b3 = PK1Om.get(i,j,d+1,k) + PK2LreO.get(j+1,d,k+1,l);
				b4 = PK1Om.get(i,j,d+1,k) + PK2MreO.get(j+1,d,k+1,l);
				b5 = PK1Om.get(i,j,d+1,k) + PK2LreR.get(j+1,d,k+1,l);
				b6 = PK1Om.get(i,j,d+1,k) + PK2MreR.get(j+1,d,k+1,l);
				b7 = PK1Om.get(i,j,d+1,k) + PK2R.get(j+1,d,k+1,l);

				b8 = PK1Os.get(i,j,d+1,k) + PK2Om.get(j+1,d,k+1,l);
				b9 = PK1Os.get(i,j,d+1,k) + PK2Os.get(j+1,d,k+1,l);
				b10 = PK1Os.get(i,j,d+1,k) + PK2MreO.get(j+1,d,k+1,l);
				b11 = PK1Os.get(i,j,d+1,k) + PK2MreR.get(j+1,d,k+1,l);
				b12 = PK1Os.get(i,j,d+1,k) + PK2R.get(j+1,d,k+1,l);

				b13 = PK1LreO.get(i,j,d+1,k) + PK2Om.get(j+1,d,k+1,l);
				b14 = PK1LreO.get(i,j,d+1,k) + PK2LreO.get(j+1,d,k+1,l);
				b15 = PK1LreO.get(i,j,d+1,k) + PK2MreO.get(j+1,d,k+1,l);

				b16 = PK1LMreR.get(i,j,d+1,k) + PK2Om.get(j+1,d,k+1,l);
				b17 = PK1LMreR.get(i,j,d+1,k) + PK2LreO.get(j+1,d,k+1,l);
				b18 = PK1LMreR.get(i,j,d+1,k) + PK2MreO.get(j+1,d,k+1,l);

				b19 = PK1LMorO.get(i,j,d+1,k) + PK2LreR.get(j+1,d,k+1,l);
				b20 = PK1LMorO.get(i,j,d+1,k) + PK2MreO.get(j+1,d,k+1,l);
				b21 = PK1LMorO.get(i,j,d+1,k) + PK2MreO.get(j+1,d,k+1,l);
				min_energy = std::min({min_energy,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,b15,b16,b17,b18,b19,b20,b21});
			}
		}
	}
	if (min_energy < INF/2){
		P.set(i,l) = min_energy;
	}
}

template<class Penalty> energy_t pseudo_loop::penalty(const Index4D &x, Penalty p, MType type) {
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
/**
 * This always chops from the interior with every single recurrence which means the code is the same and can be made generic
 */
void pseudo_loop::compute_PK1X(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l, MType type){
	energy_t min_energy = INF;
	Matrix4D &PX = PX_by_mtype(type);

	for(cand_pos_t d = i+1;d<=j;++d){
		min_energy = std::min(min_energy,PX.get(i,d,k,l) + get_WP(d+1,j) + PB_penalty);
	}
	Matrix4D &PK1X = PK1X_by_mtype(type);
	PK1X.set(i,j,k,l,min_energy);
}
/**
 * Same as PK1
 */
void pseudo_loop::compute_PK2X(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l, MType type){
	energy_t min_energy = INF;

	Matrix4D &PK1X = PK1X_by_mtype(type);
	for(cand_pos_t d = k;d<l;++d){
		min_energy = std::min(min_energy,PK1X.get(i,j,d,l) + get_WP(k,d-1));
	}
	Matrix4D &PK2X = PK2X_by_mtype(type);
	PK2X.set(i,j,k,l,min_energy);
}
/**
 * This was a case where the code was almost the same for all recurrences.
 * Each time, you would take the respective recurrence, shrink the related indices
 * by one and then add the penalties. This could be solve by the Index4D shrink function
 * which shrinks based on the Mtype. And so 10 functions become one
 */
energy_t pseudo_loop::calc_PXmloop(const Index4D &x, MType type){
	if(impossible_case(x)) return INF;

	Index4D xp(x);
	xp.shrink(type);
	Matrix4D &PXmloop00 = PXmloop00_by_mtype(type);
	return PXmloop00.get(xp.i(),xp.j(),xp.k(),xp.l())+ ap_penalty + bp_penalty;
}

energy_t pseudo_loop::calc_PXiloop(const Index4D &x, MType type){
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
    __builtin_unreachable();
}

void pseudo_loop::compute_PfromX(const Index4D &x, MType type){
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
    __builtin_unreachable();
}

// void pseudo_loop::compute_PfromXprime(const Index4D &x, MType type){
//     switch(type) {
//     case MType::L: return compute_PfromLprime(x,type);
//     case MType::M: return compute_PfromMprime(x,type);
//     case MType::R: return compute_PfromRprime(x,type);
//     case MType::Om: return compute_PfromOprime(x,type);
// 	case MType::LreO: return compute_PfromLprime(x,type);
// 	case MType::LreR: return compute_PfromLprime(x,type);
// 	case MType::MreO: return compute_PfromMprime(x,type);
// 	case MType::MreR: return compute_PfromMprime(x,type);
// 	case MType::LMreR: return compute_PfromLprime(x,type);
// 	case MType::LMorO: return compute_PfromLprime(x,type);
//     }
//     __builtin_unreachable();
// }

void pseudo_loop::compute_PXmloop00(const Index4D &x, MType type){
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
    __builtin_unreachable();
}
void pseudo_loop::compute_PXmloop01(const Index4D &x, MType type){
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
    __builtin_unreachable();
}
void pseudo_loop::compute_PXmloop10(const Index4D &x, MType type){
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
    __builtin_unreachable();
}
/**
 * As this just calls the other functions, we can reduce it to just a PX.
 * ptype closing can be done because lend and rend give the correct indices
 * for each recurrence. The functions coming after just use their respective
 * PX functions
 */
void pseudo_loop::compute_PX(const Index4D &x, MType type){
    energy_t min_energy = INF,b1=INF,b2=INF,b3=INF;
	Matrix4D &PX = PX_by_mtype(type);

	const int ptype_closing = pair[S_[x.lend(type)]][S_[x.rend(type)]];

	if (ptype_closing>0){
		b1 = calc_PXiloop(x, type);
		b2 = calc_PXmloop(x,type);

		// Hosna, July 11, 2014
		// To avoid addition of close base pairs we check for the following here
		if(type == MType::Os){
			if(x.i()==x.j() && x.k()==x.l()){
				b3=gamma2(x.l(),x.i());
			}
		} else if (x.difference(type)>TURN){
			Index4D xp(x);
			xp.shrink(type);
			Matrix4D &PfromX = PfromX_by_mtype(type);
			b3 = PfromX.get(xp) + penalty(xp, gamma2, type);
		}
	}
	min_energy = std::min({b1,b2,b3});
	if (min_energy < INF/2){
		PX.setI(x, min_energy);
	}
}
energy_t pseudo_loop::get_WB(cand_pos_t i, cand_pos_t j){
	if (i<=0 || j<=0 || i>n || j>n){
		return INF;
	}
	if (i>j) return 0;
	return (std::min(cp_penalty*(j-i+1),WBP.get(i,j)));
}

energy_t pseudo_loop::get_WP(cand_pos_t i, cand_pos_t j){
	if (i<=0 || j<=0 || i>n || j>n){
		return INF;
	}
	if (i>j) return 0; // needed as this will happen 
	return (std::min(PUP_penalty*(j-i+1),WPP.get(i,j)));
}

energy_t pseudo_loop::get_PfromLdoubleprime(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	if (!(i <= j && j < k-1 && k <= l)){
		return INF;
	}
	assert(!(i<=0 || l> n));

	energy_t b1 = PR.get(i,j,k,l) + gamma2(l,k) + PB_penalty;
	energy_t b2 = PM.get(i,j,k,l) + gamma2(j,k) + PB_penalty;
	energy_t b3 = get_PO(i,j,k,l) + gamma2(l,i) + PB_penalty;

	return std::min({b1,b2,b3});
}

energy_t pseudo_loop::get_PfromRdoubleprime(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	if (!(i <= j && j < k-1 && k <= l)){
		return INF;
	}
	assert(!(i<=0 || l> n));

	energy_t b1 = PM.get(i,j,k,l) + gamma2(j,k) + PB_penalty;
	energy_t b2 = get_PO(i,j,k,l) + gamma2(l,i) + PB_penalty;

	return std::min({b1,b2});
}

energy_t pseudo_loop::get_PfromMdoubleprime(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	if (!(i <= j && j < k-1 && k <= l)){
		return INF;
	}
	assert(!(i<=0 || l> n));

	energy_t b1 = PL.get(i,j,k,l) + gamma2(j,i) + PB_penalty;
	energy_t b2 = PR.get(i,j,k,l) + gamma2(l,k) + PB_penalty;
	energy_t b3 = POm.get(i,j,k,l) + gamma2(l,i) + PB_penalty;

	return std::min({b1,b2,b3});
}

energy_t pseudo_loop::get_PfromOdoubleprime(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	if (!(i <= j && j < k-1 && k <= l)){
		return INF;
	}
	assert(!(i<=0 || l> n));

	energy_t b1 = PL.get(i,j,k,l) + gamma2(j,i) + PB_penalty;
	energy_t b2 = PR.get(i,j,k,l) + gamma2(l,k) + PB_penalty;

	return std::min({b1,b2});
}

energy_t pseudo_loop::get_PfromLreOdoubleprime(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	if (!(i <= j && j < k-1 && k <= l)){
		return INF;
	}
	assert(!(i<=0 || l> n));

	energy_t b1 = PMreO.get(i,j,k,l) + gamma2(j,k) + PB_penalty;
	energy_t b2 = get_PO(i,j,k,l) + gamma2(l,i) + PB_penalty;

	return std::min({b1,b2});
}

energy_t pseudo_loop::get_PfromLreRdoubleprime(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	if (!(i <= j && j < k-1 && k <= l)){
		return INF;
	}
	assert(!(i<=0 || l> n));

	energy_t b1 = PMreR.get(i,j,k,l) + gamma2(j,k) + PB_penalty;
	energy_t b2 = PR.get(i,j,k,l) + gamma2(l,k) + PB_penalty;

	return std::min({b1,b2});
}

energy_t pseudo_loop::get_PfromMreOdoubleprime(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	if (!(i <= j && j < k-1 && k <= l)){
		return INF;
	}
	assert(!(i<=0 || l> n));

	energy_t b1 = PLreO.get(i,j,k,l) + gamma2(j,i) + PB_penalty;
	energy_t b2 = get_PO(i,j,k,l) + gamma2(l,i) + PB_penalty;

	return std::min({b1,b2});
}

energy_t pseudo_loop::get_PfromMreRdoubleprime(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	if (!(i <= j && j < k-1 && k <= l)){
		return INF;
	}
	assert(!(i<=0 || l> n));

	energy_t b1 = PLreR.get(i,j,k,l) + gamma2(j,i) + PB_penalty;
	energy_t b2 = PR.get(i,j,k,l) + gamma2(l,k) + PB_penalty;

	return std::min({b1,b2});
}

energy_t pseudo_loop::get_PfromLMreRdoubleprime(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	if (!(i <= j && j < k-1 && k <= l)){
		return INF;
	}
	assert(!(i<=0 || l> n));

	energy_t b1 = PMreR.get(i,j,k,l) + gamma2(j,k) + PB_penalty;

	return b1;
}

energy_t pseudo_loop::get_PfromLMorOdoubleprime(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	if (!(i <= j && j < k-1 && k <= l)){
		return INF;
	}
	assert(!(i<=0 || l> n));

	energy_t b1 = PM.get(i,j,k,l) + gamma2(j,k) + PB_penalty;
	energy_t b2 = get_PO(i,j,k,l) + gamma2(l,i) + PB_penalty;

	return std::min(b1,b2);
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

void pseudo_loop::backtrack(minimum_fold *f, seq_interval *cur_interval){

	// // this->structure = structure;
	// this->f = f;
	// switch (cur_interval->type)
	// {
	// 	case P_P:
	// 	{
	// 		cand_pos_t i = cur_interval->i;
	// 		cand_pos_t l = cur_interval->j;
	// 		if(debug) printf("At P_P with i=%d,l=%d\n",i,l);

	// 		if (i >= l){
	// 			std::cerr << "border case: This should not have happened!, P_P" << std::endl;
	// 			exit(EXIT_FAILURE);
	// 		}
	// 		energy_t min_energy = INF,b1=INF;

	// 		cand_pos_t best_d=0, best_j=0,best_k=0;
	// 		for(cand_pos_t j=i; j< l; j++){
	// 			for (cand_pos_t d=j+1; d<l; d++){
	// 				for (cand_pos_t k=d+1; k<l; k++){
	// 					b1 = PK.get(i,j,d+1,k) + PK.get(j+1,d,k+1,l);
	// 					if(b1 < min_energy){
	// 						min_energy = b1;
	// 						best_d = d;
	// 						best_j = j;
	// 						best_k= k;
	// 					}
	// 				}
	// 			}
	// 		}

	// 		insert_node(i,best_k,best_j,best_d+1,P_PK);
	// 		insert_node(best_j+1,l,best_d,best_k+1,P_PK);
	// 	}
	// 	break;

	// 	case P_PK:
	// 	{
	// 		cand_pos_t i = cur_interval->i;
	// 		cand_pos_t l = cur_interval->j;
	// 		cand_pos_t j = cur_interval->k;
	// 		cand_pos_t k = cur_interval->l;
	// 		if(debug) printf("At P_PK with i=%d,l=%d,j=%d,k=%d\n",i,l,j,k);
		
	// 		if (!(i <= j && j < k-1 && k <= l)){
	// 			std::cerr << "border cases: This should not have happened!, P_PK" << std::endl;
	// 			exit(EXIT_FAILURE);
	// 		}

	// 		if (i<=0 ||j<=0 || k<=0 || l<=0 || i>n || j>n || k> n || l> n){
	// 			std::cerr << "impossible cases: This should not have happened!, P_PK" << std::endl;
	// 			exit(EXIT_FAILURE);
	// 		}
	// 		energy_t min_energy = INF,temp=INF;
	// 		cand_pos_t best_row = -1,best_d=-1;

	// 		// branch 1
	// 		for(cand_pos_t  d=i+1; d< j; ++d){
	// 			temp = PK.get(i,d,k,l) + get_WP(d+1,j);
	// 			if (temp < min_energy){
	// 				min_energy=temp;
	// 				best_row = 1;
	// 				best_d=d;
	// 			}
	// 		}

	// 		// branch 2
	// 		for(cand_pos_t d=k+1; d< l; ++d){
	// 			temp = PK.get(i,j,d,l) + get_WP(k,d-1);
	// 			if (temp < min_energy){
	// 				min_energy=temp;
	// 				best_row = 2;
	// 				best_d = d;
	// 			}
	// 		}

	// 		// branch 3
	// 		temp = PL.get(i,j,k,l) + gamma2(j,i)+PB_penalty;
	// 		if(temp < min_energy){
	// 			min_energy = temp;
	// 			best_row = 3;
	// 			best_d = -1;
	// 		}

	// 		//branch 4
	// 		temp = PM.get(i,j,k,l) + gamma2(j,k)+PB_penalty;
	// 		if(temp < min_energy){
	// 			min_energy = temp;
	// 			best_row = 4;
	// 			best_d = -1;
	// 		}

	// 		// branch 5
	// 		temp = PR.get(i,j,k,l) + gamma2(l,k)+PB_penalty;
	// 		if(temp < min_energy){
	// 			min_energy = temp;
	// 			best_row = 5;
	// 			best_d = -1;
	// 		}

	// 		// branch 6
	// 		temp = PO.get(i,j,k,l) + gamma2(l,i)+PB_penalty;
	// 		if(temp < min_energy){
	// 			min_energy = temp;
	// 			best_row = 6;
	// 			best_d = -1;
	// 		}

	// 		switch (best_row){
	// 			case 1:
	// 				if (best_d > -1){
	// 					insert_node(i,l,best_d,k,P_PK);
	// 					insert_node(best_d+1,j,P_WP);
	// 				}
	// 				break;
	// 			case 2:
	// 				if (best_d > -1){
	// 					insert_node(i,l,j,best_d,P_PK);
	// 					insert_node(k,best_d-1,P_WP);
	// 				}
	// 				break;
	// 			case 3:
	// 				insert_node(i,l,j,k,P_PL);
	// 				break;
	// 			case 4:
	// 				insert_node(i,l,j,k,P_PM);
	// 				break;
	// 			case 5:
	// 				insert_node(i,l,j,k,P_PR);
	// 				break;
	// 			case 6:
	// 				insert_node(i,l,j,k,P_PO);
	// 				break;
	// 		}
	// 	}
	// 	break;

	// 	case P_PL:
	// 	{
	// 		cand_pos_t i = cur_interval->i;
	// 		cand_pos_t l = cur_interval->j;
	// 		cand_pos_t j = cur_interval->k;
	// 		cand_pos_t k = cur_interval->l;
	// 		if(debug) printf("At P_PL with i=%d,l=%d,j=%d,k=%d\n",i,l,j,k);

	// 		if (!(i <= j && j < k-1 && k <= l)){
	// 			std::cerr << "border cases: This should not have happened!, P_PL" << std::endl;
	// 			exit(EXIT_FAILURE);
	// 		}

	// 		if (i<=0 ||j<=0 || k<=0 || l<=0 || i>n || j>n || k> n || l> n){
	// 			std::cerr << "impossible cases: This should not have happened!, P_PL" << std::endl;
	// 			exit(EXIT_FAILURE);
	// 		}

	// 		energy_t min_energy = INF,temp=INF;
	// 		cand_pos_t best_row = -1;

	// 		if (pair[S_[i]][S_[j]]> 0){

	// 			//branch 1
	// 			temp = get_PLiloop(i,j,k,l);
	// 			if(temp < min_energy){
	// 				min_energy = temp;
	// 				best_row = 1;
	// 			}

	// 			//branch 2
	// 			temp = get_PLmloop(i,j,k,l) + bp_penalty;
	// 			if(temp < min_energy){
	// 				min_energy = temp;
	// 				best_row=2;
	// 			}

	// 			//branch 3
	// 			if (j>=(i+TURN+1)){
	// 				temp = PfromL.get(i+1,j-1,k,l) + gamma2(j,i);
	// 				if(temp < min_energy){
	// 					min_energy = temp;
	// 					best_row= 3;
	// 				}
	// 			}
	// 		}

	// 		switch (best_row){
	// 			case 1:
	// 				insert_node(i,l,j,k,P_PLiloop);
	// 				break;
	// 			case 2:
	// 				insert_node(i,l,j,k,P_PLmloop);
	// 				break;
	// 			case 3:
	// 				insert_node(i+1,l,j-1,k,P_PfromL);
	// 				// Hosna, Feb 18, 2014
	// 				// filling the structure
	// 				f[i].pair = j;
	// 				f[j].pair = i;
	// 				f[i].type = P_PL;
	// 				f[j].type = P_PL;
	// 				break;
	// 		}
	// 	}
	// 	break;

	// 	case P_PR:
	// 	{
	// 		cand_pos_t i = cur_interval->i;
	// 		cand_pos_t l = cur_interval->j;
	// 		cand_pos_t j = cur_interval->k;
	// 		cand_pos_t k = cur_interval->l;
	// 		if(debug) printf("At P_PR with i=%d,l=%d,j=%d,k=%d\n",i,l,j,k);

	// 		if (!(i <= j && j < k-1 && k <= l)){
	// 			std::cerr << "boder cases: This should not have happened!, P_PR" <<std::endl;
	// 			exit(EXIT_FAILURE);
	// 		}

	// 		if (i<0 ||j<0 || k<0 || l<0 || i>=n || j>=n || k>= n || l>= n){
	// 			std::cerr << "impossible cases: This should not have happened!, P_PR" << std::endl;
	// 			exit(EXIT_FAILURE);
	// 		}

	// 		energy_t min_energy = INF,temp=INF;
	// 		cand_pos_t best_row = -1;
	// 		if (pair[S_[k]][S_[l]]>0){
	// 			//branch 1
	// 			temp = get_PRiloop(i,j,k,l);
	// 			if(temp < min_energy){
	// 				min_energy = temp;
	// 				best_row = 1;
	// 			}

	// 			//branch 2
	// 			temp = get_PRmloop(i,j,k,l)+ bp_penalty;
	// 			if(temp < min_energy){
	// 				min_energy = temp;
	// 				best_row=2;
	// 			}

	// 			//branch 3
	// 			if (l>=(k+TURN+1)){
	// 				temp = PfromR.get(i,j,k+1,l-1) + gamma2(l,k);
	// 				if(temp < min_energy){
	// 					min_energy = temp;
	// 					best_row= 3;
	// 				}
	// 			}
	// 		}

	// 		switch (best_row){
	// 			case 1:
	// 				insert_node(i,l,j,k,P_PRiloop);
	// 				break;
	// 			case 2:
	// 				insert_node(i,l,j,k,P_PRmloop);
	// 				break;
	// 			case 3:
	// 				insert_node(i,l-1,j,k+1,P_PfromR);
	// 				f[k].pair = l;
	// 				f[l].pair = k;
	// 				f[k].type = P_PR;
	// 				f[l].type = P_PR;

	// 				break;
	// 		}

	// 	}
	// 	break;

	// 	case P_PM:
	// 	{
	// 		cand_pos_t i = cur_interval->i;
	// 		cand_pos_t l = cur_interval->j;
	// 		cand_pos_t j = cur_interval->k;
	// 		cand_pos_t k = cur_interval->l;
	// 		if(debug) printf("At P_PM with i=%d,l=%d,j=%d,k=%d\n",i,l,j,k);

	// 		if (!(i <= j && j < k-1 && k <= l)){
	// 			std::cerr << "border cases: This should not have happened!, P_PM" << std::endl;
	// 			exit(EXIT_FAILURE);
	// 		}

	// 		if (i<=0 ||j<=0 || k<=0 || l<=0 || i>n || j>n || k> n || l> n){
	// 			std::cerr << "impossible cases: This should not have happened!, P_PM" << std::endl;
	// 			exit(EXIT_FAILURE);
	// 		}

	// 		if (i==j && k ==l){
	// 			f[j].pair = k;
	// 			f[k].pair = j;
	// 			f[j].type = P_PM;
	// 			f[k].type = P_PM;
	// 			return;
	// 		}

	// 		energy_t min_energy = INF,temp=INF;
	// 		cand_pos_t best_row = -1;

	// 		if (pair[S_[j]][S_[k]]>0){
	// 			//branch 1
	// 			temp = get_PMiloop(i,j,k,l);
	// 			if(temp < min_energy){
	// 				min_energy = temp;
	// 				best_row = 1;
	// 			}
	// 			//branch 2
	// 			temp = get_PMmloop(i,j,k,l) + bp_penalty;
	// 			if(temp < min_energy){
	// 				min_energy = temp;
	// 				best_row=2;
	// 			}
	// 			//branch 3
	// 			if (k>=(j+TURN-1)){
	// 				temp = PfromM.get(i,j-1,k+1,l) + gamma2(j,k);
	// 				if(temp < min_energy){
	// 					min_energy = temp;
	// 					best_row= 3;
	// 				}
	// 			}
	// 		}

	// 		switch (best_row){
	// 			case 1:
	// 				insert_node(i,l,j,k,P_PMiloop);
	// 				break;
	// 			case 2:
	// 				insert_node(i,l,j,k,P_PMmloop);
	// 				break;
	// 			case 3:
	// 				insert_node(i,l,j-1,k+1,P_PfromM);
	// 				f[j].pair = k;
	// 				f[k].pair = j;
	// 				f[j].type = P_PM;
	// 				f[k].type = P_PM;
	// 				break;
	// 		}
	// 	}
	// 	break;

	// 	case P_PO:
	// 	{
	// 		cand_pos_t i = cur_interval->i;
	// 		cand_pos_t l = cur_interval->j;
	// 		cand_pos_t j = cur_interval->k;
	// 		cand_pos_t k = cur_interval->l;
	// 		if(debug) printf("At P_PO with i=%d,l=%d,j=%d,k=%d\n",i,l,j,k);

	// 		if (!(i <= j && j < k-1 && k <= l)){
	// 			std::cerr << "border cases: This should not have happened!, P_PO" << std::endl;
	// 			exit(EXIT_FAILURE);
	// 		}

	// 		if (i<=0 ||j<=0 || k<=0 || l<=0 || i>n || j>n || k> n || l> n){
	// 			std::cerr << "impossible cases: This should not have happened!, P_PO" << std::endl;
	// 			exit(EXIT_FAILURE);
	// 		}

	// 		energy_t min_energy = INF,temp=INF;
	// 		cand_pos_t best_row = -1;
	// 		if (pair[S_[i]][S_[l]] > 0){
	// 			//branch 1
	// 			temp = get_POiloop(i,j,k,l);
	// 			if(temp < min_energy){
	// 				min_energy = temp;
	// 				best_row = 1;
	// 			}
	// 			//branch 2
	// 			temp = get_POmloop(i,j,k,l)+bp_penalty;
	// 			if(temp < min_energy){
	// 				min_energy = temp;
	// 				best_row=2;
	// 			}
	// 			//branch 3
	// 			if (l>=(i+TURN+1)){
	// 				temp = PfromO.get(i+1,j,k,l-1) + gamma2(l,i);
	// 				if(temp < min_energy){
	// 					min_energy = temp;
	// 					best_row= 3;
	// 				}
	// 			}
	// 		}
	// 		switch (best_row){
	// 			case 1:
	// 				insert_node(i,l,j,k,P_POiloop);
	// 				break;
	// 			case 2:
	// 				insert_node(i,l,j,k,P_POmloop);
	// 				break;
	// 			case 3:
	// 				insert_node(i+1,l-1,j,k,P_PfromO);
	// 				f[i].pair = l;
	// 				f[l].pair = i;
	// 				f[i].type = P_PO;
	// 				f[l].type = P_PO;
	// 				break;
	// 		}

	// 	}
	// 	break;

	// 	case P_PfromL:
	// 	{
	// 		cand_pos_t i = cur_interval->i;
	// 		cand_pos_t l = cur_interval->j;
	// 		cand_pos_t j = cur_interval->k;
	// 		cand_pos_t k = cur_interval->l;
	// 		if(debug) printf("At P_PfromL with i=%d,l=%d,j=%d,k=%d\n",i,l,j,k);

	// 		if (!(i <= j && j < k-1 && k <= l)){
	// 			std::cerr << "This should not have happened!, P_PfromL" << std::endl;
	// 			exit(EXIT_FAILURE);
	// 		}

	// 		if (i<=0 ||j<=0 || k<=0 || l<=0 || i>n || j>n || k> n || l> n){
	// 			std::cerr << "This should not have happened!, P_PfromL" << std::endl;
	// 			exit(EXIT_FAILURE);
	// 		}

	// 		if (i==j && k==l){
	// 			return;
	// 		}

	// 		energy_t min_energy = INF,temp=INF;
	// 		cand_pos_t best_row = -1,best_d=-1;

	// 		for(cand_pos_t d=i+1; d< j; d++){
	// 			//branch 1
	// 			temp=PfromL.get(d,j,k,l)+get_WP(i,d-1);
	// 			if(temp < min_energy){
	// 				min_energy=temp;
	// 				best_row=1;
	// 				best_d = d;

	// 			}
	// 			//branch 2
	// 			temp=PfromL.get(i,d,k,l)+get_WP(d+1,j);
	// 			if(temp < min_energy){
	// 				min_energy=temp;
	// 				best_row=2;
	// 				best_d = d;
	// 			}
	// 		}
	// 		// branch 3
	// 		temp = PR.get(i,j,k,l) + gamma2(l,k)+PB_penalty;
	// 		if(temp < min_energy){
	// 			min_energy = temp;
	// 			best_row=3;
	// 			best_d = -1;
	// 		}

	// 		//branch 4
	// 		temp = PM.get(i,j,k,l) + gamma2(j,k)+PB_penalty;
	// 		if(temp < min_energy){
	// 			min_energy = temp;
	// 			best_row=4;
	// 			best_d = -1;
	// 		}

	// 		// branch 5
	// 		temp = PO.get(i,j,k,l) + gamma2(l,i)+PB_penalty;
	// 		if(temp < min_energy){
	// 			min_energy = temp;
	// 			best_row=5;
	// 			best_d = -1;
	// 		}
	// 		switch (best_row){
	// 			case 1:
	// 				if (best_d > -1){
	// 					insert_node(best_d,l,j,k,P_PfromL);
	// 					insert_node(i,best_d-1,P_WP);
	// 				}
	// 				break;
	// 			case 2:
	// 				if (best_d > -1){
	// 					insert_node(i,l,best_d,k,P_PfromL);
	// 					insert_node(best_d+1,j,P_WP);
	// 				}
	// 				break;
	// 			case 3:
	// 				insert_node(i,l,j,k,P_PR);
	// 				break;
	// 			case 4:
	// 				insert_node(i,l,j,k,P_PM);
	// 				break;
	// 			case 5:
	// 				insert_node(i,l,j,k,P_PO);
	// 				break;
	// 		}


	// 	}
	// 		break;

	// 	case P_PfromR:
	// 	{
	// 		cand_pos_t i = cur_interval->i;
	// 		cand_pos_t l = cur_interval->j;
	// 		cand_pos_t j = cur_interval->k;
	// 		cand_pos_t k = cur_interval->l;
	// 		if(debug) printf("At P_PfromR with i=%d,l=%d,j=%d,k=%d\n",i,l,j,k);

	// 		if (!(i <= j && j < k-1 && k <= l)){
	// 			std::cerr << "This should not have happened!, P_PfromR" << std::endl;
	// 			exit(EXIT_FAILURE);
	// 		}

	// 		if (i<=0 ||j<=0 || k<=0 || l<=0 || i>n || j>n || k> n || l> n){
	// 			std::cerr << "impossible case: This should not have happened!, P_PfromR" << std::endl;
	// 			exit(EXIT_FAILURE);
	// 		}

	// 		if (i==j && k==l){
	// 			return;
	// 		}

	// 		energy_t min_energy = INF,temp=INF;
	// 		cand_pos_t best_row = -1,best_d=-1;

	// 		for(cand_pos_t d=k+1; d< l; d++){
	// 			//branch 1
	// 			temp=PfromR.get(i,j,d,l)+get_WP(k,d-1);
	// 			if(temp < min_energy){
	// 				min_energy = temp;
	// 				best_row=1;
	// 				best_d = d;

	// 			}
	// 			//branch 2
	// 			temp=PfromR.get(i,j,k,d)+get_WP(d+1,l);
	// 			if(temp < min_energy){
	// 				min_energy=temp;
	// 				best_row=2;
	// 				best_d = d;
	// 			}
	// 		}

	// 		//branch 3
	// 		temp = PM.get(i,j,k,l) + gamma2(j,k)+PB_penalty;
	// 		if(temp < min_energy){
	// 			min_energy = temp;
	// 			best_row=3;
	// 			best_d = -1;
	// 		}

	// 		// branch 4
	// 		temp = PO.get(i,j,k,l) + gamma2(l,i)+PB_penalty;
	// 		if(temp < min_energy){
	// 			min_energy = temp;
	// 			best_row=4;
	// 			best_d = -1;
	// 		}
	// 		switch (best_row){
	// 			case 1:
	// 				if (best_d > -1){
	// 					insert_node(i,l,j,best_d,P_PfromR);
	// 					insert_node(k,best_d-1,P_WP);
	// 				}
	// 				break;
	// 			case 2:
	// 				if (best_d > -1){
	// 					insert_node(i,best_d,j,k,P_PfromR);
	// 					insert_node(best_d+1,l,P_WP);
	// 				}
	// 				break;
	// 			case 3:
	// 				insert_node(i,l,j,k,P_PM);
	// 				break;
	// 			case 4:
	// 				insert_node(i,l,j,k,P_PO);
	// 				break;
	// 		}


	// 	}
	// 	break;

	// 	case P_PfromM:
	// 	{
	// 		cand_pos_t i = cur_interval->i;
	// 		cand_pos_t l = cur_interval->j;
	// 		cand_pos_t j = cur_interval->k;
	// 		cand_pos_t k = cur_interval->l;
	// 		if(debug) printf("At P_PfromM with i=%d,l=%d,j=%d,k=%d\n",i,l,j,k);

	// 		if (!(i <= j && j < k-1 && k <= l)){
	// 			std::cerr << "This should not have happened!, P_PfromM" << std::endl;
	// 			exit(EXIT_FAILURE);
	// 		}

	// 		if (i<=0 ||j<=0 || k<=0 || l<=0 || i>n || j>n || k> n || l> n){
	// 			std::cerr << "This should not have happened!, P_PfromM" << std::endl;
	// 			exit(EXIT_FAILURE);
	// 		}

	// 		if (i==j && k==l){
	// 			return;
	// 		}

	// 		energy_t min_energy = INF,temp=INF;
	// 		cand_pos_t best_d=-1;

	// 		for(cand_pos_t d=i+1; d< j; d++){
	// 			//branch 1
	// 			temp=PfromMprime.get(i,d,k,l)+get_WP(d+1,j);
	// 			if(temp < min_energy){
	// 				min_energy=temp;
	// 				best_d = d;
	// 			}
	// 		}

	// 		if (best_d > -1){
	// 			insert_node(i,l,best_d,k,P_PfromMprime);
	// 			insert_node(best_d+1,j,P_WP);
	// 		}


	// 	}
	// 		break;

	// 		case P_PfromMprime:
	// 	{
	// 		cand_pos_t i = cur_interval->i;
	// 		cand_pos_t l = cur_interval->j;
	// 		cand_pos_t j = cur_interval->k;
	// 		cand_pos_t k = cur_interval->l;
	// 		if(debug) printf("At P_PfromM with i=%d,l=%d,j=%d,k=%d\n",i,l,j,k);

	// 		if (!(i <= j && j < k-1 && k <= l)){
	// 			std::cerr << "This should not have happened!, P_PfromMprime" << std::endl;
	// 			exit(EXIT_FAILURE);
	// 		}

	// 		if (i<=0 ||j<=0 || k<=0 || l<=0 || i>n || j>n || k> n || l> n){
	// 			std::cerr << "This should not have happened!, P_PfromMprime" << std::endl;
	// 			exit(EXIT_FAILURE);
	// 		}

	// 		if (i==j && k==l){
	// 			return;
	// 		}

	// 		energy_t min_energy = INF,temp=INF;
	// 		cand_pos_t best_d=-1;

	// 		for(cand_pos_t d=k+1; d< l; d++){
	// 			//branch 2
	// 			temp=get_PfromMdoubleprime(i,j,d,l)+get_WP(k,d-1);
	// 			if(temp < min_energy){
	// 				min_energy=temp;
	// 				best_d = d;
	// 			}
	// 		}

	// 		if (best_d > -1){
	// 			insert_node(i,l,j,best_d,P_PfromMdoubleprime);
	// 			insert_node(k,best_d-1,P_WP);
	// 		}

	// 	}
	// 		break;

	// 		case P_PfromMdoubleprime:
	// 	{
	// 		cand_pos_t i = cur_interval->i;
	// 		cand_pos_t l = cur_interval->j;
	// 		cand_pos_t j = cur_interval->k;
	// 		cand_pos_t k = cur_interval->l;
	// 		if(debug) printf("At P_PfromMdoubleprime with i=%d,l=%d,j=%d,k=%d\n",i,l,j,k);

	// 		if (!(i <= j && j < k-1 && k <= l)){
	// 			std::cerr << "This should not have happened!, P_PfromMdoubleprime" << std::endl;
	// 			exit(EXIT_FAILURE);
	// 		}

	// 		if (i<=0 ||j<=0 || k<=0 || l<=0 || i>n || j>n || k> n || l> n){
	// 			std::cerr << "This should not have happened!, P_PfromMdoubleprime" << std::endl;
	// 			exit(EXIT_FAILURE);
	// 		}

	// 		if (i==j && k==l){
	// 			return;
	// 		}

	// 		energy_t min_energy = INF,temp=INF;
	// 		cand_pos_t best_row = -1;

	// 		// branch 1
	// 		temp = PL.get(i,j,k,l) + gamma2(j,i)+PB_penalty;
	// 		if(temp < min_energy){
	// 			min_energy = temp;
	// 			best_row=1;
	// 		}

	// 		//branch 2
	// 		temp = PR.get(i,j,k,l) + gamma2(l,k)+PB_penalty;
	// 		if(temp < min_energy){
	// 			min_energy = temp;
	// 			best_row=2;
	// 		}

	// 		switch (best_row){
	// 			case 1:
	// 				insert_node(i,l,j,k,P_PL);
	// 				break;
	// 			case 2:
	// 				insert_node(i,l,j,k,P_PR);
	// 				break;
	// 		}


	// 	}
	// 		break;

	// 	case P_PfromO:
	// 	{
	// 		cand_pos_t i = cur_interval->i;
	// 		cand_pos_t l = cur_interval->j;
	// 		cand_pos_t j = cur_interval->k;
	// 		cand_pos_t k = cur_interval->l;
	// 		if(debug) printf("At P_PfromO with i=%d,l=%d,j=%d,k=%d\n",i,l,j,k);

	// 		if (!(i <= j && j < k-1 && k <= l)){
	// 			std::cerr << "border cases: This should not have happened!, P_PfromO" << std::endl;
	// 			exit(EXIT_FAILURE);
	// 		}

	// 		if (i<=0 ||j<=0 || k<=0 || l<=0 || i>n || j>n || k> n || l> n){
	// 			std::cerr << "impossible case: This should not have happened!, P_PfromO" << std::endl;
	// 			exit(EXIT_FAILURE);
	// 		}

	// 		if (i==j && k==l){
	// 			return;
	// 		}

	// 		energy_t min_energy = INF,temp=INF;
	// 		cand_pos_t best_row = -1,best_d=-1;

	// 		for(cand_pos_t d=i+1; d< j; d++){
	// 			//branch 1
	// 			temp=PfromO.get(d,j,k,l)+get_WP(i,d-1);
	// 			if(temp < min_energy){
	// 				min_energy=temp;
	// 				best_row=1;
	// 				best_d = d;

	// 			}
	// 		}
	// 		for(cand_pos_t d=k+1; d< l; d++){
	// 			//branch 2
	// 			temp=PfromO.get(i,j,k,d)+get_WP(d+1,l);
	// 			if(temp < min_energy){
	// 				min_energy=temp;
	// 				best_row=2;
	// 				best_d = d;
	// 			}
	// 		}
	// 		// branch 3
	// 		temp = PL.get(i,j,k,l) + gamma2(j,i)+PB_penalty;
	// 		if(temp < min_energy){
	// 			min_energy = temp;
	// 			best_row=3;
	// 			best_d = -1;
	// 		}

	// 		//branch 4
	// 		temp = PR.get(i,j,k,l) + gamma2(l,k)+PB_penalty;
	// 		if(temp < min_energy){
	// 			min_energy = temp;
	// 			best_row=4;
	// 			best_d = -1;
	// 		}

	// 		switch (best_row){
	// 			case 1:
	// 				if (best_d > -1){
	// 					insert_node(best_d,l,j,k,P_PfromO);
	// 					insert_node(i,best_d-1,P_WP);
	// 				}
	// 				break;
	// 			case 2:
	// 				if (best_d > -1){
	// 					insert_node(i,best_d,j,k,P_PfromO);
	// 					insert_node(best_d+1,l,P_WP);
	// 				}
	// 				break;
	// 			case 3:
	// 				insert_node(i,l,j,k,P_PL);
	// 				break;
	// 			case 4:
	// 				insert_node(i,l,j,k,P_PR);
	// 				break;
	// 		}


	// 	}
	// 		break;
	// 	case P_WB:
	// 	{
	// 		cand_pos_t i = cur_interval->i;
	// 		cand_pos_t l = cur_interval->j;
	// 		if(debug) printf("At P_WM with i=%d,l=%d\n",i,l);

	// 		if (i<=0 ||l<=0 || i>n || l >n){
	// 			std::cerr << "impossible cases: This should not have happened!, P_WB" << std::endl;
	// 			exit(EXIT_FAILURE);
	// 		}

	// 		if (i>l){
	// 			return;
	// 		}

	// 		energy_t min_energy = INF,temp=INF;
	// 		cand_pos_t best_row = -1;
	// 		//branch 1
	// 		temp = WBP.get(i,l);
	// 		if (temp < min_energy){
	// 			min_energy = temp;
	// 			best_row=1;
	// 		}

	// 		temp = cp_penalty*(l-i+1);
	// 		if (temp < min_energy){
	// 			min_energy = temp;
	// 			best_row=2;
	// 		}

	// 		switch (best_row){
	// 			case 1:
	// 				insert_node(i,l,P_WBP);
	// 				break;
	// 			case 2:
	// 				// do nothing.
	// 				break;
	// 		}

	// 	}
	// 		break;
	// 	case P_WBP:
	// 	{
	// 		cand_pos_t i = cur_interval->i;
	// 		cand_pos_t l = cur_interval->j;
	// 		if(debug) printf("At P_WBP with i=%d,l=%d\n",i,l);
	// 		if (i>l){
	// 			std::cerr << "border case: This should not have happened!, P_WBP" << std::endl;
	// 			exit(EXIT_FAILURE);
	// 		}

	// 		if (i<=0 ||l<=0 || i>n || l >n){
	// 			std::cerr << "impossible cases: This should not have happened!, P_WBP" << std::endl;
	// 			exit(EXIT_FAILURE);
	// 		}

	// 		energy_t min_energy = INF,temp=INF;
	// 		cand_pos_t best_row = -1,best_d=-1;

	// 		for(cand_pos_t d=i; d< l; d++){
	// 			//branch 1
	// 			temp = get_WB(i,d-1) + V->get_energy(d,l) +beta2P(l,d) + PPS_penalty;
	// 			if (temp < min_energy){
	// 				min_energy = temp;
	// 				best_row = 1;
	// 				best_d = d;
	// 			}

	// 			//branch 2
	// 			temp = get_WB(i,d-1) + P.get(d,l) + PSM_penalty + PPS_penalty;
	// 			if (temp < min_energy){
	// 				min_energy = temp;
	// 				best_row = 2;
	// 				best_d = d;
	// 			}
	// 		}
	// 		temp = WBP.get(i,l-1) + cp_penalty;
	// 		if(temp< min_energy){
	// 			min_energy = temp;
	// 			best_row = 3;
	// 		}

	// 		switch (best_row){
	// 			case 1:
	// 				insert_node(i,best_d-1,P_WB);
	// 				insert_node(best_d,l,LOOP);
	// 				break;
	// 			case 2:
	// 				insert_node(i,best_d-1,P_WB);
	// 				insert_node(best_d,l,P_P);
	// 				break;
	// 			case 3:
	// 				insert_node(i,l-1,P_WBP);
	// 				break;
	// 		}
	// 	}
	// 		break;

	// 	case P_WP:
	// 	{
	// 		cand_pos_t i = cur_interval->i;
	// 		cand_pos_t l = cur_interval->j;
	// 		if(debug) printf("At P_WP with i=%d,l=%d\n",i,l);

	// 		if (i<=0 ||l<=0 || i>n || l >n){
	// 			std::cerr << "impossible cases: This should not have happened!, P_WP" << std::endl;
	// 			exit(EXIT_FAILURE);
	// 		}

	// 		if (i>l){
	// 			return;
	// 		}

	// 		energy_t min_energy = INF,temp=INF;
	// 		cand_pos_t best_row = -1;
	// 		//branch 1
	// 		temp = WPP.get(i,l);
	// 		if (temp < min_energy){
	// 			min_energy = temp;
	// 			best_row=1;
	// 		}

	// 		temp = PUP_penalty*(l-i+1);
	// 		if (temp < min_energy){
	// 			min_energy = temp;
	// 			best_row=2;
	// 		}

	// 		switch (best_row){
	// 			case 1:
	// 				insert_node(i,l,P_WPP);
	// 				break;
	// 			case 2:
	// 				// do nothing.
	// 				break;
	// 		}

	// 	}
	// 		break;
	// 	case P_WPP:
	// 	{
	// 		cand_pos_t i = cur_interval->i;
	// 		cand_pos_t l = cur_interval->j;
	// 		if(debug) printf("At P_WPP with i=%d,l=%d\n",i,l);
	// 		if (i>l){
	// 			std::cerr << "border case: This should not have happened!, P_WPP" << std::endl;
	// 			exit(EXIT_FAILURE);
	// 		}

	// 		if (i<=0 ||l<=0 || i>n || l >n){
	// 			std::cerr << "impossible cases: This should not have happened!, P_WPP" << std::endl;
	// 			exit(EXIT_FAILURE);
	// 		}

	// 		energy_t min_energy = INF,temp=INF;
	// 		cand_pos_t best_row = -1,best_d=-1;

	// 		for(cand_pos_t d=i; d< l; d++){
	// 			//branch 1
	// 			temp = get_WP(i,d-1) + V->get_energy(d,l) + gamma2(l,d) + PPS_penalty;
	// 			if (temp < min_energy){
	// 				min_energy = temp;
	// 				best_row = 1;
	// 				best_d = d;
	// 			}

	// 			//branch 2
	// 			temp = get_WP(i,d-1) + P.get(d,l) + PSP_penalty + PPS_penalty;
	// 			if (temp < min_energy){
	// 				min_energy = temp;
	// 				best_row = 2;
	// 				best_d = d;
	// 			}
	// 		}
	// 		temp = WPP.get(i,l-1) + PUP_penalty;
	// 		if(temp<min_energy){
	// 			min_energy = temp;
	// 			best_row = 3;
	// 		}

	// 		switch (best_row){
	// 			case 1:
	// 				insert_node(i,best_d-1,P_WP);
	// 				insert_node(best_d,l,LOOP);
	// 				break;
	// 			case 2:
	// 				insert_node(i,best_d-1,P_WP);
	// 				insert_node(best_d,l,P_P);
	// 				break;
	// 			case 3:
	// 				insert_node(i,l-1,P_WPP);
	// 				break;
	// 		}
	// 	}
	// 		break;
	// 	case P_PLiloop:
	// 	{
	// 		cand_pos_t i = cur_interval->i;
	// 		cand_pos_t l = cur_interval->j;
	// 		cand_pos_t j = cur_interval->k;
	// 		cand_pos_t k = cur_interval->l;
	// 		if(debug) printf("At P_PLiloop with i=%d,l=%d,j=%d,k=%d\n",i,l,j,k);

	// 		if (!(i < j && j < k-1 && k < l)){
	// 			std::cerr << "border cases: This should not have happened!, P_PLiloop" << std::endl;
	// 			exit(EXIT_FAILURE);
	// 		}

	// 		if (i<=0 || j<=0 || k<=0 || l<=0 || i>n || j>n || k> n || l> n){
	// 			std::cerr << "impossbible cases: This should not have happened!, P_PLiloop" << std::endl;
	// 			exit(EXIT_FAILURE);
	// 		}

	// 		f[i].pair = j;
	// 		f[j].pair = i;
	// 		f[i].type = P_PLiloop;
	// 		f[j].type = P_PLiloop;

	// 		energy_t min_energy = INF,temp=INF;
	// 		cand_pos_t best_row = -1,best_d=-1,best_dp=-1;
	// 		if (pair[S_[i]][S_[j]] > 0){
	// 			//branch 1
	// 			temp = PL.get(i+1,j-1,k,l) + get_e_stP(i,j);
	// 			if (temp < min_energy){
	// 				min_energy = temp;
	// 				best_row=1;
	// 			}
	// 			//branch 2
	// 			cand_pos_t max_d = std::min(j,i+MAXLOOP);
	// 			for(cand_pos_t d= i+1; d<max_d; ++d){
	// 				cand_pos_t min_dp = std::max(d+TURN,j-MAXLOOP);
	// 				for(cand_pos_t dp = j-1; dp > min_dp; --dp){
	// 					temp = get_e_intP(i,d,dp,j) + PL.get(d,dp,k,l);
	// 					if(temp < min_energy){
	// 						min_energy = temp;
	// 						best_d = d;
	// 						best_dp = dp;
	// 						best_row = 2;
	// 					}
	// 				}
    //         }	
				
	// 		}

	// 		switch (best_row){
	// 			case 1:
	// 				insert_node(i+1,l,j-1,k,P_PL);
	// 				break;
	// 			case 2:
	// 				insert_node(best_d,l,best_dp,k,P_PL);
	// 				break;
	// 		}
	// 	}
	// 		break;

	// 	case P_PLmloop:
	// 	{
	// 		cand_pos_t i = cur_interval->i;
	// 		cand_pos_t l = cur_interval->j;
	// 		cand_pos_t j = cur_interval->k;
	// 		cand_pos_t k = cur_interval->l;
	// 		if(debug) printf("At P_PLmloop with i=%d,l=%d,j=%d,k=%d\n",i,l,j,k);

	// 		if (!(i <= j && j < k-1 && k <= l)){
	// 			std::cerr << "border cases: This should not have happened!, P_PLmloop" << std::endl;
	// 			exit(EXIT_FAILURE);
	// 		}

	// 		if (i<=0 || j<=0 || k<=0 || l<=0 || i>n || j>n || k> n || l> n){
	// 			std::cerr << "impossible cases: This should not have happened!, P_PLmloop" << std::endl;
	// 			exit(EXIT_FAILURE);
	// 		}

	// 		f[i].pair = j;
	// 		f[j].pair = i;
	// 		f[i].type = P_PLmloop;
	// 		f[j].type = P_PLmloop;

	// 		energy_t branch1 = PLmloop10.get(i+1,j-1,k,l)+ ap_penalty + beta2P(j,i);
    //         energy_t branch2 = PLmloop01.get(i+1,j-1,k,l)+ ap_penalty + beta2P(j,i);

    //         cand_pos_t best_row = -1;
    //         if (branch1 < branch2) { best_row = 1;} else {best_row = 2;}

	// 		switch (best_row){
	// 			case 1:
	// 				insert_node(i+1,l,j-1,k,P_PLmloop10);
	// 				break;
	// 			case 2:
	// 				insert_node(i+1,l,j-1,k,P_PLmloop01);
	// 				break;
	// 		}

	// 	}
	// 	break;
	// 	case P_PLmloop00:
	// 	{
	// 		cand_pos_t i = cur_interval->i;
	// 		cand_pos_t l = cur_interval->j;
	// 		cand_pos_t j = cur_interval->k;
	// 		cand_pos_t k = cur_interval->l;
	// 		if(debug) printf("At P_PLmloop00 with i=%d,l=%d,j=%d,k=%d\n",i,l,j,k);

	// 		if (!(i <= j && j < k-1 && k <= l)){
	// 			std::cerr << "border cases: This should not have happened!, P_PLmloop00" << std::endl;
	// 			exit(EXIT_FAILURE);
	// 		}

	// 		if (i<=0 || j<=0 || k<=0 || l<=0 || i>n || j>n || k> n || l> n){
	// 			std::cerr << "impossible cases: This should not have happened!, P_PLmloop00" << std::endl;
	// 			exit(EXIT_FAILURE);
	// 		}


	// 		energy_t min_energy = PL.get(i,j,k,l)+beta2P(j,i),temp=INF;
	// 		cand_pos_t best_row = 1, best_d=-1;
    //         for(cand_pos_t d = i; d<=j; ++d){
    //             if (d>i){
    //                 temp=get_WB(i,d-1)+PLmloop00.get(d,j,k,l);
    //                 if (temp < min_energy){
    //                     min_energy = temp;
    //                     best_row = 2;
    //                     best_d = d;
    //                 }

    //             }
    //             if(d<j){
    //                 temp=PLmloop00.get(i,d,k,l)+get_WB(d+1,j);
    //                 if (temp < min_energy){
    //                     min_energy = temp;
    //                     best_row = 3;
    //                     best_d = d;
    //                 }
    //             }
    //         }
	// 		switch (best_row){
    //             case 1:
    //                 insert_node(i,l,j,k,P_PL);
    //                 break;
	// 			case 2:
	// 				insert_node(best_d,l,j,k,P_PLmloop00);
	// 				insert_node(i,best_d-1,P_WB);
	// 				break;
	// 			case 3:
	// 				insert_node(i,l,best_d,k,P_PLmloop00);
	// 				insert_node(best_d+1,j,P_WB);
	// 				break;
	// 		}

	// 	}
	// 		break;
	// 	case P_PLmloop01:
	// 	{
	// 		cand_pos_t i = cur_interval->i;
	// 		cand_pos_t l = cur_interval->j;
	// 		cand_pos_t j = cur_interval->k;
	// 		cand_pos_t k = cur_interval->l;
	// 		if(debug) printf("At P_PLmloop01 with i=%d,l=%d,j=%d,k=%d\n",i,l,j,k);

	// 		if (!(i <= j && j < k-1 && k <= l)){
	// 			std::cerr << "border cases: This should not have happened!, P_PLmloop01" << std::endl;
	// 			exit(EXIT_FAILURE);
	// 		}

	// 		if (i<=0 || j<=0 || k<=0 || l<=0 || i>n || j>n || k> n || l> n){
	// 			std::cerr << "impossible cases: This should not have happened!, P_PLmloop01" << std::endl;
	// 			exit(EXIT_FAILURE);
	// 		}

	// 		energy_t min_energy = INF,temp=INF;
	// 		cand_pos_t best_d=-1;
	// 		for(cand_pos_t d = i; d < j; ++d){
    //             temp = PLmloop00.get(i,d,k,l) + WBP.get(d+1,j);
    //             if (temp < min_energy){
    //                 min_energy = temp;
    //                 best_d = d;
    //             }
    //         }
	// 		insert_node(i,l,best_d,k,P_PLmloop00);
	// 		insert_node(best_d+1,j,P_WBP);

	// 	}
	// 		break;
	// 	case P_PLmloop10:
	// 	{
	// 		cand_pos_t i = cur_interval->i;
	// 		cand_pos_t l = cur_interval->j;
	// 		cand_pos_t j = cur_interval->k;
	// 		cand_pos_t k = cur_interval->l;
	// 		if(debug) printf("At P_PLmloop10 with i=%d,l=%d,j=%d,k=%d\n",i,l,j,k);

	// 		if (!(i <= j && j < k-1 && k <= l)){
	// 			std::cerr<< "border cases: This should not have happened!, P_PLmloop10" << std::endl;
	// 			exit(EXIT_FAILURE);
	// 		}

	// 		if (i<=0 || j<=0 || k<=0 || l<=0 || i>n || j>n || k> n || l> n){
	// 			std::cerr<< "impossible cases: This should not have happened!, P_PLmloop10" << std::endl;
	// 			exit(EXIT_FAILURE);
	// 		}

	// 		energy_t min_energy = INF,temp=INF;
	// 		cand_pos_t best_d=-1,best_row=-1;
	// 		for(cand_pos_t d = i+1; d <= j; ++d){
    //             temp = WBP.get(i,d-1) + PLmloop00.get(d,j,k,l);
    //             if (temp < min_energy){
    //                 min_energy = temp;
    //                 best_row = 1;
    //                 best_d = d;
    //             }
    //             if(d<j){
    //                 temp = PLmloop10.get(i,d,k,l) + get_WB(d+1,j);
    //                 if (temp < min_energy){
    //                     min_energy = temp;
    //                     best_row = 2;
    //                     best_d = d;
    //                 }
    //             }
    //         }
	// 		switch (best_row){
	// 			case 1:
	// 				insert_node(i, best_d-1, P_WBP);
    //                 insert_node(best_d,l,j,k, P_PLmloop00);
	// 				break;
	// 			case 2:
	// 				insert_node(i,l,best_d,k, P_PLmloop10);
	// 				insert_node(best_d+1,j,P_WB);
	// 				break;
	// 		}
	// 	}
	// 		break;

	// 	case P_PRiloop:
	// 	{
	// 		cand_pos_t i = cur_interval->i;
	// 		cand_pos_t l = cur_interval->j;
	// 		cand_pos_t j = cur_interval->k;
	// 		cand_pos_t k = cur_interval->l;
	// 		if(debug) printf("At P_PRiloop with i=%d,l=%d,j=%d,k=%d\n",i,l,j,k);

	// 		if (!(i <= j && j < k-1 && k <= l)){
	// 			std::cerr << "border cases: This should not have happened!, P_PRiloop" << std::endl;
	// 			exit(EXIT_FAILURE);
	// 		}

	// 		if (i<=0 || j<=0 || k<=0 || l<=0 || i>n || j>n || k> n || l> n){
	// 			std::cerr << "impossible cases: This should not have happened!, P_PRiloop" << std::endl;
	// 			exit(EXIT_FAILURE);
	// 		}

	// 		f[k].pair = l;
	// 		f[l].pair = k;
	// 		f[k].type = P_PRiloop;
	// 		f[l].type = P_PRiloop;

	// 		energy_t min_energy = INF,temp=INF;
	// 		cand_pos_t best_row=-1,best_d=-1,best_dp=-1;
	// 		if (pair[S_[k]][S_[l]]>0){
	// 			//branch 1
	// 			temp = PR.get(i,j,k+1,l-1) + get_e_stP(k,l);
	// 			if(temp < min_energy){
	// 				min_energy = temp;
	// 				best_row = 1;
	// 			}
	// 			//branch 2
	// 			cand_pos_t max_d = std::min(l,k+MAXLOOP);
	// 			for(cand_pos_t d= k+1; d<max_d; ++d){
	// 				cand_pos_t min_dp = std::max(d+TURN,l-MAXLOOP);
	// 				for(cand_pos_t dp=l-1; dp > min_dp; --dp){
	// 					temp = get_e_intP(k,d,dp,l) + PR.get(i,j,d,dp);
	// 					if(temp < min_energy){
	// 						min_energy = temp;
	// 						best_d = d;
	// 						best_dp = dp;
	// 						best_row = 2;
	// 					}
	// 				}
	// 			}
	// 		}

	// 		switch (best_row)
	// 		{
	// 			case 1:
	// 				insert_node(i,l-1,j,k+1,P_PR);
	// 				break;
	// 			case 2:
	// 				insert_node(i,best_dp,j,best_d,P_PR);
	// 				break;
	// 		}

	// 	}
	// 	break;

	// 	case P_PRmloop:
	// 	{
	// 		cand_pos_t i = cur_interval->i;
	// 		cand_pos_t l = cur_interval->j;
	// 		cand_pos_t j = cur_interval->k;
	// 		cand_pos_t k = cur_interval->l;
	// 		if(debug) printf("At P_PRmloop with i=%d,l=%d,j=%d,k=%d\n",i,l,j,k);

	// 		if (!(i <= j && j < k-1 && k <= l)){
	// 			std::cerr << "border cases: This should not have happened!, P_PRmloop" << std::endl;
	// 			exit(EXIT_FAILURE);
	// 		}
	// 		if (i<=0 || j<=0 || k<=0 || l<=0 || i>n || j>n || k> n || l> n){
	// 			std::cerr << "impossible cases: This should not have happened!, P_PRmloop" << std::endl;
	// 			exit(EXIT_FAILURE);
	// 		}

	// 		f[k].pair = l;
	// 		f[l].pair = k;
	// 		f[k].type = P_PRmloop;
	// 		f[l].type = P_PRmloop;

	// 		cand_pos_t best_row=-1;

	// 		energy_t branch1 = PRmloop10.get(i,j,k+1,l-1) + ap_penalty + beta2P(l,k);
    //         energy_t branch2 = PRmloop01.get(i,j,k+1,l-1) + ap_penalty + beta2P(l,k);

    //         if (branch1 < branch2) {best_row = 1; } else {best_row = 2;}

	// 		switch (best_row){
	// 			case 1:
	// 				insert_node(i,l-1,j,k+1,P_PRmloop10);
	// 				break;
	// 			case 2:
	// 				insert_node(i,l-1,j,k+1,P_PRmloop01);
	// 				break;
	// 		}
	// 	}
	// 	break;

	// 	case P_PRmloop00: // Mateo - This was omitted in the sparse code; not sure why.
	// 	{
	// 		cand_pos_t i = cur_interval->i;
	// 		cand_pos_t l = cur_interval->j;
	// 		cand_pos_t j = cur_interval->k;
	// 		cand_pos_t k = cur_interval->l;
	// 		if(debug) printf("At P_PRmloop00 with i=%d,l=%d,j=%d,k=%d\n",i,l,j,k);

	// 		if (!(i <= j && j < k-1 && k <= l)){
	// 			std::cerr << "border cases: This should not have happened!, P_PRmloop00" << std::endl;
	// 			exit(EXIT_FAILURE);
	// 		}

	// 		if (i<=0 || j<=0 || k<=0 || l<=0 || i>n || j>n || k> n || l> n){
	// 			std::cerr << "impossible cases: This should not have happened!, P_PRmloop00" << std::endl;
	// 			exit(EXIT_FAILURE);
	// 		}

	// 		energy_t min_energy = INF,temp=INF;
	// 		cand_pos_t best_d=-1;
    //         min_energy = PR.get(i,j,k,l)+beta2P(l,k);
	// 		cand_pos_t best_row=1;
	// 		for(cand_pos_t d=k; d<=l; ++d){
	// 			if(d>k){
	// 				temp=get_WB(k,d-1)+PRmloop00.get(i,j,d,l);
	// 				if (temp < min_energy){
	// 					min_energy = temp;
	// 					best_row = 2;
	// 					best_d = d;
    //             	}
	// 			}
	// 			if (d<l){
	// 				temp = PRmloop00.get(i,j,k,d)+get_WB(d+1,l);
	// 				if (temp < min_energy){
	// 					min_energy = temp;
	// 					best_row = 3;
	// 					best_d = d;
    //             	}
	// 			}
	// 		}
	// 		switch (best_row){
    //             case 1:
    //                 insert_node(i,j,k,l,P_PR);
    //                 break;
	// 			case 2:
	// 				insert_node(i,j,best_d,l,P_PRmloop00);
	// 				insert_node(k,best_d-1,P_WB);
	// 				break;
	// 			case 3:
	// 				insert_node(i,j,k,best_d,P_PRmloop00);
	// 				insert_node(best_d+1,l,P_WB);
	// 				break;
	// 		}

	// 	}
	// 	break;

	// 	case P_PRmloop01:
	// 	{
	// 		cand_pos_t i = cur_interval->i;
	// 		cand_pos_t l = cur_interval->j;
	// 		cand_pos_t j = cur_interval->k;
	// 		cand_pos_t k = cur_interval->l;
	// 		if(debug) printf("At P_PRmloop01 with i=%d,l=%d,j=%d,k=%d\n",i,l,j,k);

	// 		if (!(i <= j && j < k-1 && k <= l)){
	// 			std::cerr << "border cases: This should not have happened!, P_PRmloop01" << std::endl;
	// 			exit(EXIT_FAILURE);
	// 		}

	// 		if (i<=0 || j<=0 || k<=0 || l<=0 || i>n || j>n || k> n || l> n){
	// 			std::cerr << "impossible cases: This should not have happened!, P_PRmloop01" << std::endl;
	// 			exit(EXIT_FAILURE);
	// 		}

	// 		energy_t min_energy = INF,temp=INF;
	// 		cand_pos_t best_d=-1;
    //         min_energy = PRmloop01.get(i,j,k,l-1) + cp_penalty;
    //         cand_pos_t best_row = 1;

    //         for(cand_pos_t d = k; d < l; ++d){
    //             temp = PRmloop00.get(i,j,k,d) + WBP.get(d+1,l);
    //             if (temp < min_energy){
    //                 min_energy = temp;
    //                 best_row = 2;
    //                 best_d = d;
    //             }
    //         }
	// 		switch (best_row) {
    //             case 1:
    //                 insert_node(i,l-1,j,k,P_PRmloop01);
    //                 break;
    //             case 2:
    //                 insert_node(best_d+1,l,P_WBP);
    //                 insert_node(i,best_d,j,k,P_PRmloop00);
    //                 break;
    //         }

	// 	}
	// 	break;

	// 	case P_PRmloop10:
	// 	{
	// 		cand_pos_t i = cur_interval->i;
	// 		cand_pos_t l = cur_interval->j;
	// 		cand_pos_t j = cur_interval->k;
	// 		cand_pos_t k = cur_interval->l;
	// 		if(debug) printf("At P_PRmloop10 with i=%d,l=%d,j=%d,k=%d\n",i,l,j,k);

	// 		if (!(i <= j && j < k-1 && k <= l)){
	// 			std::cerr << "border cases: This should not have happened!, P_PRmloop10" << std::endl;
	// 			exit(EXIT_FAILURE);
	// 		}

	// 		if (i<=0 || j<=0 || k<=0 || l<=0 || i>n || j>n || k> n || l> n){
	// 			std::cerr << "impossible cases: This should not have happened!, P_PRmloop10" << std::endl;
	// 			exit(EXIT_FAILURE);
	// 		}

	// 		energy_t min_energy = INF,temp=INF; 
	// 		cand_pos_t best_d=-1;
	// 		min_energy = PRmloop10.get(i,j,k+1,l)+ cp_penalty;
    //         cand_pos_t best_row = 1;

    //         for(cand_pos_t d = k+1; d <= l; ++d){
    //             temp = WBP.get(k,d-1) + PRmloop00.get(i,j,d,l);
    //             if (temp < min_energy){
    //                 min_energy = temp;
    //                 best_row = 2;
    //                 best_d = d;
    //             }
    //         }
	// 		switch (best_row) {
    //             case 1:
    //                 insert_node(i,l,j,k+1,P_PRmloop10);
    //                 break;
    //             case 2:
    //                 insert_node(k,best_d-1,P_WBP);
    //                 insert_node(i,l,j,best_d,P_PRmloop00);
    //                 break;
    //         }

	// 	}
	// 	break;

	// 	case P_PMiloop:
	// 	{
	// 		cand_pos_t i = cur_interval->i;
	// 		cand_pos_t l = cur_interval->j;
	// 		cand_pos_t j = cur_interval->k;
	// 		cand_pos_t k = cur_interval->l;
	// 		if(debug) printf("At P_PMiloop with i=%d,l=%d,j=%d,k=%d\n",i,l,j,k);

	// 		if (!(i <= j && j < k-1 && k <= l)){
	// 			std::cerr << "border cases: This should not have happened!, P_PMiloop" << std::endl;
	// 			exit(EXIT_FAILURE);
	// 		}

	// 		if (i<=0 || j<=0 || k<=0 || l<=0 || i>n || j>n || k> n || l> n){
	// 			std::cerr << "impossible cases: This should not have happened!, P_PMiloop" << std::endl;
	// 			exit(EXIT_FAILURE);
	// 		}

	// 		f[j].pair = k;
	// 		f[k].pair = j;
	// 		f[j].type = P_PMiloop;
	// 		f[k].type = P_PMiloop;
			
	// 		energy_t min_energy = INF,temp=INF;
	// 		cand_pos_t best_d=-1,best_dp=-1,best_row=-1;

	// 		if (pair[S_[j]][S_[k]]>0){
	// 				// branch 1
	// 			temp = PM.get(i,j-1,k+1,l) + get_e_stP(j-1,k+1);
	// 			if (temp < min_energy){
	// 				min_energy = temp;
	// 				best_row=1;
	// 			}

	// 			cand_pos_t max_d = std::max(i,j-MAXLOOP);
	// 			for(cand_pos_t d= j-1; d>max_d; --d){
	// 				cand_pos_t min_dp = std::min(l,k+MAXLOOP);
	// 				for (cand_pos_t dp=k+1; dp <min_dp; ++dp) {
	// 					temp = get_e_intP(d,j,k,dp) + PM.get(i,d,dp,l);
	// 					if(temp < min_energy){
	// 						min_energy = temp;
	// 						best_d = d;
	// 						best_dp = dp;
	// 						best_row = 2;
	// 					}
	// 				}
	// 			}
	// 		}

	// 		switch (best_row)
	// 		{
	// 			case 1:
	// 				insert_node(i,l,j-1,k+1,P_PM);
	// 				break;
	// 			case 2:
	// 				insert_node(i,l,best_d,best_dp,P_PM);
	// 				break;
	// 		}
	// 	}
	// 		break;
	// 	case P_PMmloop:
	// 	{
	// 		cand_pos_t i = cur_interval->i;
	// 		cand_pos_t l = cur_interval->j;
	// 		cand_pos_t j = cur_interval->k;
	// 		cand_pos_t k = cur_interval->l;
	// 		if(debug) printf("At P_PMmloop with i=%d,l=%d,j=%d,k=%d\n",i,l,j,k);

	// 		if (!(i <= j && j < k-1 && k <= l)){
	// 			std::cerr << "border cases: This should not have happened!, P_PMmloop" << std::endl;
	// 			exit(EXIT_FAILURE);
	// 		}

	// 		if (i<=0 || j<=0 || k<=0 || l<=0 || i>n || j>n || k> n || l> n){
	// 			std::cerr << "impossible cases: This should not have happened!, P_PMmloop" << std::endl;
	// 			exit(EXIT_FAILURE);
	// 		}

	// 		f[j].pair = k;
	// 		f[k].pair = j;
	// 		f[j].type = P_PMmloop;
	// 		f[k].type = P_PMmloop;

	// 		energy_t branch1 = PMmloop10.get(i,j-1,k+1,l)+ap_penalty+beta2P(j,k);
    //         energy_t branch2 = PMmloop01.get(i,j-1,k+1,l)+ap_penalty+beta2P(j,k);

    //         cand_pos_t best_row = -1;
    //         if (branch1 < branch2) { best_row = 1;}  else { best_row = 2;}

	// 		switch (best_row){
	// 			case 1:
	// 				insert_node(i,l,j-1,k+1,P_PMmloop10);
	// 				break;
	// 			case 2:
	// 				insert_node(i,l,j-1,k+1,P_PMmloop01);
	// 				break;
	// 		}
	// 	}
	// 	break;
	// 	case P_PMmloop00:
	// 	{
	// 		cand_pos_t i = cur_interval->i;
	// 		cand_pos_t l = cur_interval->j;
	// 		cand_pos_t j = cur_interval->k;
	// 		cand_pos_t k = cur_interval->l;
	// 		if(debug) printf("At P_PMmloop00 with i=%d,l=%d,j=%d,k=%d\n",i,l,j,k);

	// 		if (!(i <= j && j < k-1 && k <= l)){
	// 			std::cerr << "border cases: This should not have happened!, P_PMmloop00" << std::endl;
	// 			exit(EXIT_FAILURE);
	// 		}

	// 		if (i<=0 || j<=0 || k<=0 || l<=0 || i>n || j>n || k> n || l> n){
	// 			std::cerr << "impossible cases: This should not have happened!, P_PMmloop00" << std::endl;
	// 			exit(EXIT_FAILURE);
	// 		}

	// 		f[j].pair = k;
	// 		f[k].pair = j;
	// 		f[j].type = P_PMmloop;
	// 		f[k].type = P_PMmloop;

	// 		energy_t temp=INF;
	// 		energy_t min_energy = PM.get(i,j,k,l)+beta2P(j,k);
	// 		cand_pos_t best_row = 1, best_d = -1;

    //         for(cand_pos_t d=i; d<j; ++d){
    //             temp=get_WB(d+1,j)+PMmloop00.get(i,d,k,l);
    //             if (temp < min_energy){
    //                 min_energy = temp;
    //                 best_row = 2;
    //                 best_d = d;
    //             }
    //         }
    //         for(cand_pos_t d=k+1; d<=l; d++){
    //             temp=PMmloop00.get(i,j,d,l)+get_WB(k,d-1);
    //             if (temp < min_energy){
    //                 min_energy = temp;
    //                 best_row = 3;
    //                 best_d = d;
    //             }
    //         }
			

	// 		switch (best_row){
    //             case 1:
    //             	insert_node(i,l,j,k,P_PM);
    //             	break;
    //             case 2:
    //             	insert_node(i,l,best_d,k,P_PMmloop00);
    //             	insert_node(best_d+1,j,P_WB);
    //             	break;
    //             case 3:
    //             	insert_node(i,l,j,best_d,P_PMmloop00);
    //             	insert_node(k,best_d-1,P_WB);
    //             	break;
    //         }

	// 	}
	// 		break;

	// 	case P_PMmloop01:
	// 	{
	// 		cand_pos_t i = cur_interval->i;
	// 		cand_pos_t l = cur_interval->j;
	// 		cand_pos_t j = cur_interval->k;
	// 		cand_pos_t k = cur_interval->l;
	// 		if(debug) printf("At P_PMmloop01 with i=%d,l=%d,j=%d,k=%d\n",i,l,j,k);

	// 		if (!(i <= j && j < k-1 && k <= l)){
	// 			std::cerr << "border cases: This should not have happened!, P_PMmloop01" << std::endl;
	// 			exit(EXIT_FAILURE);
	// 		}

	// 		if (i<=0 || j<=0 || k<=0 || l<=0 || i>n || j>n || k> n || l> n){
	// 			std::cerr << "impossible cases: This should not have happened!, P_PMmloop01" << std::endl;
	// 			exit(EXIT_FAILURE);
	// 		}

	// 		energy_t temp=INF;
	// 		energy_t min_energy = PMmloop01.get(i,j,k+1,l)+cp_penalty;
	// 		cand_pos_t best_row = 1, best_d = -1;

    //         for(cand_pos_t d = k+1; d <= l; ++d){
    //             temp = PMmloop00.get(i,j,d,l) + WBP.get(k,d-1);
    //             if (temp < min_energy){
    //                 min_energy = temp;
    //                 best_row = 2;
    //                 best_d = d;
    //             }
    //         }
	// 		switch (best_row){
    //             case 1:
    //             	insert_node(i,l,j,k+1,P_PMmloop01);
    //             	break;
    //             case 2:
    //             	insert_node(i,l,j,best_d,P_PMmloop00);
    //             	insert_node(k,best_d-1,P_WBP);
    //             	break;
    //         }

	// 	}
	// 		break;
	// 	case P_PMmloop10:
	// 	{
	// 		cand_pos_t i = cur_interval->i;
	// 		cand_pos_t l = cur_interval->j;
	// 		cand_pos_t j = cur_interval->k;
	// 		cand_pos_t k = cur_interval->l;
	// 		if(debug) printf("At P_PMmloop10 with i=%d,l=%d,j=%d,k=%d\n",i,l,j,k);

	// 		if (!(i <= j && j < k-1 && k <= l)){
	// 			std::cerr << "border cases: This should not have happened!, P_PMmloop10" << std::endl;
	// 			exit(EXIT_FAILURE);
	// 		}

	// 		if (i<=0 || j<=0 || k<=0 || l<=0 || i>n || j>n || k> n || l> n){
	// 			std::cerr << "impossible cases: This should not have happened!, P_PMmloop10" << std::endl;
	// 			exit(EXIT_FAILURE);
	// 		}

	// 		energy_t temp=INF;
	// 		energy_t min_energy = PMmloop10.get(i,j-1,k,l)+cp_penalty;
	// 		cand_pos_t best_row = 1, best_d = -1;

    //         for(cand_pos_t d = i+1; d < j; ++d){
    //             temp = WBP.get(d,j) + PMmloop00.get(i,d-1,k,l);
    //             if (temp < min_energy){
    //                 min_energy = temp;
    //                 best_row = 2;
    //                 best_d = d;
    //             }
    //         }

    //         switch (best_row){
    //             case 1:
    //             	insert_node(i,l,j-1,k,P_PMmloop10);
    //             	break;
    //             case 2:
    //             	insert_node(i,l,best_d-1,k,P_PMmloop00);
    //             	insert_node(best_d,j,P_WBP);
    //             	break;
    //         }

	// 	}
	// 		break;
	// 	case P_POiloop:
	// 	{
	// 		cand_pos_t i = cur_interval->i;
	// 		cand_pos_t l = cur_interval->j;
	// 		cand_pos_t j = cur_interval->k;
	// 		cand_pos_t k = cur_interval->l;
	// 		if(debug) printf("At P_POiloop with i=%d,l=%d,j=%d,k=%d\n",i,l,j,k);

	// 		if (i<=0 || j<=0 || k<=0 || l<=0 || i>n || j>n || k> n || l> n){
	// 			std::cerr << "impossible cases: This should not have happened!, P_POiloop" << std::endl;
	// 			exit(EXIT_FAILURE);
	// 		}

	// 		if (!(i <= j && j < k-1 && k <= l)){
	// 			std::cerr << "border cases: This should not have happened!, P_POiloop" << std::endl;
	// 			exit(EXIT_FAILURE);
	// 		}

	// 		f[i].pair = l;
	// 		f[l].pair = i;
	// 		f[i].type = P_POiloop;
	// 		f[l].type = P_POiloop;


	// 		energy_t min_energy = INF;
	// 		cand_pos_t best_d=-1,best_dp=-1,best_row=-1;
	// 		if (pair[S_[i]][S_[l]]>0){
	// 			//branch 1
	// 			energy_t temp = PO.get(i+1,j,k,l-1) + get_e_stP(i,l);
	// 			if(temp < min_energy){
	// 				min_energy = temp;
	// 				best_row=1;
	// 			}
	// 			cand_pos_t max_d = std::min(j,i+MAXLOOP);
	// 			for(cand_pos_t d= i+1; d<max_d; ++d){
	// 				cand_pos_t min_dp = std::max(l-MAXLOOP,k);
	// 				for (cand_pos_t dp=l-1; dp >min_dp; --dp) {
	// 					energy_t branch2 = get_e_intP(i,d,dp,l) + PO.get(d,j,dp,k);
	// 					if(branch2 < min_energy){
	// 						min_energy = branch2;
	// 						best_row = 2;
	// 						best_d = d;
	// 						best_dp = dp;
	// 					}
	// 				}

    // 			}
	// 		}

	// 		switch (best_row)
	// 		{
	// 			case 1:
	// 				insert_node(i+1,l-1,j,k,P_PO);
	// 				break;
	// 			case 2:
	// 				insert_node(best_d,k,j,best_dp,P_PO);
	// 				break;
	// 		}

	// 	}
	// 		break;
	// 	case P_POmloop:
	// 	{
	// 		cand_pos_t i = cur_interval->i;
	// 		cand_pos_t l = cur_interval->j;
	// 		cand_pos_t j = cur_interval->k;
	// 		cand_pos_t k = cur_interval->l;
	// 		if(debug) printf("At P_POmloop with i=%d,l=%d,j=%d,k=%d\n",i,l,j,k);

	// 		if (!(i <= j && j < k-1 && k <= l)){
	// 			std::cerr << "border cases: This should not have happened!, P_POmloop" << std::endl;
	// 			exit(EXIT_FAILURE);
	// 		}

	// 		if (i<=0 || j<=0 || k<=0 || l<=0 || i>n || j>n || k> n || l> n){
	// 			std::cerr << "impossible cases: This should not have happened!, P_POmloop" << std::endl;
	// 			exit(EXIT_FAILURE);
	// 		}
	// 		f[i].pair = l;
	// 		f[l].pair = i;
	// 		f[i].type = P_POmloop;
	// 		f[l].type = P_POmloop;

    //         energy_t branch1 = POmloop10.get(i+1,j,k,l-1) + ap_penalty + beta2P(l,i);
    //         energy_t branch2 = POmloop01.get(i+1,j,k,l-1) + ap_penalty + beta2P(l,i);

    //         cand_pos_t best_row = -1;
    //         if (branch1 < branch2) { best_row = 1;} else { best_row = 2;}

	// 		switch (best_row){
	// 			case 1:
	// 				insert_node(i+1,l-1,j,k,P_POmloop10);
	// 				break;
	// 			case 2:
	// 				insert_node(i+1,l-1,j,k,P_POmloop01);
	// 				break;
	// 		}
	// 	}
	// 	break;

	// 	case P_POmloop00:
	// 	{
	// 		cand_pos_t i = cur_interval->i;
	// 		cand_pos_t l = cur_interval->j;
	// 		cand_pos_t j = cur_interval->k;
	// 		cand_pos_t k = cur_interval->l;
	// 		if(debug) printf("At P_POmloop00 with i=%d,l=%d,j=%d,k=%d\n",i,l,j,k);

	// 		if (!(i <= j && j < k-1 && k <= l)){
	// 			std::cerr << "border cases: This should not have happened!, P_POmloop00" << std::endl;
	// 			exit(EXIT_FAILURE);
	// 		}

	// 		if (i<=0 || j<=0 || k<=0 || l<=0 || i>n || j>n || k> n || l> n){
	// 			std::cerr << "impossible cases: This should not have happened!, P_POmloop00" << std::endl;
	// 			exit(EXIT_FAILURE);
	// 		}

	// 		energy_t min_energy = PO.get(i,j,k,l)+beta2P(l,i),temp=INF;
	// 		cand_pos_t best_row=1,best_d=-1;

	// 		for(cand_pos_t d=i+1; d<=j; ++d){
	// 			temp=get_WB(i,d-1)+POmloop00.get(d,j,k,l);
	// 			if (temp < min_energy){
	// 				min_energy = temp;
	// 				best_row = 2;
	// 				best_d = d;
	// 			}
	// 		}
	// 		for(cand_pos_t d=k; d<l; ++d){
	// 			temp=POmloop00.get(i,j,k,d)+get_WB(d+1,l);
	// 			if (temp < min_energy){
	// 				min_energy = temp;
	// 				best_row = 3;
	// 				best_d = d;
	// 			}
	// 		}

	// 		switch  (best_row){
	// 			case 1:
	// 				insert_node(i,l,j,k,P_PO);
	// 				break;
	// 			case 2:
	// 				insert_node(best_d,l,j,k,P_POmloop00);
	// 				insert_node(i,best_d-1,P_WBP);
	// 				break;
	// 			case 3:
	// 				insert_node(i,best_d,j,k,P_POmloop00);
	// 				insert_node(best_d+1,l,P_WB);
	// 				break;
	// 		}
	// 	}
	// 		break;

	// 	case P_POmloop01:
	// 	{
	// 		cand_pos_t i = cur_interval->i;
	// 		cand_pos_t l = cur_interval->j;
	// 		cand_pos_t j = cur_interval->k;
	// 		cand_pos_t k = cur_interval->l;
	// 		if(debug) printf("At P_POmloop01 with i=%d,l=%d,j=%d,k=%d\n",i,l,j,k);

	// 		if (!(i <= j && j < k-1 && k <= l)){
	// 			std::cerr << "border cases: This should not have happened!, P_POmloop01" << std::endl;
	// 			exit(EXIT_FAILURE);
	// 		}

	// 		if (i<=0 || j<=0 || k<=0 || l<=0 || i>n || j>n || k> n || l> n){
	// 			std::cerr << "impossible cases: This should not have happened!, P_POmloop01" << std::endl;
	// 			exit(EXIT_FAILURE);
	// 		}

	// 		energy_t min_energy = INF,temp=INF;
	// 		cand_pos_t best_d=-1;
    //         for(cand_pos_t d = k; d < l; d++){
    //             temp = POmloop00.get(i,j,k,d) + WBP.get(d+1,l);
    //             if (temp < min_energy){
    //                 min_energy = temp;
	// 				best_d = d;
    //             }
    //         }
	// 		insert_node(i,best_d,j,k,P_POmloop00);
    //         insert_node(best_d+1,l,P_WBP);
	// 	}
	// 		break;
	// 	case P_POmloop10:
	// 	{
	// 		cand_pos_t i = cur_interval->i;
	// 		cand_pos_t l = cur_interval->j;
	// 		cand_pos_t j = cur_interval->k;
	// 		cand_pos_t k = cur_interval->l;
	// 		if(debug) printf("At P_POmloop10 with i=%d,l=%d,j=%d,k=%d\n",i,l,j,k);

	// 		if (!(i <= j && j < k-1 && k <= l)){
	// 			std::cerr << "border cases: This should not have happened!, P_POmloop10" << std::endl;
	// 			exit(EXIT_FAILURE);
	// 		}

	// 		if (i<=0 || j<=0 || k<=0 || l<=0 || i>n || j>n || k> n || l> n){
	// 			std::cerr << "impossible cases: This should not have happened!, P_POmloop10" << std::endl;
	// 			exit(EXIT_FAILURE);
	// 		}

	// 		energy_t min_energy = INF,temp=INF;
	// 		cand_pos_t best_d=-1,best_row=-1;
	// 		for(cand_pos_t d=i+1; d<=j;++d){
    //             temp=WBP.get(i,d-1)+POmloop00.get(d,j,k,l);
    //             if (temp < min_energy){
    //                 min_energy = temp;
    //                 best_row = 1;
    //                 best_d = d;
    //             }

    //         }
    //         for(cand_pos_t d = k+1; d < l; ++d){
    //             temp = POmloop10.get(i,j,k,d) + get_WB(d+1,l);
    //             if (temp < min_energy){
    //                 min_energy = temp;
    //                 best_row = 2;
    //                 best_d = d;
    //             }
    //         }

	// 		switch  (best_row){
    //             case 1:
    //                 insert_node(best_d,l,j,k,P_POmloop00);
    //                 insert_node(i,best_d-1,P_WBP);
    //                 break;
    //             case 2:
    //                 insert_node(i,best_d,j,k,P_POmloop10);
    //                 insert_node(best_d+1,l,P_WB);
    //                 break;
	// 		}
	// 	}
	// 		break;
	// }
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
    tmp->type = type;
    tmp->next = stack_interval;
    stack_interval = tmp;
}

void pseudo_loop::set_stack_interval(seq_interval *stack_interval){
	this->stack_interval = stack_interval;
}
