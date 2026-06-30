#include "pseudo_loop.hh"
#include "h_globals.hh"
#include <stdio.h>
#include <string>
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <algorithm>
#include <cassert>

pseudo_loop::pseudo_loop(std::string seq, int dangle) : seq(seq), params_(vrna_params(NULL))
{
	n = seq.length();
	params_->model_details.dangles = dangle;
    S_ = encode_sequence(seq.c_str(),0);
	S1_ = encode_sequence(seq.c_str(),1);
    allocate_space();
}

void pseudo_loop::allocate_space()
{

	TriangleMatrix::new_index(index,n+1);
	Matrix4D::construct_index(index3D,n);

	fres.resize(n+1,-2);

	W.resize(n+1,0);
	cand_pos_t total_length = ((n+1) *(n+2))/2;
	V.resize(total_length);
	WM.init(n+1,index);
	WMv.init(n+1,index);
	WMp.init(n+1,index);
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
}

pseudo_loop::~pseudo_loop()
{
	free(params_);
	free(S_);
	free(S1_);
}

double pseudo_loop::ccj (){
	for (cand_pos_t i = n; i>=1; --i){	
		for (cand_pos_t j =i; j<=n; ++j){
			// compute_energy (i,j);
			compute_energies(i,j);
			// compute_WMv_WMp(i,j);
			// compute_energy_WM(i,j);
		}
	}
	for (cand_pos_t j= TURN+1; j <= n; j++){
		energy_t m1 = INF, m2 = INF, m3 = INF;
		m1 = W[j-1];
		for (cand_pos_t k=1; k<=j-TURN-1; ++k){
			energy_t acc = (k>1) ? W[k-1]: 0;
			m2 = std::min(m2,acc + E_ext_Stem(get_energy(k,j),get_energy(k+1,j),get_energy(k,j-1),get_energy(k+1,j-1),S_,params_,k,j,n));
			m3 = std::min(m3,acc + P.get(k,j) + PS_penalty);
		}
		W[j] = std::min({m1,m2,m3});
	}

    double energy = W[n]/100.0;

	// backtrack
	backtrack();

	fill_structure(fres,structure);
	this->structure = structure.substr(1,n);
    return energy;
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

void pseudo_loop::compute_WBP(cand_pos_t i, cand_pos_t l){
	energy_t min_energy= INF, b1 = INF, b2=INF, b3 = INF,tmp =INF;

	for(cand_pos_t d=i; d< l; ++d){
		tmp = calc_WB(i,d-1) + get_energy(d,l) + bp_penalty + PPS_penalty;
		b1 = std::min(b1,tmp);
		tmp = calc_WB(i,d-1) + P.get(d,l) + PSM_penalty + PPS_penalty;
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
		tmp = calc_WP(i,d-1) + get_energy(d,l) + gamma2(l,d) + PPS_penalty;
		b1 = std::min(b1,tmp);
		tmp = calc_WP(i,d-1) + P.get(d,l) + PSP_penalty + PPS_penalty;
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
				b20 = PK1LMorO.get(i,j,d+1,k) + PK2MreR.get(j+1,d,k+1,l);
				b21 = PK1LMorO.get(i,j,d+1,k) + PK2R.get(j+1,d,k+1,l);
				min_energy = std::min({min_energy,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,b15,b16,b17,b18,b19,b20,b21});
			}
		}
	}
	if (min_energy < INF/2){
		P.set(i,l) = min_energy;
	}
}
/**
 * This always chops from the interior with every single recurrence which means the code is the same and can be made generic
 */
void pseudo_loop::compute_PK1X(const Index4D &x, MType type){
	energy_t min_energy = INF;
	const cand_pos_t i = x.i(), j = x.j(), k = x.k(), l = x.l();
	Matrix4D &PX = PX_by_mtype(type);

	for(cand_pos_t d = i+1;d<=j;++d){
		min_energy = std::min(min_energy,PX.get(i,d,k,l) + calc_WP(d+1,j) + PB_penalty);
	}
	Matrix4D &PK1X = PK1X_by_mtype(type);
	PK1X.set(i,j,k,l,min_energy);
}
/**
 * Same as PK1
 */
void pseudo_loop::compute_PK2X(const Index4D &x, MType type){
	energy_t min_energy = INF;
	const cand_pos_t i = x.i(), j = x.j(), k = x.k(), l = x.l();
	Matrix4D &PK1X = PK1X_by_mtype(type);

	for(cand_pos_t d = k;d<l;++d){
		min_energy = std::min(min_energy,PK1X.get(i,j,d,l) + calc_WP(k,d-1));
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
    UNREACHABLE();
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
    UNREACHABLE();
}
void pseudo_loop::compute_PfromXprime(const Index4D &x, MType type){
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

energy_t pseudo_loop::calc_PfromXdoubleprime(cand_pos_t i,cand_pos_t j, cand_pos_t k, cand_pos_t l, MType type){
    switch(type) {
    case MType::L: return calc_PfromLdoubleprime(i,j,k,l);
    case MType::M: return calc_PfromMdoubleprime(i,j,k,l);
    case MType::R: return calc_PfromRdoubleprime(i,j,k,l);
    case MType::Om: return calc_PfromOdoubleprime(i,j,k,l);
	case MType::Os: return INF;
	case MType::LreO: return calc_PfromLreOdoubleprime(i,j,k,l);
	case MType::LreR: return calc_PfromLreRdoubleprime(i,j,k,l);
	case MType::MreO: return calc_PfromMreOdoubleprime(i,j,k,l);
	case MType::MreR: return calc_PfromMreRdoubleprime(i,j,k,l);
	case MType::LMreR: return calc_PfromLMreRdoubleprime(i,j,k,l);
	case MType::LMorO: return calc_PfromLMorOdoubleprime(i,j,k,l);
    }
    UNREACHABLE();
}
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
    UNREACHABLE();
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
    UNREACHABLE();
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
    UNREACHABLE();
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
///////////////// Traceback ////////////////////////////////

void pseudo_loop::backtrack(){
   Trace_W(1,n,W[n]);
   return;
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
	UNREACHABLE();
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
	UNREACHABLE();
}

void pseudo_loop::Trace_PXmloop(const Index4D &x, MType type, energy_t e){
	if (debug) std::cout << "PXmloop at " << x.i() << " and " << x.j() << " and " << x.k() << " and " << x.l() << " with type: " << type << " and en: " << e << std::endl;
	assert(impossible_case(x));

	Index4D xp(x);
	xp.shrink(type);
	Matrix4D &PXmloop00 = PXmloop00_by_mtype(type);
	energy_t tmp = PXmloop00.get(xp.i(),xp.j(),xp.k(),xp.l())+ ap_penalty + bp_penalty;
	if(e==tmp){
		return Trace_PXmloop00(xp,type,PXmloop00.get(xp.i(),xp.j(),xp.k(),xp.l()));
	}
	UNREACHABLE();
}

void pseudo_loop::Trace_PXiloop(const Index4D &x, MType type, energy_t e){
	if (debug) std::cout << "PXiloop at " << x.i() << " and " << x.j() << " and " << x.k() << " and " << x.l() << " with type: " << type << " and en: " << e << std::endl;
	switch(type) {
    case MType::L: return Trace_PLiloop(x,type,e);
    case MType::M: return Trace_PMiloop(x,type,e);
    case MType::R: return Trace_PRiloop(x,type,e);
    case MType::Om: return Trace_POiloop(x,type,e);
	case MType::Os: return Trace_POiloop(x,type,e);
	case MType::LreO: return Trace_PLiloop(x,type,e);
	case MType::LreR: return Trace_PLiloop(x,type,e);
	case MType::MreO: return Trace_PMiloop(x,type,e);
	case MType::MreR: return Trace_PMiloop(x,type,e);
	case MType::LMreR: return Trace_PLiloop(x,type,e);
	case MType::LMorO: return Trace_PLiloop(x,type,e);
    }
    UNREACHABLE();
}

void pseudo_loop::Trace_PfromX(const Index4D &x, MType type,energy_t e){
    switch(type) {
    case MType::L: return Trace_PfromL(x.i(),x.j(),x.k(),x.l(),type,e);
    case MType::M: return Trace_PfromM(x.i(),x.j(),x.k(),x.l(),type,e);
    case MType::R: return Trace_PfromR(x.i(),x.j(),x.k(),x.l(),type,e);
    case MType::Om: return Trace_PfromO(x.i(),x.j(),x.k(),x.l(),type,e);
	case MType::Os: return;
	case MType::LreO: return Trace_PfromL(x.i(),x.j(),x.k(),x.l(),type,e);
	case MType::LreR: return Trace_PfromL(x.i(),x.j(),x.k(),x.l(),type,e);
	case MType::MreO: return Trace_PfromM(x.i(),x.j(),x.k(),x.l(),type,e);
	case MType::MreR: return Trace_PfromM(x.i(),x.j(),x.k(),x.l(),type,e);
	case MType::LMreR: return Trace_PfromL(x.i(),x.j(),x.k(),x.l(),type,e);
	case MType::LMorO: return Trace_PfromL(x.i(),x.j(),x.k(),x.l(),type,e);
    }
    UNREACHABLE();
}

void pseudo_loop::Trace_PfromXprime(cand_pos_t i,cand_pos_t j,cand_pos_t k, cand_pos_t l, MType type,energy_t e){
    switch(type) {
    case MType::L: return Trace_PfromLprime(i,j,k,l,type,e);
    case MType::M: return Trace_PfromMprime(i,j,k,l,type,e);
    case MType::R: return Trace_PfromRprime(i,j,k,l,type,e);
    case MType::Om: return Trace_PfromOprime(i,j,k,l,type,e);
	case MType::Os: return;
	case MType::LreO: return Trace_PfromLprime(i,j,k,l,type,e);
	case MType::LreR: return Trace_PfromLprime(i,j,k,l,type,e);
	case MType::MreO: return Trace_PfromMprime(i,j,k,l,type,e);
	case MType::MreR: return Trace_PfromMprime(i,j,k,l,type,e);
	case MType::LMreR: return Trace_PfromLprime(i,j,k,l,type,e);
	case MType::LMorO: return Trace_PfromLprime(i,j,k,l,type,e);
    }
    UNREACHABLE();
}

void pseudo_loop::Trace_PfromXdoubleprime(cand_pos_t i,cand_pos_t j,cand_pos_t k, cand_pos_t l, MType type, energy_t e){
    switch(type) {
    case MType::L: return Trace_PfromLdoubleprime(i,j,k,l,e);
    case MType::M: return Trace_PfromMdoubleprime(i,j,k,l,e);
    case MType::R: return Trace_PfromRdoubleprime(i,j,k,l,e);
    case MType::Om: return Trace_PfromOdoubleprime(i,j,k,l,e);
	case MType::Os: return;
	case MType::LreO: return Trace_PfromLreOdoubleprime(i,j,k,l,e);
	case MType::LreR: return Trace_PfromLreRdoubleprime(i,j,k,l,e);
	case MType::MreO: return Trace_PfromMreOdoubleprime(i,j,k,l,e);
	case MType::MreR: return Trace_PfromMreRdoubleprime(i,j,k,l,e);
	case MType::LMreR: return Trace_PfromLMreRdoubleprime(i,j,k,l,e);
	case MType::LMorO: return Trace_PfromLMorOdoubleprime(i,j,k,l,e);
    }
    UNREACHABLE();
}

void pseudo_loop::Trace_PXmloop00(const Index4D &x, MType type, energy_t e){
	if (debug) std::cout << "PXmloop00 at " << x.i() << " and " << x.j() << " and " << x.k() << " and " << x.l() << " with type: " << type << " and en: " << e << std::endl;
    switch(type) {
    case MType::L: return Trace_PLmloop00(x,type,e);
    case MType::M: return Trace_PMmloop00(x,type,e);
    case MType::R: return Trace_PRmloop00(x,type,e);
    case MType::Om: return Trace_POmloop00(x,type,e);
	case MType::Os: return Trace_POmloop00(x,type,e);
	case MType::LreO: return Trace_PLmloop00(x,type,e);
	case MType::LreR: return Trace_PLmloop00(x,type,e);
	case MType::MreO: return Trace_PMmloop00(x,type,e);
	case MType::MreR: return Trace_PMmloop00(x,type,e);
	case MType::LMreR: return Trace_PLmloop00(x,type,e);
	case MType::LMorO: return Trace_PLmloop00(x,type,e);
    }
    UNREACHABLE();
}
void pseudo_loop::Trace_PXmloop10(const Index4D &x, MType type, energy_t e){
	if (debug) std::cout << "PXmloop10 at " << x.i() << " and " << x.j() << " and " << x.k() << " and " << x.l() << " with type: " << type << " and en: " << e << std::endl;
    switch(type) {
    case MType::L: return Trace_PLmloop10(x,type,e);
    case MType::M: return Trace_PMmloop10(x,type,e);
    case MType::R: return Trace_PRmloop10(x,type,e);
    case MType::Om: return Trace_POmloop10(x,type,e);
	case MType::Os: return Trace_POmloop10(x,type,e);
	case MType::LreO: return Trace_PLmloop10(x,type,e);
	case MType::LreR: return Trace_PLmloop10(x,type,e);
	case MType::MreO: return Trace_PMmloop10(x,type,e);
	case MType::MreR: return Trace_PMmloop10(x,type,e);
	case MType::LMreR: return Trace_PLmloop10(x,type,e);
	case MType::LMorO: return Trace_PLmloop10(x,type,e);
    }
    UNREACHABLE();
}
void pseudo_loop::Trace_PXmloop01(const Index4D &x, MType type, energy_t e){
	if (debug) std::cout << "PXmloop01 at " << x.i() << " and " << x.j() << " and " << x.k() << " and " << x.l() << " with type: " << type << " and en: " << e << std::endl;
    switch(type) {
    case MType::L: return Trace_PLmloop01(x,type,e);
    case MType::M: return Trace_PMmloop01(x,type,e);
    case MType::R: return Trace_PRmloop01(x,type,e);
    case MType::Om: return Trace_POmloop01(x,type,e);
	case MType::Os: return Trace_POmloop01(x,type,e);
	case MType::LreO: return Trace_PLmloop01(x,type,e);
	case MType::LreR: return Trace_PLmloop01(x,type,e);
	case MType::MreO: return Trace_PMmloop01(x,type,e);
	case MType::MreR: return Trace_PMmloop01(x,type,e);
	case MType::LMreR: return Trace_PLmloop01(x,type,e);
	case MType::LMorO: return Trace_PLmloop01(x,type,e);
    }
    UNREACHABLE();
}
void pseudo_loop::Trace_PX(cand_pos_t i,cand_pos_t j,cand_pos_t k, cand_pos_t l, MType type, energy_t e){
	if (debug) std::cout << "PX at " << i << " and " << j << " and " << k << " and " << l << " with type: " << type << " and en: " << e << std::endl;
	const Index4D x(i,j,k,l);
	const int ptype_closing = pair[S_[x.lend(type)]][S_[x.rend(type)]];
	fres[x.lend(type)] = x.rend(type);
	fres[x.rend(type)] = x.lend(type);

	if (ptype_closing>0){
		energy_t tmp = calc_PXmloop(x,type);
		if(e==tmp){
			Trace_PXmloop(x,type,tmp);
			return;
		}
		if(type == MType::Os){
			if(x.i()==x.j() && x.k()==x.l()){
				if(e==gamma2(x.l(),x.i())) return;
			}
		} else if (x.difference(type)>TURN){
			Index4D xp(x);
			xp.shrink(type);
			Matrix4D &PfromX = PfromX_by_mtype(type);
			tmp = PfromX.get(xp) + penalty(xp, gamma2, type);
			if(e==tmp){
				Trace_PfromX(xp,type,PfromX.get(xp));
				return;
			}
		}
		tmp = calc_PXiloop(x, type);
		if(e==tmp){
			Trace_PXiloop(x,type,tmp);
			return;
		}
	}
	UNREACHABLE();
}







////////////////// Util Functions ///////////////////////////

/**
 * @brief Gives the W(i,j) energy. The type of dangle model being used affects this energy. 
 * The type of dangle is also changed to reflect this.
 * 
 * @param vij The V(i,j) energy
 * @param vi1j The V(i+1,j) energy
 * @param vij1 The V(i,j-1) energy
 * @param vi1j1 The V(i+1,j-1) energy
*/
energy_t pseudo_loop::E_ext_Stem(const energy_t& vij,const energy_t& vi1j,const energy_t& vij1,const energy_t& vi1j1,const short* S, vrna_param_t* params, const cand_pos_t i,const cand_pos_t j, cand_pos_t n){

	energy_t e = INF, en = INF;
  	pair_type tt  = pair[S[i]][S[j]];
	
	en = vij; // i j

	if (en != INF) {
		if (params->model_details.dangles == 2){
			base_type si1 = i>1 ? S[i-1] : -1;
			base_type sj1 = j<n ? S[j+1] : -1;
			en += E_ExtLoop(tt, si1, sj1, params);
		}
		else{
			en += E_ExtLoop(tt, -1, -1, params);
		}
		e = std::min(e, en);
	}
	if(params->model_details.dangles  == 1){
        tt  = pair[S[i+1]][S[j]];
		en = (j-i-1>TURN) ? vi1j : INF; //i+1 j

		if (en != INF) {
			base_type si1 = S[i];
			en += E_ExtLoop(tt, si1, -1, params);
		}
		e = std::min(e,en);
        
        tt  = pair[S[i]][S[j-1]];
		en = (j-1-i>TURN) ? vij1 : INF; // i j-1
		if (en != INF) {
			base_type sj1 = S[j];
			en += E_ExtLoop(tt, -1, sj1, params);
		}
		e = std::min(e,en);

        tt  = pair[S[i+1]][S[j-1]];
		en = (j-1-i-1>TURN) ? vi1j1 : INF; // i+1 j-1
		if (en != INF) {
			base_type si1 = S[i];
			base_type sj1 = S[j];
			en += E_ExtLoop(tt, si1, sj1, params);
		}
		e = std::min(e,en);
	}
	return e;
}

/**
 * @brief Gives the WM(i,j) energy. The type of dangle model being used affects this energy. 
 * The type of dangle is also changed to reflect this.
 * 
 * I am adding +1 to all S as I haven't shifted the variables over to 1->n instead of 0->n-1
 * 
 * @param vij The V(i,j) energy
 * @param vi1j The V(i+1,j) energy
 * @param vij1 The V(i,j-1) energy
 * @param vi1j1 The V(i+1,j-1) energy
*/
energy_t pseudo_loop::E_MLStem(const energy_t& vij,const energy_t& vi1j,const energy_t& vij1,const energy_t& vi1j1,const short* S, vrna_param_t* params,cand_pos_t i, cand_pos_t j, cand_pos_t n){

	energy_t e = INF,en=INF;

	pair_type type = pair[S[i]][S[j]];

	
	en = vij; // i j
	if (en != INF) {
		if (params->model_details.dangles == 2){
			base_type mm5 = i>1 ? S[i-1] : -1;
			base_type mm3 = j<n ? S[j+1] : -1;
			en += E_MLstem(type, mm5, mm3, params);
		}
		else{
			en += E_MLstem(type, -1, -1, params);
		}
		e = std::min(e, en);
	}
	
	if(params->model_details.dangles == 1){
		const base_type mm5 = S[i], mm3 = S[j];


		en = (j-i-1 >TURN) ? vi1j : INF; // i+1 j
		if (en != INF) {
			en += params->MLbase;

			type = pair[S[i+1]][S[j]];
			en += E_MLstem(type, mm5, -1, params);

			e = std::min(e, en);
		}
    	

		en = (j-1-i>TURN) ? vij1 : INF; // i j-1
		if (en != INF) {
			en += params->MLbase;

			type = pair[S[i]][S[j-1]];
			en += E_MLstem(type, -1, mm3, params);

			e = std::min(e, en);
		}
    	
		en = (j-1-i-1>TURN) ? vi1j1 : INF; // i+1 j-1
		if (en != INF) {
			en += 2 * params->MLbase;

			type = pair[S[i+1]][S[j-1]];
			en += E_MLstem(type, mm5, mm3, params);
	
			e = std::min(e, en);
		}
	}


    return e;
}

/**
* @brief Computes the multiloop V contribution. This gives back essentially VM(i,j).
* 
* Added plus 1 to all S's as I haven't changed it over to 1->n from 0->n-1
* 
* @param dmli1 Row of WM2 from one iteration ago
* @param dmli2 Row of WM2 from two iterations ago 
*/
energy_t pseudo_loop::E_MbLoop(const energy_t WM2ij, const energy_t WM2ip1j, const energy_t WM2ijm1, const energy_t WM2ip1jm1, const short* S, vrna_param_t* params, cand_pos_t i, cand_pos_t j){

	energy_t e = INF,en = INF;
  	pair_type tt  = pair[S[j]][S[i]];
	
	/* double dangles */
	switch(params->model_details.dangles){
		case 2:
			e = WM2ij;

			if (e != INF) {

				base_type si1 = S[i+1];
				base_type sj1 = S[j-1];

				e += E_MLstem(tt, sj1, si1, params) + params->MLclosing;
			}

			break;

		case 1:
			/**
			* ML pair D0
			*  new closing pair (i,j) with mb part [i+1,j-1]  
			*/
			e = WM2ij;

			if (e != INF) {
				e += E_MLstem(tt, -1, -1, params) + params->MLclosing;
			}
      		/** 
			* ML pair 5
			* new closing pair (i,j) with mb part [i+2,j-1] 
			*/
			en = WM2ip1j;

			if (en != INF) {

				base_type si1 =  S[i+1];

				en += E_MLstem(tt, -1, si1, params) + params->MLclosing + params->MLbase;	
			}
      		e   = std::min(e, en);
			
			/** 
			* ML pair 3
			* new closing pair (i,j) with mb part [i+1, j-2] 
			*/
			en = WM2ijm1;

			if (en != INF) {
				base_type sj1 = S[j-1];

				en += E_MLstem(tt, sj1, -1, params) + params->MLclosing + params->MLbase;
			}
			e   = std::min(e, en);
			/** 
			* ML pair 53
			* new closing pair (i,j) with mb part [i+2.j-2]
			*/
			en = WM2ip1jm1;			

			if (en != INF) {

				base_type si1 = S[i+1];
				base_type sj1 = S[j-1];

				en += E_MLstem(tt, sj1, si1, params) + params->MLclosing + 2 * params->MLbase;
			}
			e = std::min(e, en);
      		break;
		case 0:
			e = WM2ij;

			if (e != INF) {
				e += E_MLstem(tt, -1, -1, params) + params->MLclosing;
			}
			break; 
	}
	return e;
}
void pseudo_loop::compute_WMv_WMp(cand_pos_t i, cand_pos_t j){
	if(j-i+1<4) return;

	WMv.set(i,j) = std::min(E_MLStem(get_energy(i,j),get_energy(i+1,j),get_energy(i,j-1),get_energy(i+1,j-1),S_,params_,i,j,n),WMv.get(i,j-1) + params_->MLbase);
	WMp.set(i,j) = std::min(P.get(i,j)+PSM_penalty+b_penalty,WMp.get(i,j-1) + params_->MLbase);
}

void pseudo_loop::compute_energy_WM (cand_pos_t i, cand_pos_t j)
// compute de MFE of a partial multi-loop closed at (i,j), the restricted case
{
    if(j-i+1<4) return;
	energy_t m1 = INF,m2=INF,m3=INF,m4=INF,m5=INF;
	
	for (cand_pos_t k=j-TURN-1; k >= i; --k)
	{
		energy_t wm_kj = E_MLStem(get_energy(k,j),get_energy(k+1,j),get_energy(k,j-1),get_energy(k+1,j-1),S_,params_,k,j,n);
		energy_t wmb_kj = P.get(k,j)+PSM_penalty+b_penalty;
		m1 = std::min(m1,static_cast<energy_t>((k-i)*params_->MLbase) + wm_kj);
		m2 = std::min(m2,static_cast<energy_t>((k-i)*params_->MLbase) + wmb_kj);
		m3 =  std::min(m3,WM.get(i,k-1) + wm_kj);
		m4 =  std::min(m4,WM.get(i,k-1) + wmb_kj);

	}
	m5 = std::min(m5,WM.get(i,j-1) + params_->MLbase);
	WM.set(i,j) = std::min({m1,m2,m3,m4,m5});  
}

energy_t pseudo_loop::compute_energy_VM(cand_pos_t i, cand_pos_t j)
// compute the MFE of a multi-loop closed at (i,j), the restricted case
{
    energy_t min = INF;
    for (cand_pos_t k = i+1; k <= j-3; ++k)
    {
        energy_t WM2ij = WM.get(i+1,k-1) + WMv.get(k,j-1);
		WM2ij = std::min(WM2ij,WM.get(i+1,k-1) + WMp.get(k,j-1));
		WM2ij = std::min(WM2ij,static_cast<energy_t>((k-i-1)*params_->MLbase) + WMp.get(k,j-1));

        energy_t WM2ip1j = WM.get(i+2,k-1) + WMv.get(k,j-1);
		WM2ip1j = std::min(WM2ip1j,WM.get(i+2,k-1) + WMp.get(k-1,j-1));
		WM2ip1j = std::min(WM2ip1j,static_cast<energy_t>((k-(i+1)-1)*params_->MLbase) + WMp.get(k,j-1));

        energy_t WM2ijm1 = WM.get(i+1,k-1) + WMv.get(k,j-2);
		WM2ijm1 = std::min(WM2ijm1, WM.get(i+1,k-1) + WMp.get(k,j-2));
		WM2ijm1 = std::min(WM2ijm1,static_cast<energy_t>((k-i-1)*params_->MLbase) + WMp.get(k,j-2));

        energy_t WM2ip1jm1 = WM.get(i+2,k-1) + WMv.get(k,j-2);
		WM2ip1jm1 = std::min(WM2ip1jm1,WM.get(i+2,k-1) + WMp.get(k,j-2));
		WM2ip1jm1 = std::min(WM2ip1jm1,static_cast<energy_t>((k-(i+1)-1)*params_->MLbase) + WMp.get(k,j-2));

        min = std::min(min,E_MbLoop(WM2ij,WM2ip1j,WM2ijm1,WM2ip1jm1,S_,params_,i,j));
    }
    return min;
}

/**
 * @brief This code returns the hairpin energy for a given base pair.
 * @param i The left index in the base pair
 * @param j The right index in the base pair
*/
energy_t pseudo_loop::HairpinE(const std::string& seq, cand_pos_t i, cand_pos_t j) {
	
	const int ptype_closing = pair[S_[i]][S_[j]];

	if (ptype_closing==0) return INF;

	return E_Hairpin(j-i-1,ptype_closing,S1_[i+1],S1_[j-1],&seq.c_str()[i-1], const_cast<vrna_param_t*>(params_));
}

/**
 * @brief non-restricted version
*/
energy_t pseudo_loop::compute_internal(cand_pos_t i, cand_pos_t j){
	energy_t v_iloop = INF;
	cand_pos_t max_k = std::min(j-TURN-2,i+MAXLOOP+1);
	const int ptype_closing = pair[S_[i]][S_[j]];
	for ( cand_pos_t k=i+1; k<=max_k; ++k) {
		cand_pos_t min_l=std::max(k+TURN+1 + MAXLOOP+2, k+j-i) - MAXLOOP-2;
		for (int l=j-1; l>=min_l; --l) {
			energy_t v_iloop_kl = E_IntLoop(k-i-1,j-l-1,ptype_closing,rtype[pair[S_[k]][S_[l]]],S1_[i+1],S1_[j-1],S1_[k-1],S1_[l+1],const_cast<vrna_param_t*>(params_)) + get_energy(k,l);
			v_iloop = std::min(v_iloop,v_iloop_kl);	
		} 
	}
	return v_iloop;
}

energy_t pseudo_loop::compute_stack(cand_pos_t i, cand_pos_t j){

	const int ptype_closing = pair[S_[i]][S_[j]];
	cand_pos_t k = i+1;
    cand_pos_t l = j-1;
    return E_IntLoop(k-i-1,j-l-1,ptype_closing,rtype[pair[S_[k]][S_[l]]],S1_[i+1],S1_[j-1],S1_[k-1],S1_[l+1],const_cast<vrna_param_t*>(params_)) + get_energy(k,l);
}

void pseudo_loop::compute_energy (cand_pos_t i, cand_pos_t j)
// compute the V(i,j) value, if the structure must be restricted
{
    energy_t min, min_en[3];
    cand_pos_t k, min_rank;
    char type;

    min_rank = -1;
    min = INF/2;
    min_en[0] = INF;
    min_en[1] = INF;
    min_en[2] = INF;



	min_en[0] = HairpinE(seq,i,j);

	min_en[1] = compute_internal(i,j);
	min_en[2] = compute_energy_VM(i,j);
    

    for (k=0; k<3; k++)
    {
        if (min_en[k] < min)
        {
            min = min_en[k];
            min_rank = k;
        }
    }

    switch (min_rank)
    {
        case  0: type = HAIRP; break;
        case  1: type = INTER; break;
        case  2: type = MULTI; break;
        default: type = NONE;
    }

    if (min < INF/2) {
        int ij = index[i]+j-i;
        V[ij].energy = min;
        V[ij].type = type;
    }
}

energy_t pseudo_loop::compute_int(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){

	const pair_type ptype_closing = pair[S_[i]][S_[j]];
    return E_IntLoop(k-i-1,j-l-1,ptype_closing,rtype[pair[S_[k]][S_[l]]],S1_[i+1],S1_[j-1],S1_[k-1],S1_[l+1],const_cast<vrna_param_t*>(params_));
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