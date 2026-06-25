#include "W_final.hh"
#include "h_struct.hh"
#include "h_externs.hh"

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <algorithm>
#include <iostream>
#include <stack>
#include <list>
#define debug 0

// Hosna June 20th, 2007
// calls the constructor for s_min_folding
// to create all the matrixes required for simfold
// and then calls allocate_space in here to allocate
// space for WMB and V_final
W_final::W_final(std::string seq, int dangle) : params_(vrna_params(NULL))
{
	seq_ = seq;
	this->n = seq.length();
	make_pair_matrix();
	params_->model_details.dangles = dangle;
    S_ = encode_sequence(seq.c_str(),0);
	S1_ = encode_sequence(seq.c_str(),1);
	W.resize(n+1,0);
	space_allocation();
}


W_final::~W_final()
{
	delete P;
	delete V;
	delete [] f;
	free(params_);
	free(S_);
	free(S1_);
}

// Hosna June 20th, 2007
// allocates space for WMB object and V_final
void W_final::space_allocation(){

	// From simfold
	f = new minimum_fold [n+1];

    V = new s_energy_matrix (seq_, n,S_,S1_,params_);
	structure = std::string (n+1,'.');

	// Hosna: June 20th 2007
    P = new pseudo_loop (seq_,V,this,S_,S1_,params_);

}

double W_final::ccj(){

	for (cand_pos_t i = n; i>=1; --i){	
		for (cand_pos_t j =i; j<=n; ++j){
			V->compute_energy (i,j);
			P->compute_energies(i,j);
			V->compute_WMv_WMp(i,j,P->get_energy(i,j));
			V->compute_energy_WM(i,j,P->P);
		}
	}
	for (cand_pos_t j= TURN+1; j <= n; j++){
		energy_t m1 = INF, m2 = INF, m3 = INF;
		m1 = W[j-1];
		for (cand_pos_t k=1; k<=j-TURN-1; ++k){
			energy_t acc = (k>1) ? W[k-1]: 0;
			m2 = std::min(m2,acc + E_ext_Stem(V->get_energy(k,j),V->get_energy(k+1,j),V->get_energy(k,j-1),V->get_energy(k+1,j-1),S_,params_,k,j,n));
			m3 = std::min(m3,acc + P->get_energy(k,j) + PS_penalty);
		}
		W[j] = std::min({m1,m2,m3});
	}

    double energy = W[n]/100.0;

	P->set_fold(f);
	// backtrack
	backtrack();

	fill_structure();
	this->structure = structure.substr(1,n);
    return energy;
}

/**
 * @brief Gives the W(i,j) energy. The type of dangle model being used affects this energy. 
 * The type of dangle is also changed to reflect this.
 * 
 * @param vij The V(i,j) energy
 * @param vi1j The V(i+1,j) energy
 * @param vij1 The V(i,j-1) energy
 * @param vi1j1 The V(i+1,j-1) energy
*/
energy_t W_final::E_ext_Stem(const energy_t& vij,const energy_t& vi1j,const energy_t& vij1,const energy_t& vi1j1,const short* S, vrna_param_t* params, const cand_pos_t i,const cand_pos_t j, cand_pos_t n){

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

void W_final::backtrack(){
   Trace_W(1,n,W[n]);
   return;
}

void W_final::Trace_W(cand_pos_t i, cand_pos_t j, energy_t e){
	if (debug) printf("W at %d and %d with %d\n", i, j, e);
	if (j<=i) return;

	energy_t acc = INF;

	// this case is for j unpaired, so I have to check that.
	energy_t tmp = W[j-1];
	if (e==tmp){
		Trace_W(i,j-1,W[j-1]);
		return;
	}
	for (cand_pos_t i=1; i<=j-1; i++){
		acc = (i>1) ? W[i-1] : 0;
		base_type si1 = i>1 ? S_[i-1] : -1;
		base_type sj1 = j<n ? S_[j+1] : -1;
		tmp = acc + V->get_energy(i,j) + ((params_->model_details.dangles == 2) ? E_ExtLoop(pair[S_[i]][S_[j]],si1,sj1,params_) : E_ExtLoop(pair[S_[i]][S_[j]],-1,-1,params_));
		if(e==tmp){
			Trace_W(1,i-1,W[i-1]);
			Trace_V(i,j,V->get_energy(i,j));
			return;
		}
		if(params_->model_details.dangles ==1){
			tmp = acc + V->get_energy(i+1,j) + E_ExtLoop(pair[S_[i+1]][S_[j]],S_[i],-1,params_);
			if(e==tmp){
				Trace_W(1,i-1,W[i-1]);
				Trace_V(i+1,j,V->get_energy(i+1,j));
				return;
			}
			tmp = acc + V->get_energy(i,j-1) + E_ExtLoop(pair[S_[i]][S_[j-1]],-1,S_[j],params_);
			if(e==tmp){
				Trace_W(1,i-1,W[i-1]);
				Trace_V(i,j-1,V->get_energy(i,j-1));
				return;
			}
			tmp = acc + V->get_energy(i+1,j-1) + E_ExtLoop(pair[S_[i+1]][S_[j-1]],S_[i],S_[j],params_);
			if(e==tmp){
				Trace_W(1,i-1,W[i-1]);
				Trace_V(i+1,j-1,V->get_energy(i+1,j-1));
				return;
			}
		}
	}
	for (cand_pos_t i=1; i<=j-1; i++){
		acc = (i-1>0) ? W[i-1]: 0;
		tmp = acc + P->get_energy(i,j)+ PS_penalty;
		if(e==tmp){
			Trace_W(1,i-1,W[i-1]);
			P->Trace_P(i,j,P->get_energy(i,j));
			return;
		}
	}
	__builtin_unreachable();
}
void W_final::Trace_V(cand_pos_t i, cand_pos_t j, energy_t e){
	if (debug) printf("V at %d and %d with %d\n", i, j, e);
	f[i].pair = j;
	f[j].pair = i;
	char type = V->get_type (i,j);

	switch(type){
		case HAIRP:{
			return;
		}
		case INTER:{
			cand_pos_t max_k = std::min(j-TURN-2,i+MAXLOOP+1);
			for (cand_pos_t k = i+1; k <= max_k; ++k){
				cand_pos_t min_l=std::max(k+TURN+1 + MAXLOOP+2, k+j-i) - MAXLOOP-2;
				for (cand_pos_t l = j-1; l >= min_l; --l)
				{
					energy_t tmp = V->compute_int(i,j,k,l,params_);
					if (e == tmp)
					{
						Trace_V(k,l,V->get_energy(k,l));
						return;
					}
				}
				
			}
		}
		case MULTI: {
			energy_t tmp = INF;
			for (cand_pos_t k = i+1; k <= j-1; k++){
				tmp = V->get_energy_WM(i+1,k-1) + std::min(V->get_energy_WMv(k,j-1),V->get_energy_WMp(k,j-1)) + E_MLstem(pair[S_[j]][S_[i]],-1,-1,params_) + params_->MLclosing;
				if (e==tmp){
					Trace_WM(i+1,k-1,V->get_energy_WM(i+1,k-1));
					if(tmp == V->get_energy_WM(i+1,k-1) + V->get_energy_WMv(k,j-1)){
						Trace_WMv(k,j-1,V->get_energy_WMv(k, j-1));
					} else {
						Trace_WMp(k,j-1,V->get_energy_WMv(k, j-1));
					}
					return;
				}
				tmp = V->get_energy_WM(i+2,k-1) + std::min(V->get_energy_WMv(k,j-1),V->get_energy_WMp(k,j-1)) + E_MLstem(pair[S_[j]][S_[i]],-1,S_[i+1],params_) + params_->MLclosing + params_->MLbase;
				if (e==tmp)
				{
					Trace_WM(i+2,k-1,V->get_energy_WM (i+2,k-1));
					if(tmp == V->get_energy_WM(i+2,k-1) + V->get_energy_WMv(k,j-1)){
						Trace_WMv(k,j-1,V->get_energy_WMv(k, j-1));
					} else {
						Trace_WMp(k,j-1,V->get_energy_WMv(k, j-1));
					}
					return;
				}
				tmp = V->get_energy_WM (i+1,k-1) + std::min(V->get_energy_WMv(k,j-2),V->get_energy_WMp(k,j-2)) + E_MLstem(pair[S_[j]][S_[i]],S_[j-1],-1,params_) + params_->MLclosing + params_->MLbase;
				if (e==tmp)
				{
					Trace_WM(i+1,k-1,V->get_energy_WM (i+1,k-1));
					if(tmp == V->get_energy_WM(i+1,k-1) + V->get_energy_WMv(k,j-2)){
						Trace_WMv(k,j-2,V->get_energy_WMv(k, j-2));
					} else {
						Trace_WMp(k,j-2,V->get_energy_WMv(k, j-2));
					}
					return;
				}
				tmp = V->get_energy_WM (i+2,k-1) + std::min(V->get_energy_WMv(k,j-2),V->get_energy_WMp(k,j-2)) + E_MLstem(pair[S_[j]][S_[i]],S_[j-1],S_[i+1],params_) + params_->MLclosing + 2*params_->MLbase;
				if (e==tmp)
				{
					Trace_WM(i+2,k-1,V->get_energy_WM (i+2,k-1));
					if(tmp == V->get_energy_WM(i+2,k-1) + V->get_energy_WMv(k,j-2)){
						Trace_WMv(k,j-2,V->get_energy_WMv(k, j-2));
					} else {
						Trace_WMp(k,j-2,V->get_energy_WMv(k, j-2));
					}
					return;
				}
				tmp = static_cast<energy_t>((k-i-1)*params_->MLbase + V->get_energy_WMp(k,j-1))+ E_MLstem(pair[S_[j]][S_[i]],-1,-1,params_) + params_->MLclosing;
				if (e==tmp){
					Trace_WMp(k,j-1,V->get_energy_WMp(k,j-1));
					return;
				}

				if((k-(i+1)-1) >=0) tmp = static_cast<energy_t>((k-(i+1)-1)*params_->MLbase) + V->get_energy_WMp(k,j-1) + E_MLstem(pair[S_[j]][S_[i]],-1,S_[i+1],params_) + params_->MLclosing + params_->MLbase;
				if (e==tmp){
					Trace_WMp(k,j-1,V->get_energy_WMp(k, j-1));
					return;
				}
				tmp = static_cast<energy_t>((k-i-1)*params_->MLbase) + V->get_energy_WMp(k,j-2) + E_MLstem(pair[S_[j]][S_[i]],S_[j-1],-1,params_) + params_->MLclosing + params_->MLbase;
				if (e==tmp){
					Trace_WMp(k,j-2,V->get_energy_WMp(k, j-2));
					return;
				}
				
				if((k-(i+1)-1) >=0) tmp = static_cast<energy_t>((k-(i+1)-1)*params_->MLbase) + V->get_energy_WMp(k,j-2) + E_MLstem(pair[S_[j]][S_[i]],S_[j-1],S_[i+1],params_) + params_->MLclosing + 2*params_->MLbase;
				if (e==tmp){
					Trace_WMp(k,j-2,V->get_energy_WMp(k, j-2));
					return;
				}					
			}
		}
	}
	__builtin_unreachable();
}
void W_final::Trace_WM(cand_pos_t i, cand_pos_t j, energy_t e){
	if (debug) printf("WM at %d and %d with %d\n", i, j, e);
	energy_t tmp = INF;

	tmp = V->get_energy_WM(i,j-1)+params_->MLbase;
	if(e==tmp){
		Trace_WM(i,j-1,V->get_energy_WM(i,j-1));
		return;
	}
	for (cand_pos_t k=i; k <= j-TURN-1; k++){	
		tmp = static_cast<energy_t>((k-i)*params_->MLbase) + V->get_energy_WMv (k,j);
		if(e==tmp){
			Trace_WMv(k,j,V->get_energy_WMv(k,j));
			return;
		}
		tmp = static_cast<energy_t>((k-i)*params_->MLbase) + V->get_energy_WMp(k,j);
		if(e==tmp){
			Trace_WMp(k,j,V->get_energy_WMp(k,j));
			return;
		}
		tmp = V->get_energy_WM(i,k-1) + V->get_energy_WMv(k,j);
		if(e==tmp){
			Trace_WM(i,k-1,V->get_energy_WM(i, k-1));
			Trace_WMv(k,j,V->get_energy_WMv(k,j));
			return;
		}
		tmp = V->get_energy_WM(i,k-1) + V->get_energy_WMp(k,j);
		if(e==tmp){
			Trace_WM(i,k-1,V->get_energy_WM(i, k-1));
			Trace_WMp(k,j,V->get_energy_WMp(k,j));
			return;
		}
	}
	__builtin_unreachable();
}
void W_final::Trace_WMv(cand_pos_t i, cand_pos_t j, energy_t e){
	if (debug) printf("WMv at %d and %d with %d\n", i, j, e);
	cand_pos_t si = S_[i];
	cand_pos_t sj = S_[j];
	cand_pos_t si1 = (i>1) ? S_[i-1] : -1;
	cand_pos_t sj1 = (j<n) ? S_[j+1] : -1;
	pair_type tt = pair[S_[i]][S_[j]];
	energy_t tmp = V->get_energy(i,j) + ((params_->model_details.dangles == 2) ? E_MLstem(tt,si1,sj1,params_) : E_MLstem(tt,-1,-1,params_));
	if(e==tmp){
		Trace_V(i,j,V->get_energy(i,j));
		return;
	}

	if(params_->model_details.dangles == 1){
		tt = pair[S_[i+1]][S_[j]];
		energy_t tmp = V->get_energy(i+1,j) + E_MLstem(tt,si,-1,params_) + params_->MLbase;
		if(e==tmp){
			Trace_V(i+1,j,V->get_energy(i+1,j));
			return;
		}
		tt = pair[S_[i]][S_[j-1]];
		tmp = V->get_energy(i,j-1) + E_MLstem(tt,-1,sj,params_) + params_->MLbase;
		if(e==tmp){
			Trace_V(i,j-1,V->get_energy(i,j-1));
			return;
		}
		tt = pair[S_[i+1]][S_[j-1]];
		tmp = V->get_energy(i+1,j-1) + E_MLstem(tt,si,sj,params_) + 2*params_->MLbase;
		if(e==tmp){
			Trace_V(i+1,j-1,V->get_energy(i+1,j-1));
			return;
		}
	}

	tmp = V->get_energy_WMv(i,j-1) + params_->MLbase;
	if(e==tmp){
		Trace_WMv(i,j-1,V->get_energy_WMv(i,j-1));
		return;
	}
	__builtin_unreachable();
}
void W_final::Trace_WMp(cand_pos_t i, cand_pos_t j, energy_t e){
	if (debug) printf("WMp at %d and %d with %d\n", i, j, e);
	energy_t tmp = P->get_energy(i,j) + PSM_penalty + b_penalty;
	if(e==tmp){
		P->Trace_P(i,j,V->get_energy_WMp(i,j));
		return;
	}
	tmp = V->get_energy_WMp(i,j-1) + params_->MLbase;
	if(e==tmp){
		Trace_WMp(i,j-1,V->get_energy_WMp(i,j-1));
		return;
	}
	__builtin_unreachable();
}

void W_final::fill_structure()
{
    std::stack < brack_type > st;

    st.push(brack_type('<','>'));

    st.push(brack_type('{','}'));

    st.push(brack_type('[',']'));

    st.push(brack_type('(',')'));

    cand_pos_t isInABand=0;
    // cand_pos_t num_crossing_bands=0;

    std::list <band_elem > bands;
    bands.push_back(band_elem('|','|',0,0,0,0));

    for (cand_pos_t i = 1; i <= n; i++){
        cand_pos_t j = f[i].pair;
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