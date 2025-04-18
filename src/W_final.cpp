#include "W_final.hh"
#include "h_struct.hh"
#include "h_externs.hh"

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <algorithm>
#include <iostream>


// Hosna June 20th, 2007
// calls the constructor for s_min_folding
// to create all the matrixes required for simfold
// and then calls allocate_space in here to allocate
// space for WMB and V_final
W_final::W_final(std::string seq, int dangle) : params_(scale_parameters())
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
    P = new pseudo_loop (seq_,V,S_,S1_,params_);

}

double W_final::ccj(){

	for (int i = n; i >=1; --i){	
		for (int j =i; j<=n; ++j){
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
			m3 = std::min(m3,acc + std::min({P->get_energy(k,j),P->get_energy(k+1,j),P->get_energy(k,j-1),P->get_energy(k+1,j-1)}) + PS_penalty);
		}
		W[j] = std::min({m1,m2,m3});
	}

    double energy = W[n]/100.0;


	// backtrack
    // first add (1,n) on the stack
    stack_interval = new seq_interval;
    stack_interval->i = 1;
    stack_interval->j = n;
    stack_interval->energy = W[n];
    stack_interval->type = FREE;
    stack_interval->next = NULL;

    seq_interval *cur_interval = stack_interval;

    while ( cur_interval != NULL)
    {
        stack_interval = stack_interval->next;
        backtrack(cur_interval);
        delete cur_interval;    // this should make up for the new in the insert_node
        cur_interval = stack_interval;
    }
	fill_structure();
	this->structure = structure.substr(1,n);
    return energy;
}

/**
 * @brief Gives the W(i,j) energy. The type of dangle model being used affects this energy. 
 * The type of dangle is also changed to reflect this.
 * 
 * Until the changes to fres, I am adding +1 to the ptype closing and Si and Sj's to make them match - Mateo 2024
 * 
 * @param vij The V(i,j) energy
 * @param vi1j The V(i+1,j) energy
 * @param vij1 The V(i,j-1) energy
 * @param vi1j1 The V(i+1,j-1) energy
*/
energy_t W_final::E_ext_Stem(const energy_t& vij,const energy_t& vi1j,const energy_t& vij1,const energy_t& vi1j1,const short* S, paramT* params, const cand_pos_t i,const cand_pos_t j, cand_pos_t n){

	energy_t e = INF, en = INF;
  	pair_type tt  = pair[S[i]][S[j]];
	
	en = vij; // i j

	if (en != INF) {
		if (params->model_details.dangles == 2){
			base_type si1 = i>1 ? S[i-1] : -1;
			base_type sj1 = j<n ? S[j+1] : -1;
			en += vrna_E_ext_stem(tt, si1, sj1, params);
		}
		else{
			en += vrna_E_ext_stem(tt, -1, -1, params);
		}

		e = MIN2(e, en);
		
	}
	if(params->model_details.dangles  == 1){
        tt  = pair[S[i+1]][S[j]];
		en = (j-i-1>TURN) ? vi1j : INF; //i+1 j

		if (en != INF) {

			base_type si1 = S[i];
			en += vrna_E_ext_stem(tt, si1, -1, params);
		}

		e = MIN2(e,en);
        
        tt  = pair[S[i]][S[j-1]];
		en = (j-1-i>TURN) ? vij1 : INF; // i j-1
		if (en != INF) {

			base_type sj1 = S[j];

			en += vrna_E_ext_stem(tt, -1, sj1, params);
		}
		e = MIN2(e,en);

        tt  = pair[S[i+1]][S[j-1]];
		en = (j-1-i-1>TURN) ? vi1j1 : INF; // i+1 j-1

		if (en != INF) {

			base_type si1 = S[i];
			base_type sj1 = S[j];

			en += vrna_E_ext_stem(tt, si1, sj1, params);
		}
		e = MIN2(e,en);
	}
	return e;
}

void W_final::backtrack(seq_interval *cur_interval){
    char type;

	switch (cur_interval->type){
		case LOOP:
		{
			int i = cur_interval->i;
			int j = cur_interval->j;
			if (i >= j)
				return;
			f[i].pair = j;
			f[j].pair = i;

			structure[i] = '(';
			structure[j] = ')';		

			type = V->get_type (i,j);
			
			// Hosna, March 8, 2012
			// changing nested ifs to switch for optimality
			switch (type){
				case HAIRP:
				{
					f[i].type = HAIRP;
					f[j].type = HAIRP;
				}
					break;
				case INTER:
				{
					f[i].type = INTER;
					f[j].type = INTER;
					// detect the other closing pair
					cand_pos_t best_ip=j, best_jp=i;
					energy_t min = INF;
					cand_pos_t max_ip = std::min(j-TURN-2,i+MAXLOOP+1);
					for (cand_pos_t k = i+1; k <= max_ip; ++k)
					{
						cand_pos_t min_l=std::max(k+TURN+1 + MAXLOOP+2, k+j-i) - MAXLOOP-2;
						for (cand_pos_t l = j-1; l >= min_l; --l)
						{
							energy_t tmp = V->compute_int(i,j,k,l,params_);
							if (tmp < min)
							{
								min = tmp;
								best_ip = k;
								best_jp = l;
							}
						}
						
					}

					if (best_ip < best_jp)
						insert_node (best_ip, best_jp, LOOP);
					else
					{
						fprintf(stderr,"NOT GOOD RESTR INTER, i=%d, j=%d, best_ip=%d, best_jp=%d\n", i, j, best_ip, best_jp);
						exit (0);
					}
				}
					break;
				case MULTI:
				{
					f[i].type = MULTI;
					f[j].type = MULTI;
					cand_pos_t best_k = -1, best_row = -1;
					energy_t tmp= INF, min = INF;
					for (cand_pos_t k = i+1; k <= j-1; k++){
						tmp = V->get_energy_WM (i+1,k-1) + std::min(V->get_energy_WMv(k, j-1),V->get_energy_WMp(k, j-1)) + E_MLstem(pair[S_[j]][S_[i]],-1,-1,params_) + params_->MLclosing;
						if (tmp < min){
							min = tmp;
							best_k = k;
							best_row = 1;
						}


						tmp = V->get_energy_WM (i+2,k-1) + std::min(V->get_energy_WMv(k, j-1),V->get_energy_WMp(k, j-1)) + E_MLstem(pair[S_[j]][S_[i]],-1,S_[i+1],params_) + params_->MLclosing + params_->MLbase;
						if (tmp < min)
						{
							min = tmp;
							best_k = k;
							best_row = 2;
						}

						tmp = V->get_energy_WM (i+1,k-1) + std::min(V->get_energy_WMv(k, j-2),V->get_energy_WMp(k, j-2)) + E_MLstem(pair[S_[j]][S_[i]],S_[j-1],-1,params_) + params_->MLclosing + params_->MLbase;
						if (tmp < min)
						{
							min = tmp;
							best_k = k;
							best_row = 3;
						}

						tmp = V->get_energy_WM (i+2,k-1) + std::min(V->get_energy_WMv(k, j-2),V->get_energy_WMp(k, j-2)) + E_MLstem(pair[S_[j]][S_[i]],S_[j-1],S_[i+1],params_) + params_->MLclosing + 2*params_->MLbase;
						
						if (tmp < min)
						{
							min = tmp;
							best_k = k;
							best_row = 4;
						}
						

						tmp = static_cast<energy_t>((k-i-1)*params_->MLbase + V->get_energy_WMp(k,j-1))+ E_MLstem(pair[S_[j]][S_[i]],-1,-1,params_) + params_->MLclosing;
						if (tmp < min){
							min = tmp;
							best_k = k;
							best_row = 5;
						}

						if((k-(i+1)-1) >=0) tmp = static_cast<energy_t>((k-(i+1)-1)*params_->MLbase) + V->get_energy_WMp(k,j-1) + E_MLstem(pair[S_[j]][S_[i]],-1,S_[i+1],params_) + params_->MLclosing + params_->MLbase;
						if (tmp < min){
							min = tmp;
							best_k = k;
							best_row = 6;
						}
						tmp = static_cast<energy_t>((k-i-1)*params_->MLbase) + V->get_energy_WMp(k,j-2) + E_MLstem(pair[S_[j]][S_[i]],S_[j-1],-1,params_) + params_->MLclosing + params_->MLbase;
						if (tmp < min){
							min = tmp;
							best_k = k;
							best_row = 7;
						}
						
						if((k-(i+1)-1) >=0) tmp = static_cast<energy_t>((k-(i+1)-1)*params_->MLbase) + V->get_energy_WMp(k,j-2) + E_MLstem(pair[S_[j]][S_[i]],S_[j-1],S_[i+1],params_) + params_->MLclosing + 2*params_->MLbase;
						if (tmp < min){
							min = tmp;
							best_k = k;
							best_row = 8;
						}					
					  }
					switch (best_row)
					  {
					  case 1:
						insert_node (i+1, best_k-1, M_WM);
						insert_node (best_k, j-1, M_WM);
						break;
					  case 2:
						insert_node (i+2, best_k-1, M_WM);
						insert_node (best_k, j-1, M_WM);
						break;
					  case 3:
						insert_node (i+1, best_k-1, M_WM);
						insert_node (best_k, j-2, M_WM);
						break;
					  case 4:
						insert_node (i+2, best_k-1, M_WM);
						insert_node (best_k, j-2, M_WM);
						break;
					  case 5:
						insert_node (best_k, j-1, M_WM);
						break;
					  case 6:
						insert_node (best_k, j-1, M_WM);
						break;
					  case 7:
						insert_node (best_k, j-2, M_WM);
						break;
					  case 8:
						insert_node (best_k, j-2, M_WM);
						break;
					  }
				}
					break;
			}
		}
			break;
		case FREE:
		{
			cand_pos_t j = cur_interval->j;

			if (j==1) return;

			energy_t min = INF, tmp = INF, acc = INF, energy_ij = INF;
			cand_pos_t best_row = -1, best_i = -1;

			// this case is for j unpaired, so I have to check that.
			tmp = W[j-1];
			if (tmp < min){
				min = tmp;
				best_row = 0;
			}
			for (cand_pos_t i=1; i<=j-1; i++)    // no TURN
			{
				acc = (i>1) ? W[i-1] : 0;
				energy_ij = V->get_energy(i,j);

				if (energy_ij < INF)
				{	
					if(params_->model_details.dangles == 2){
						base_type si1 = i>1 ? S_[i-1] : -1;
						base_type sj1 = j<n ? S_[j+1] : -1;
						tmp = energy_ij + E_ExtLoop(pair[S_[i]][S_[j]],si1,sj1,params_) + acc;
					} else 
						tmp = energy_ij + E_ExtLoop(pair[S_[i]][S_[j]],-1,-1,params_) + acc; 
					if (tmp < min){
						min = tmp;
						best_i = i;
						best_row = 1;
					}
					
				}
				if(params_->model_details.dangles ==1){

					energy_ij = V->get_energy(i+1,j);
					if (energy_ij < INF){
						tmp = energy_ij + E_ExtLoop(pair[S_[i+1]][S_[j]],S_[i],-1,params_) + acc;
						
						if (tmp < min){
							min = tmp;
							best_i = i;
							best_row = 2;
						}
					}
					
					energy_ij = V->get_energy(i,j-1);
					if (energy_ij < INF){
						tmp = energy_ij + E_ExtLoop(pair[S_[i]][S_[j-1]],-1,S_[j],params_) + acc;
						if (tmp < min){
							min = tmp;
							best_i = i;
							best_row = 3;
						}
					}
					
					energy_ij = V->get_energy(i+1,j-1);
					if (energy_ij < INF){
						tmp = energy_ij + E_ExtLoop(pair[S_[i+1]][S_[j-1]],S_[i],S_[j],params_) + acc;
						if (tmp < min){
							min = tmp;
							best_i = i;
							best_row = 4;
						}
					}
				}
			}
			// Hosna June 30, 2007
			// The following would not take care of when
			// we have some unpaired bases before the start of the P
			for (cand_pos_t i=1; i<=j-1; i++){
				acc = (i-1>0) ? W[i-1]: 0;

				energy_ij = P->get_energy(i,j);

				if (energy_ij < INF)
				{
					tmp = energy_ij + PS_penalty + acc;

					if (tmp < min)
					{
						min = tmp;
						best_row = 5;
						best_i = i;
					}
				}

				if(params_->model_details.dangles ==1){
					energy_ij = P->get_energy(i+1,j);
					if (energy_ij < INF){
						tmp = energy_ij + PS_penalty + acc;
						if (tmp < min){
							min = tmp;
							best_row = 6;
							best_i = i;
						}
					}
				
					energy_ij = P->get_energy(i,j-1);
					if (energy_ij < INF){
						tmp = energy_ij + PS_penalty + acc;
						if (tmp < min){
							min = tmp;
							best_row = 7;
							best_i = i;
						}
					}
				
					energy_ij = P->get_energy(i+1,j-1);
					if (energy_ij < INF){
						tmp = energy_ij + PS_penalty + acc;
						if (tmp < min){
							min = tmp;
							best_row = 8;
							best_i = i;
						}
					}
				}
			}
			switch (best_row)
			{
				case 0:
					//printf("W(%d) case 0: inserting Free (0,%d)\n",j,j-1);
					insert_node (1, j-1, FREE); break;
				case 1:
					//printf("W(%d) case 1: inserting Loop(%d,%d) and Free (0,%d)\n",j,best_i,j,best_i-1);
					insert_node (best_i, j, LOOP);
					if (best_i-1 > 1)     // it was TURN instead of 0  - not sure if TURN shouldn't be here
						insert_node (1, best_i-1, FREE);
					break;
				case 2:
					//printf("W(%d) case 2: inserting Loop(%d,%d) and Free (0,%d)\n",j,best_i+1,j,best_i);
					insert_node (best_i+1, j, LOOP);
					if (best_i >= 1)// Hosna, March 26, 2012, was best_i-1 instead of best_i
						insert_node (1, best_i, FREE);
					break;
				case 3:
					//printf("W(%d) case 3: inserting Loop(%d,%d) and Free (0,%d)\n",j,best_i,j-1,best_i-1);
					insert_node (best_i, j-1, LOOP);
					if (best_i-1 > 1)
						insert_node (1, best_i-1, FREE);
					break;
				case 4:
					//printf("W(%d) case 4: inserting Loop(%d,%d) and Free (0,%d)\n",j,best_i+1,j-1,best_i);
					insert_node (best_i+1, j-1, LOOP);
					if (best_i >= 1) // Hosna, March 26, 2012, was best_i-1 instead of best_i
						insert_node (1, best_i, FREE);
					break;
				// Hosna: June 28, 2007
				// the last branch of W, which is P_i,j
				case 5:
					//printf("W(%d) case 5: inserting WMB(%d,%d) and Free (0,%d)\n",j,best_i,j,best_i-1);
					insert_node (best_i, j, P_P);
					if (best_i-1 > 1)     // it was TURN instead of 0  - not sure if TURN shouldn't be here
						insert_node (1, best_i-1, FREE);
					break;
				case 6:
					//printf("W(%d) case 6: inserting WMB(%d,%d) and Free (0,%d)\n",j,best_i+1,j,best_i);
					insert_node (best_i+1, j, P_P);
					if (best_i >= 1) // Hosna, March 26, 2012, was best_i-1 instead of best_i
						insert_node (1, best_i, FREE);
					break;
				case 7:
					//printf("W(%d) case 7: inserting WMB(%d,%d) and Free (0,%d)\n",j,best_i,j-1,best_i-1);
					insert_node (best_i, j-1, P_P);
					if (best_i-1 > 1)
						insert_node (1, best_i-1, FREE);
					break;
				case 8:
					//printf("W(%d) case 8: inserting WMB(%d,%d) and Free (0,%d)\n",j,best_i+1,j-1,best_i);
					insert_node (best_i+1, j-1, P_P);
					if (best_i >= 1) // Hosna, March 26, 2012, was best_i-1 instead of best_i
						insert_node (1, best_i, FREE);
					break;
			}
		}
		break;
		case M_WM:
		{
			cand_pos_t i = cur_interval->i;
			cand_pos_t j = cur_interval->j;
			energy_t min = INF;
			cand_pos_t best_k = j, best_row;

			min = V->get_energy_WM(i,j-1)+params_->MLbase;
			best_row = 5;
			
			  

			for (cand_pos_t k=i; k <= j-TURN-1; k++)
			{	energy_t m1 = INF,m2 = INF;
				m1 = static_cast<energy_t>((k-i)*params_->MLbase) + V->get_energy_WMv (k, j);
				if (m1 < min){
					min = m1;
					best_k = k;
					best_row = 1;
				}
				m2 = static_cast<energy_t>((k-i)*params_->MLbase) + V->get_energy_WMp(k, j);
				if (m2 < min){
					min = m2;
					best_k = k;
					best_row = 2;
				}
				energy_t m3 = V->get_energy_WM (i, k-1) + V->get_energy_WMv (k, j);
				if (m3 < min){
					min = m3;
					best_k = k;
					best_row = 3;
				}
				energy_t m4 = V->get_energy_WM (i, k-1) + V->get_energy_WMp(k, j);
				if (m4 < min){
					min = m4;
					best_k = k;
					best_row = 4;
				}
			}
			switch (best_row)
			{
				case 1: insert_node (best_k, j, M_WMv); break;
				case 2: insert_node (best_k, j, M_WMp); break;
				case 3:
					insert_node (i, best_k-1, M_WM);
					insert_node (best_k, j, M_WMv);
				break;
				case 4:
					insert_node (i, best_k-1, M_WM);
					insert_node (best_k+1, j, M_WMp);
				break;
				case 5:
					insert_node (i,j-1,M_WM); break;
			}
		}
		break;
		case M_WMv:
		{
			cand_pos_t i = cur_interval->i;
			cand_pos_t j = cur_interval->j;
			energy_t min = INF;
			cand_pos_t best_row;
			cand_pos_t si = S_[i];
			cand_pos_t sj = S_[j];
			cand_pos_t si1 = (i>1) ? S_[i-1] : -1;
			cand_pos_t sj1 = (j<n) ? S_[j+1] : -1;
			pair_type tt = pair[S_[i]][S_[j]];
			min = V->get_energy(i,j) + ((params_->model_details.dangles == 2) ? E_MLstem(tt,si1,sj1,params_) : E_MLstem(tt,-1,-1,params_));
			best_row = 1;

			if(params_->model_details.dangles == 1){
				tt = pair[S_[i+1]][S_[j]];
				energy_t tmp = V->get_energy(i+1,j) + E_MLstem(tt,si,-1,params_) + params_->MLbase;
				if(tmp<min){
					min = tmp;
					best_row = 2;
				}
				tt = pair[S_[i]][S_[j-1]];
				tmp = V->get_energy(i,j-1) + E_MLstem(tt,-1,sj,params_) + params_->MLbase;
				if(tmp<min){
					min = tmp;
					best_row = 3;
				}
				tt = pair[S_[i+1]][S_[j-1]];
				tmp = V->get_energy(i+1,j-1) + E_MLstem(tt,si,sj,params_) + 2*params_->MLbase;
				if(tmp<min){
					min = tmp;
					best_row = 4;
				}
			}

			energy_t tmp = V->get_energy_WMv(i,j-1) + params_->MLbase;
			if(tmp< min){
				min = tmp;
				best_row = 5;
			}
			switch (best_row){
				case 1: insert_node (i, j, LOOP); break;
				case 2: insert_node (i+1, j, LOOP); break;
				case 3: insert_node (i, j-1, LOOP); break;
				case 4: insert_node (i+1, j-1, LOOP); break;
				case 5: insert_node (i, j-1, M_WMv); break;
			}
		}
		break;
		case M_WMp:
		{
			int i = cur_interval->i;
			int j = cur_interval->j;
			int min = INF;
			int best_row;

			min = P->get_energy(i,j) + PSM_penalty + b_penalty;
			best_row = 1;
			energy_t tmp = V->get_energy_WMp(i,j-1) + params_->MLbase;
			if(tmp< min){
				min = tmp;
				best_row = 2;
			}
			switch (best_row){
				// case 1: insert_node (i, j, P_WMB); break;
				case 2: insert_node (i, j-1, M_WMp); break;
			}
		}
		break;
		case P_PK:
		case P_PL:
		case P_PR:
		case P_PM:
		case P_PO:
		case P_PfromL:
		case P_PfromR:
		case P_PfromM:
		case P_PfromO:
		case P_PLiloop:
		case P_PLiloop5:
		case P_PLmloop:
		case P_PLmloop0:
		case P_PLmloop1:
		case P_PRiloop:
		case P_PRiloop5:
		case P_PRmloop:
		case P_PRmloop0:
		case P_PRmloop1:
		case P_PMiloop:
		case P_PMiloop5:
		case P_PMmloop:
		case P_PMmloop0:
		case P_PMmloop1:
		case P_POiloop:
		case P_POiloop5:
		case P_POmloop:
		case P_POmloop0:
		case P_POmloop1:
		case P_WB:
		case P_WBP:
		case P_WP:
		case P_WPP:
		case P_P:
		{
			P->set_stack_interval(stack_interval);
			P->backtrack(f,cur_interval);
			stack_interval = P->get_stack_interval();
			f = P->get_minimum_fold();
		}
			break;
		default:
			printf("Should not be here!\n");
	}


}

void W_final::insert_node (int i, int j, char type)
  // insert at the beginning
{
    seq_interval *tmp;
    tmp = new seq_interval;
    tmp->i = i;
    tmp->j = j;
    tmp->type = type;
    tmp->next = stack_interval;
    stack_interval = tmp;
}

void W_final::fill_structure()
{

    std::vector<cand_pos_t> pair;
	for(cand_pos_t i=1;i<=n;++i){
		cand_pos_t j = f[i].pair;
		if(structure[i] == '(' || structure[i] == ')' || j<0) continue;

		if(j<i){
			if(!pair.empty() && i==pair.back()){
				pair.pop_back();
				structure[j] = '(';
				structure[i] = ')';
			}
		}
		else{
			if(pair.empty()){
				pair.push_back(j);
			}
			else{
				if(j<pair.back()) pair.push_back(j);
				else{
					structure[i] = '[';
					structure[j] = ']';
				}
			}
		}
	}

}