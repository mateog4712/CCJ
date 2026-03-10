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
	PK.init(n,index3D);
	PL.init(n,index3D);
	PR.init(n,index3D);
	PM.init(n,index3D);
	PO.init(n,index3D);
	PfromL.init(n,index3D);
	PfromR.init(n,index3D);
	PfromM.init(n,index3D);
	PfromMprime.init(n,index3D);
	PfromO.init(n,index3D);

	PLmloop00.init(n,index3D);
	PLmloop01.init(n,index3D);
	PLmloop10.init(n,index3D);

	PRmloop00.init(n,index3D);
	PRmloop01.init(n,index3D);
	PRmloop10.init(n,index3D);

	PMmloop00.init(n,index3D);
	PMmloop01.init(n,index3D);
	PMmloop10.init(n,index3D);

	POmloop00.init(n,index3D);
	POmloop01.init(n,index3D);
	POmloop10.init(n,index3D);
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

			compute_PLmloop00(i,j,k,l);

			compute_PLmloop01(i,j,k,l);

			compute_PLmloop10(i,j,k,l);

			compute_PRmloop00(i,j,k,l);

			compute_PRmloop01(i,j,k,l);

			compute_PRmloop10(i,j,k,l);

			compute_PMmloop00(i,j,k,l);

			compute_PMmloop01(i,j,k,l);

			compute_PMmloop10(i,j,k,l);

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

			compute_PfromMprime(i,j,k,l);

			compute_PfromO(i,j,k,l);

			compute_PK(i,j,k,l);

		}
	}

}

void pseudo_loop::compute_WBP(int i, int l){
	energy_t min_energy= INF, b1 = INF, b2=INF, b3 = INF,tmp =INF;

	cand_pos_t il = index[i]+l-i;
	cand_pos_t ilm1 = index[i]+l-1-i;
	for(cand_pos_t d=i; d< l; ++d){
		tmp = get_WB(i,d-1) + V->get_energy(d,l) + beta2P(l,d) + PPS_penalty;
		b1 = std::min(b1,tmp);
		tmp = get_WB(i,d-1) + P.get(d,l) + PSM_penalty + PPS_penalty;
		b2 = std::min(b2,tmp);
	}
	b3 = WBP[ilm1]+cp_penalty;
	min_energy = std::min({b1,b2,b3});
	if (min_energy < INF/2){
		WBP[il] = min_energy;
	}
}

void pseudo_loop::compute_WPP(cand_pos_t i, cand_pos_t l){
	energy_t min_energy = INF, b1 = INF, b2=INF, b3 =INF, tmp = INF;

	cand_pos_t il = index[i]+l-i;
	cand_pos_t ilm1 = index[i]+l-1-i;
	for(cand_pos_t d=i; d<l; ++d){
		tmp = get_WP(i,d-1) + V->get_energy(d,l) + gamma2(l,d) + PPS_penalty;
		b1 = std::min(b1,tmp);
		tmp = get_WP(i,d-1) + P.get(d,l) + PSP_penalty + PPS_penalty;
		b2 = std::min(b2,tmp);
	}
	b3 = WPP[ilm1]+PUP_penalty;
	min_energy = std::min({b1,b2,b3});
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
				b1 = PK.get(i,j,d+1,k) + PK.get(j+1,d,k+1,l);
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

	for(cand_pos_t d=i+1; d< j; ++d){
		energy_t tmp = PK.get(i,d,k,l) + get_WP(d+1,j);  //12G1
		b1=std::min(b1,tmp);
	}

	for(cand_pos_t d=k+1; d< l; ++d){
		energy_t tmp = PK.get(i,j,d,l) + get_WP(k,d-1); //1G21
		b2=std::min(b2,tmp);
	}

	b3 = PL.get(i,j,k,l) + gamma2(j,i)+PB_penalty;
	b4 = PM.get(i,j,k,l) + gamma2(j,k)+PB_penalty;
	b5 = PR.get(i,j,k,l) + gamma2(l,k)+PB_penalty;
	b6 = PO.get(i,j,k,l) + gamma2(l,i)+PB_penalty;
	min_energy = std::min({min_energy,b1,b2,b3,b4,b5,b6});
	if (min_energy < INF/2){
		PK.set(i,j,k,l,min_energy);
	}

}

void pseudo_loop::compute_PL(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	energy_t min_energy = INF,b1=INF,b2=INF,b3=INF;

	const int ptype_closing = pair[S_[i]][S_[j]];

	if (ptype_closing>0){

		b1 = get_PLiloop(i,j,k,l);

		b2 = get_PLmloop(i,j,k,l) + bp_penalty;

		// Hosna, July 11, 2014
		// To avoid addition of close base pairs we check for the following here
		if (j>=(i+TURN+1)){
			b3 = PfromL.get(i+1,j-1,k,l) + gamma2(j,i);
		}
	}
	min_energy = std::min({min_energy,b1,b2,b3});
	if (min_energy < INF/2){
		PL.set(i,j,k,l,min_energy);
	}
}

void pseudo_loop::compute_PR(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	energy_t min_energy = INF,b1=INF,b2=INF,b3=INF;

	const int ptype_closing = pair[S_[k]][S_[l]];

	if (ptype_closing>0){
		b1 = get_PRiloop(i,j,k,l);

		b2 = get_PRmloop(i,j,k,l)+ bp_penalty;

		// Hosna, July 11, 2014
		// To avoid addition of close base pairs we check for the following here
		if (l>=(k+TURN+1)){
			b3 = PfromR.get(i,j,k+1,l-1) + gamma2(l,k);
		}
	}
	min_energy = std::min({min_energy,b1,b2,b3});
	if (min_energy < INF/2){
		PR.set(i,j,k,l,min_energy);
	}
}

void pseudo_loop::compute_PM(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	energy_t min_energy = INF,b1=INF,b2=INF,b3=INF,b4=INF;

	const int ptype_closing = pair[S_[j]][S_[k]];
	if (ptype_closing>0){
		b1 = get_PMiloop(i,j,k,l);

		b2 = get_PMmloop(i,j,k,l) + bp_penalty;

		// Hosna, July 11, 2014
		// To avoid addition of close base pairs we check for the following here
		if (k>=(j+TURN-1)){
			b3 = PfromM.get(i,j-1,k+1,l) + gamma2(j,k);
		}

		if(i==j && k==l){
			b4=gamma2(i,l);
		}
	}
	min_energy = std::min({min_energy,b1,b2,b3,b4});
	if (min_energy < INF/2){
		PM.set(i,j,k,l,min_energy);
	}
}

void pseudo_loop::compute_PO(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	energy_t min_energy = INF,b1=INF,b2=INF,b3=INF;

	const int ptype_closing = pair[S_[i]][S_[l]];

	if (ptype_closing>0){
		b1 = get_POiloop(i,j,k,l);

		b2 = get_POmloop(i,j,k,l)+ bp_penalty;

		// Hosna, July 11, 2014
		// To avoid addition of close base pairs we check for the following here
		if (l>=(i+TURN+1)){
			b3 = PfromO.get(i+1,j,k,l-1) + gamma2(l,i);
		}
	}
	min_energy = std::min({min_energy,b1,b2,b3});
	if (min_energy < INF/2){
		PO.set(i,j,k,l,min_energy);
	}
}

void pseudo_loop::compute_PfromL(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	energy_t min_energy = INF,b1=INF,b2=INF,b3=INF,b4=INF,b5=INF;

	for(cand_pos_t d=i+1; d< j; ++d){
		energy_t tmp = PfromL.get(d,j,k,l)+get_WP(i,d-1);
		b1 = std::min(b1,tmp);
		tmp = PfromL.get(i,d,k,l)+get_WP(d+1,j);
		b2 = std::min(b2,tmp);
	}

	b3 = PR.get(i,j,k,l) + gamma2(l,k) + PB_penalty; //;

	b4 = PM.get(i,j,k,l) + gamma2(j,k)+ PB_penalty;//;

	b5 = PO.get(i,j,k,l) + gamma2(l,i)+ PB_penalty;//;

	min_energy = std::min({min_energy,b1,b2,b3,b4,b5});
	if (min_energy < INF/2){
		PfromL.set(i,j,k,l,min_energy);
	}
}

void pseudo_loop::compute_PfromR(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	energy_t min_energy = INF,b1=INF,b2=INF,b3=INF,b4=INF;

	for(cand_pos_t d=k+1; d< l; ++d){
		energy_t tmp = PfromR.get(i,j,d,l)+get_WP(k,d-1);
		b1 = std::min(b1,tmp);
		tmp = PfromR.get(i,j,k,d)+get_WP(d+1,l);
		b2 = std::min(b2,tmp);
	}

	b3 = PM.get(i,j,k,l) + gamma2(j,k)+ PB_penalty;//;

	b4 = PO.get(i,j,k,l) + gamma2(l,i)+ PB_penalty;//;

	min_energy = std::min({min_energy,b1,b2,b3,b4});
	if (min_energy < INF/2){
		PfromR.set(i,j,k,l,min_energy);
	}
}

void pseudo_loop::compute_PfromM(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	energy_t min_energy = INF;

	for(cand_pos_t d=i+1; d<j; ++d){
		energy_t tmp = PfromMprime.get(i,d,k,l)+get_WP(d+1,j);
		min_energy = std::min(min_energy,tmp);
	}

	if (min_energy < INF/2){
		PfromM.set(i,j,k,l,min_energy);
	}
}

void pseudo_loop::compute_PfromMprime(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	energy_t min_energy = INF;

	for(cand_pos_t d=k+1; d<l; ++d){
		energy_t tmp=get_PfromMdoubleprime(i,j,d,l)+get_WP(k,d-1);
		min_energy = std::min(min_energy,tmp);
	}
	
	if (min_energy < INF/2){
		PfromM.set(i,j,k,l,min_energy);
	}
}

void pseudo_loop::compute_PfromO(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	energy_t min_energy = INF,b1=INF,b2=INF,b3=INF,b4=INF;

	for(cand_pos_t d=i+1; d< j; ++d){
		energy_t tmp=PfromO.get(d,j,k,l)+get_WP(i,d-1);
		b1 = std::min(b1,tmp);
	}
	for(cand_pos_t d=k+1; d<l; ++d){
		energy_t tmp=PfromO.get(i,j,k,d)+get_WP(d+1,l);
		b2 = std::min(b2,tmp);
	}

	b3 = PL.get(i,j,k,l) + gamma2(j,i) + PB_penalty;

	b4 = PR.get(i,j,k,l) + gamma2(l,k)+PB_penalty;

	min_energy = std::min({min_energy,b1,b2,b3,b4});
	if (min_energy < INF/2){
		PfromO.set(i,j,k,l,min_energy);
	}

}

void pseudo_loop::compute_PLmloop00(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	energy_t min_energy = INF,tmp=INF;

	min_energy = PL.get(i,j,k,l)+beta2P(j,i);
    for(cand_pos_t d = i; d<=j; ++d){
        if (d>i){
            tmp=get_WB(i,d-1)+PLmloop00.get(d,j,k,l);
            min_energy = std::min(min_energy,tmp);
        }
        if(d<j){
            tmp=PLmloop00.get(i,d,k,l)+get_WB(d+1,j);
            min_energy = std::min(min_energy,tmp);
        }

    }
	if (min_energy < INF/2){
		PLmloop00.set(i,j,k,l,min_energy);
	}
}

void pseudo_loop::compute_PLmloop01(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	energy_t min_energy = INF,tmp=INF;

	for(cand_pos_t d = i; d < j; ++d){
        tmp = PLmloop00.get(i,d,k,l) + WBP.get(d+1,j);
        min_energy = std::min(min_energy,tmp);
    }
	if (min_energy < INF/2){
		PLmloop01.set(i,j,k,l,min_energy);
	}

}

void pseudo_loop::compute_PLmloop10(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	energy_t min_energy = INF,tmp=INF;

	for(cand_pos_t d = i+1; d <= j; ++d){
        tmp = WBP.get(i,d-1) + PLmloop00.get(d,j,k,l);
        min_energy = std::min(min_energy,tmp);
        if(d<j){
            tmp = PLmloop10.get(i,d,k,l) + get_WB(d+1,j);
            min_energy = std::min(min_energy,tmp);
        }
    }
	if (min_energy < INF/2){
		PLmloop10.set(i,j,k,l,min_energy);
	}

}

void pseudo_loop::compute_PRmloop00(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	energy_t min_energy = INF,tmp=INF;

	min_energy = PR.get(i,j,k,l)+beta2P(l,k);
	for(cand_pos_t d=k; d<=l; ++d){
        if(d>k){
            tmp=get_WB(k,d-1)+PRmloop00.get(i,j,d,l);
            min_energy = std::min(min_energy,tmp);
        }
        if (d<l){
            tmp = PRmloop00.get(i,j,k,d)+get_WB(d+1,l);
            min_energy = std::min(min_energy,tmp);
        }
    }
	if (min_energy < INF/2){
		PRmloop00.set(i,j,k,l,min_energy);
	}

}


void pseudo_loop::compute_PRmloop01(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	energy_t min_energy = INF,tmp=INF;

	min_energy = PRmloop01.get(i,j,k,l-1)+cp_penalty;
    for(cand_pos_t d = k; d < l; d++){
        tmp = PRmloop00.get(i,j,k,d) + WBP.get(d+1,l);
        min_energy = std::min(min_energy,tmp);
    }
	if (min_energy < INF/2){
		PRmloop01.set(i,j,k,l,min_energy);
	}

}

void pseudo_loop::compute_PRmloop10(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	energy_t min_energy = INF,tmp=INF;

	min_energy = PRmloop10.get(i,j,k+1,l)+ cp_penalty;
    for(cand_pos_t d = k+1; d <= l; ++d){
        tmp = WBP.get(k,d-1) + PRmloop00.get(i,j,d,l);
        min_energy = std::min(min_energy,tmp);
    }
	if (min_energy < INF/2){
		PRmloop10.set(i,j,k,l,min_energy);
	}

}

void pseudo_loop::compute_PMmloop00(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	energy_t min_energy = INF,tmp=INF;

	min_energy = PM.get(i,j,k,l)+beta2P(j,k);
    for(cand_pos_t d=i; d<j; ++d){
        tmp=PMmloop00.get(i,d,k,l)+get_WB(d+1,j);
        min_energy = std::min(min_energy,tmp);
    }
    for(cand_pos_t d=k+1; d<=l; ++d){
        tmp=PMmloop00.get(i,j,d,l)+get_WB(k,d-1);
        min_energy = std::min(min_energy,tmp);
    }

	if (min_energy < INF/2){
		PMmloop00.set(i,j,k,l,min_energy);
	}
}


void pseudo_loop::compute_PMmloop01(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	energy_t min_energy = INF,tmp=INF;

	min_energy = PMmloop01.get(i,j,k+1,l)+cp_penalty;
	for(cand_pos_t d = k; d < l; ++d){
        tmp = POmloop00.get(i,j,k,d) + WBP.get(d+1,l);
        min_energy = std::min(min_energy,tmp);
    }

	if (min_energy < INF/2){
		PMmloop01.set(i,j,k,l,min_energy);
	}
}

void pseudo_loop::compute_PMmloop10(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	energy_t min_energy = INF,tmp=INF;

	min_energy = PMmloop10.get(i,j-1,k,l)+cp_penalty;
	for(cand_pos_t d=i+1; d<=j;++d){
        tmp=WBP.get(i,d-1)+POmloop00.get(d,j,k,l);
        min_energy = std::min(min_energy,tmp);
    }
    for(cand_pos_t d = k+1; d < l; ++d){
        tmp = POmloop10.get(i,j,k,d) + get_WB(d+1,l);
        min_energy = std::min(min_energy,tmp);
    }

	if (min_energy < INF/2){
		PMmloop10.set(i,j,k,l,min_energy);
	}
}

void pseudo_loop::compute_POmloop00(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	energy_t min_energy = INF,tmp=INF;

	min_energy = PO.get(i,j,k,l)+beta2P(l,i);
    for(cand_pos_t d=i+1; d<=j; ++d){
        tmp = get_WB(i,d-1)+POmloop00.get(d,j,k,l);
        min_energy = std::min(min_energy,tmp);
    }
    for(cand_pos_t d=k; d<l; ++d){
        tmp = POmloop00.get(i,j,k,d)+get_WB(d+1,l);
        min_energy = std::min(min_energy,tmp);
    }

	if (min_energy < INF/2){
		POmloop00.set(i,j,k,l,min_energy);
	}

}


void pseudo_loop::compute_POmloop01(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	energy_t min_energy = INF;

	for(cand_pos_t d = k; d < l; ++d){
        energy_t tmp = POmloop00.get(i,j,k,d) + WBP.get(d+1,l);
        min_energy = std::min(min_energy,tmp);
    }

	if (min_energy < INF/2){
		POmloop01.set(i,j,k,l,min_energy);
	}

}

void pseudo_loop::compute_POmloop10(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	energy_t min_energy = INF,tmp=INF;

	for(cand_pos_t d=i+1; d<=j;++d){
        tmp=WBP.get(i,d-1)+POmloop00.get(d,j,k,l);
        min_energy = std::min(min_energy,tmp);
    }
    for(cand_pos_t d = k+1; d < l; ++d){
        tmp = POmloop10.get(i,j,k,d) + get_WB(d+1,l);
        min_energy = std::min(min_energy,tmp);
    }

	if (min_energy < INF/2){
		POmloop10.set(i,j,k,l,min_energy);
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

energy_t pseudo_loop::get_PfromMdoubleprime(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	if (!(i <= j && j < k-1 && k <= l)){
		return INF;
	}
	assert(!(i<=0 || l> n));

	if (i==j && k==l){
		const int ptype_closing = pair[S_[i]][S_[l]];
		if(ptype_closing == 0) return INF;
		else return 0;
	}

	energy_t b1 = PL.get(i,j,k,l) + gamma2(j,i)+ PB_penalty;

	energy_t b2 = PR.get(i,j,k,l) + gamma2(l,k)+PB_penalty;

	return std::min(b1,b2);
}

// I will decide at a later time if this should be grabbing from a matrix or calculated every time it is called
energy_t pseudo_loop::get_PLiloop(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	if (!(i <= j && j < k-1 && k <= l)){
		return INF;
	}
	assert(!(i<=0 || l> n));

	const int ptype_closing = pair[S_[i]][S_[j]];
	if(ptype_closing == 0) return INF;

	energy_t tmp=INF;
	energy_t min_energy = PL.get(i+1,j-1,k,l) + get_e_stP(i,j);

	for(cand_pos_t d= i+1; d<std::min(j,i+MAXLOOP); ++d){
		for(cand_pos_t dp = j-1; dp > std::max(d+TURN,j-MAXLOOP); --dp){
			if (!can_pair(d,dp)) continue;
			tmp = get_e_intP(i,d,dp,j) + PL.get(d,dp,k,l);
			min_energy = std::min(min_energy,tmp);
		}
	}
	return min_energy;
}

energy_t pseudo_loop::get_PLmloop(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	if (!(i <= j && j < k-1 && k <= l)){
		return INF;
	}
	assert(!(i<=0 || l> n));

	energy_t b1 = PLmloop10.get(i+1,j-1,k,l)+ ap_penalty + beta2P(j,i);
    energy_t b2 = PLmloop01.get(i+1,j-1,k,l)+ ap_penalty + beta2P(j,i);

	return std::min(b1,b2);
}

energy_t pseudo_loop::get_PRiloop(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	if (!(i <= j && j < k-1 && k <= l)){
		return INF;
	}
	assert(!(i<=0 || l> n));

	const int ptype_closing = pair[S_[k]][S_[l]];
	if(ptype_closing == 0) return INF;

	energy_t tmp = INF;
	energy_t min_energy = PR.get(i,j,k+1,l-1) + get_e_stP(k,l);

	for(cand_pos_t d= k+1; d<std::min(l,k+MAXLOOP); ++d){
		for(cand_pos_t dp=l-1; dp > std::max(d+TURN,l-MAXLOOP); --dp){
			if (!can_pair(d,dp)) continue;
			tmp = get_e_intP(k,d,dp,l) + PR.get(i,j,d,dp);
			min_energy = std::min(min_energy,tmp);
		}
	}
	return min_energy;
}

energy_t pseudo_loop::get_PRmloop(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	if (!(i <= j && j < k-1 && k <= l)){
		return INF;
	}
	assert(!(i<=0 || l> n));

	energy_t b1 = PRmloop10.get(i,j,k+1,l-1) + ap_penalty + beta2P(l,k);
    energy_t b2 = PRmloop01.get(i,j,k+1,l-1) + ap_penalty + beta2P(l,k);

	return std::min(b1,b2);
}

energy_t pseudo_loop::get_PMiloop(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	if (!(i <= j && j < k-1 && k <= l)){
		return INF;
	}
	assert(!(i<=0 || l> n));

	const int ptype_closing = pair[S_[j]][S_[k]];
	if(ptype_closing == 0) return INF;

	energy_t tmp=INF;

	energy_t min_energy = PM.get(i,j-1,k+1,l) + get_e_stP(j-1,k+1);

	for(cand_pos_t d= j-1; d>std::max(i,j-MAXLOOP); --d){
		for (cand_pos_t dp=k+1; dp <std::min(l,k+MAXLOOP); ++dp) {
			if (!can_pair(d,dp)) continue;
			tmp = get_e_intP(d,j,k,dp) + PM.get(i,d,dp,l);
			min_energy = std::min(min_energy,tmp);
		}
	}
	return min_energy;
}

energy_t pseudo_loop::get_PMmloop(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	if (!(i <= j && j < k-1 && k <= l)){
		return INF;
	}
	assert(!(i<=0 || l> n));

	energy_t b1 = PMmloop10.get(i,j-1,k+1,l)+ap_penalty+beta2P(j,k);
    energy_t b2 = PMmloop01.get(i,j-1,k+1,l)+ap_penalty+beta2P(j,k);

	return std::min(b1,b2);
}

energy_t pseudo_loop::get_POiloop(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	if (!(i <= j && j < k-1 && k <= l)){
		return INF;
	}
	assert(!(i<=0 || l> n));

	const int ptype_closing = pair[S_[i]][S_[l]];
	if(ptype_closing == 0) return INF;

	energy_t tmp = INF;
	energy_t min_energy = PO.get(i+1,j,k,l-1) + get_e_stP(i,l);

	for(cand_pos_t d= i+1; d<std::min(j,i+MAXLOOP); ++d){
		for (cand_pos_t dp=l-1; dp >std::max(l-MAXLOOP,k); --dp) {
			if (!can_pair(d,dp)) continue;
			tmp = get_e_intP(i,d,dp,l) + PO.get(d,j,dp,k);
			min_energy = std::min(min_energy,tmp);
		}
	}
	return min_energy;
}

energy_t pseudo_loop::get_POmloop(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l){
	if (!(i <= j && j < k-1 && k <= l)){
		return INF;
	}
	assert(!(i<=0 || l> n));

	energy_t b1 = POmloop10.get(i+1,j,k,l-1) + ap_penalty + beta2P(l,i);
    energy_t b2 = POmloop01.get(i+1,j,k,l-1) + ap_penalty + beta2P(l,i);

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
						b1 = PK.get(i,j,d+1,k) + PK.get(j+1,d,k+1,l);
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
			for(cand_pos_t  d=i+1; d< j; ++d){
				temp = PK.get(i,d,k,l) + get_WP(d+1,j);
				if (temp < min_energy){
					min_energy=temp;
					best_row = 1;
					best_d=d;
				}
			}

			// branch 2
			for(cand_pos_t d=k+1; d< l; ++d){
				temp = PK.get(i,j,d,l) + get_WP(k,d-1);
				if (temp < min_energy){
					min_energy=temp;
					best_row = 2;
					best_d = d;
				}
			}

			// branch 3
			temp = PL.get(i,j,k,l) + gamma2(j,i)+PB_penalty;
			if(temp < min_energy){
				min_energy = temp;
				best_row = 3;
				best_d = -1;
			}

			//branch 4
			temp = PM.get(i,j,k,l) + gamma2(j,k)+PB_penalty;
			if(temp < min_energy){
				min_energy = temp;
				best_row = 4;
				best_d = -1;
			}

			// branch 5
			temp = PR.get(i,j,k,l) + gamma2(l,k)+PB_penalty;
			if(temp < min_energy){
				min_energy = temp;
				best_row = 5;
				best_d = -1;
			}

			// branch 6
			temp = PO.get(i,j,k,l) + gamma2(l,i)+PB_penalty;
			if(temp < min_energy){
				min_energy = temp;
				best_row = 6;
				best_d = -1;
			}

			switch (best_row){
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
				temp = get_PLmloop(i,j,k,l) + bp_penalty;
				if(temp < min_energy){
					min_energy = temp;
					best_row=2;
				}

				//branch 3
				if (j>=(i+TURN+1)){
					temp = PfromL.get(i+1,j-1,k,l) + gamma2(j,i);
					if(temp < min_energy){
						min_energy = temp;
						best_row= 3;
					}
				}
			}

			switch (best_row){
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
				temp = get_PRmloop(i,j,k,l)+ bp_penalty;
				if(temp < min_energy){
					min_energy = temp;
					best_row=2;
				}

				//branch 3
				if (l>=(k+TURN+1)){
					temp = PfromR.get(i,j,k+1,l-1) + gamma2(l,k);
					if(temp < min_energy){
						min_energy = temp;
						best_row= 3;
					}
				}
			}

			switch (best_row){
				case 1:
					insert_node(i,l,j,k,P_PRiloop);
					break;
				case 2:
					insert_node(i,l,j,k,P_PRmloop);
					break;
				case 3:
					insert_node(i,l-1,j,k+1,P_PfromR);
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
				temp = get_PMmloop(i,j,k,l) + bp_penalty;
				if(temp < min_energy){
					min_energy = temp;
					best_row=2;
				}
				//branch 3
				if (k>=(j+TURN-1)){
					temp = PfromM.get(i,j-1,k+1,l) + gamma2(j,k);
					if(temp < min_energy){
						min_energy = temp;
						best_row= 3;
					}
				}
			}

			switch (best_row){
				case 1:
					insert_node(i,l,j,k,P_PMiloop);
					break;
				case 2:
					insert_node(i,l,j,k,P_PMmloop);
					break;
				case 3:
					insert_node(i,l,j-1,k+1,P_PfromM);
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
				temp = get_POmloop(i,j,k,l)+bp_penalty;
				if(temp < min_energy){
					min_energy = temp;
					best_row=2;
				}
				//branch 3
				if (l>=(i+TURN+1)){
					temp = PfromO.get(i+1,j,k,l-1) + gamma2(l,i);
					if(temp < min_energy){
						min_energy = temp;
						best_row= 3;
					}
				}
			}
			switch (best_row){
				case 1:
					insert_node(i,l,j,k,P_POiloop);
					break;
				case 2:
					insert_node(i,l,j,k,P_POmloop);
					break;
				case 3:
					insert_node(i+1,l-1,j,k,P_PfromO);
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
				temp=PfromL.get(d,j,k,l)+get_WP(i,d-1);
				if(temp < min_energy){
					min_energy=temp;
					best_row=1;
					best_d = d;

				}
				//branch 2
				temp=PfromL.get(i,d,k,l)+get_WP(d+1,j);
				if(temp < min_energy){
					min_energy=temp;
					best_row=2;
					best_d = d;
				}
			}
			// branch 3
			temp = PR.get(i,j,k,l) + gamma2(l,k)+PB_penalty;
			if(temp < min_energy){
				min_energy = temp;
				best_row=3;
				best_d = -1;
			}

			//branch 4
			temp = PM.get(i,j,k,l) + gamma2(j,k)+PB_penalty;
			if(temp < min_energy){
				min_energy = temp;
				best_row=4;
				best_d = -1;
			}

			// branch 5
			temp = PO.get(i,j,k,l) + gamma2(l,i)+PB_penalty;
			if(temp < min_energy){
				min_energy = temp;
				best_row=5;
				best_d = -1;
			}
			switch (best_row){
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
				temp=PfromR.get(i,j,d,l)+get_WP(k,d-1);
				if(temp < min_energy){
					min_energy = temp;
					best_row=1;
					best_d = d;

				}
				//branch 2
				temp=PfromR.get(i,j,k,d)+get_WP(d+1,l);
				if(temp < min_energy){
					min_energy=temp;
					best_row=2;
					best_d = d;
				}
			}

			//branch 3
			temp = PM.get(i,j,k,l) + gamma2(j,k)+PB_penalty;
			if(temp < min_energy){
				min_energy = temp;
				best_row=3;
				best_d = -1;
			}

			// branch 4
			temp = PO.get(i,j,k,l) + gamma2(l,i)+PB_penalty;
			if(temp < min_energy){
				min_energy = temp;
				best_row=4;
				best_d = -1;
			}
			switch (best_row){
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
			cand_pos_t best_d=-1;

			for(cand_pos_t d=i+1; d< j; d++){
				//branch 1
				temp=PfromMprime.get(i,d,k,l)+get_WP(d+1,j);
				if(temp < min_energy){
					min_energy=temp;
					best_d = d;
				}
			}

			if (best_d > -1){
				insert_node(i,l,best_d,k,P_PfromMprime);
				insert_node(best_d+1,j,P_WP);
			}


		}
			break;

			case P_PfromMprime:
		{
			cand_pos_t i = cur_interval->i;
			cand_pos_t l = cur_interval->j;
			cand_pos_t j = cur_interval->k;
			cand_pos_t k = cur_interval->l;
			if(debug) printf("At P_PfromM with i=%d,l=%d,j=%d,k=%d\n",i,l,j,k);

			if (!(i <= j && j < k-1 && k <= l)){
				std::cerr << "This should not have happened!, P_PfromMprime" << std::endl;
				exit(EXIT_FAILURE);
			}

			if (i<=0 ||j<=0 || k<=0 || l<=0 || i>n || j>n || k> n || l> n){
				std::cerr << "This should not have happened!, P_PfromMprime" << std::endl;
				exit(EXIT_FAILURE);
			}

			if (i==j && k==l){
				return;
			}

			energy_t min_energy = INF,temp=INF;
			cand_pos_t best_d=-1;

			for(cand_pos_t d=k+1; d< l; d++){
				//branch 2
				temp=get_PfromMdoubleprime(i,j,d,l)+get_WP(k,d-1);
				if(temp < min_energy){
					min_energy=temp;
					best_d = d;
				}
			}

			if (best_d > -1){
				insert_node(i,l,j,best_d,P_PfromMdoubleprime);
				insert_node(k,best_d-1,P_WP);
			}

		}
			break;

			case P_PfromMdoubleprime:
		{
			cand_pos_t i = cur_interval->i;
			cand_pos_t l = cur_interval->j;
			cand_pos_t j = cur_interval->k;
			cand_pos_t k = cur_interval->l;
			if(debug) printf("At P_PfromMdoubleprime with i=%d,l=%d,j=%d,k=%d\n",i,l,j,k);

			if (!(i <= j && j < k-1 && k <= l)){
				std::cerr << "This should not have happened!, P_PfromMdoubleprime" << std::endl;
				exit(EXIT_FAILURE);
			}

			if (i<=0 ||j<=0 || k<=0 || l<=0 || i>n || j>n || k> n || l> n){
				std::cerr << "This should not have happened!, P_PfromMdoubleprime" << std::endl;
				exit(EXIT_FAILURE);
			}

			if (i==j && k==l){
				return;
			}

			energy_t min_energy = INF,temp=INF;
			cand_pos_t best_row = -1;

			// branch 1
			temp = PL.get(i,j,k,l) + gamma2(j,i)+PB_penalty;
			if(temp < min_energy){
				min_energy = temp;
				best_row=1;
			}

			//branch 2
			temp = PR.get(i,j,k,l) + gamma2(l,k)+PB_penalty;
			if(temp < min_energy){
				min_energy = temp;
				best_row=2;
			}

			switch (best_row){
				case 1:
					insert_node(i,l,j,k,P_PL);
					break;
				case 2:
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
				temp=PfromO.get(d,j,k,l)+get_WP(i,d-1);
				if(temp < min_energy){
					min_energy=temp;
					best_row=1;
					best_d = d;

				}
			}
			for(cand_pos_t d=k+1; d< l; d++){
				//branch 2
				temp=PfromO.get(i,j,k,d)+get_WP(d+1,l);
				if(temp < min_energy){
					min_energy=temp;
					best_row=2;
					best_d = d;
				}
			}
			// branch 3
			temp = PL.get(i,j,k,l) + gamma2(j,i)+PB_penalty;
			if(temp < min_energy){
				min_energy = temp;
				best_row=3;
				best_d = -1;
			}

			//branch 4
			temp = PR.get(i,j,k,l) + gamma2(l,k)+PB_penalty;
			if(temp < min_energy){
				min_energy = temp;
				best_row=4;
				best_d = -1;
			}

			switch (best_row){
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
			temp = WBP.get(i,l);
			if (temp < min_energy){
				min_energy = temp;
				best_row=1;
			}

			temp = cp_penalty*(l-i+1);
			if (temp < min_energy){
				min_energy = temp;
				best_row=2;
			}

			switch (best_row){
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
			cand_pos_t best_row = -1,best_d=-1;

			for(cand_pos_t d=i; d< l; d++){
				//branch 1
				temp = get_WB(i,d-1) + V->get_energy(d,l) +beta2P(l,d) + PPS_penalty;
				if (temp < min_energy){
					min_energy = temp;
					best_row = 1;
					best_d = d;
				}

				//branch 2
				temp = get_WB(i,d-1) + P.get(d,l) + PSM_penalty + PPS_penalty;
				if (temp < min_energy){
					min_energy = temp;
					best_row = 2;
					best_d = d;
				}
			}
			temp = WBP.get(i,l-1) + cp_penalty;
			if(temp< min_energy){
				min_energy = temp;
				best_row = 3;
			}

			switch (best_row){
				case 1:
					insert_node(i,best_d-1,P_WB);
					insert_node(best_d,l,LOOP);
					break;
				case 2:
					insert_node(i,best_d-1,P_WB);
					insert_node(best_d,l,P_P);
					break;
				case 3:
					insert_node(i,l-1,P_WBP);
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
			temp = WPP.get(i,l);
			if (temp < min_energy){
				min_energy = temp;
				best_row=1;
			}

			temp = PUP_penalty*(l-i+1);
			if (temp < min_energy){
				min_energy = temp;
				best_row=2;
			}

			switch (best_row){
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
		{
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
			cand_pos_t best_row = -1,best_d=-1;

			for(cand_pos_t d=i; d< l; d++){
				//branch 1
				temp = get_WP(i,d-1) + V->get_energy(d,l) + gamma2(l,d) + PPS_penalty;
				if (temp < min_energy){
					min_energy = temp;
					best_row = 1;
					best_d = d;
				}

				//branch 2
				temp = get_WP(i,d-1) + P.get(d,l) + PSP_penalty + PPS_penalty;
				if (temp < min_energy){
					min_energy = temp;
					best_row = 2;
					best_d = d;
				}
			}
			temp = WPP.get(i,l-1) + PUP_penalty;
			if(temp<min_energy){
				min_energy = temp;
				best_row = 3;
			}

			switch (best_row){
				case 1:
					insert_node(i,best_d-1,P_WP);
					insert_node(best_d,l,LOOP);
					break;
				case 2:
					insert_node(i,best_d-1,P_WP);
					insert_node(best_d,l,P_P);
					break;
				case 3:
					insert_node(i,l-1,P_WPP);
					break;
			}
		}
			break;
		case P_PLiloop:
		{
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

			f[i].pair = j;
			f[j].pair = i;
			f[i].type = P_PLiloop;
			f[j].type = P_PLiloop;

			energy_t min_energy = INF,temp=INF;
			cand_pos_t best_row = -1,best_d=-1,best_dp=-1;
			if (pair[S_[i]][S_[j]] > 0){
				//branch 1
				temp = PL.get(i+1,j-1,k,l) + get_e_stP(i,j);
				if (temp < min_energy){
					min_energy = temp;
					best_row=1;
				}
				//branch 2
				for(cand_pos_t d= i+1; d<std::min(j,i+MAXLOOP); ++d){
					for(cand_pos_t dp = j-1; dp > std::max(d+TURN,j-MAXLOOP); --dp){
						temp = get_e_intP(i,d,dp,j) + PL.get(d,dp,k,l);
						if(temp < min_energy){
							min_energy = temp;
							best_d = d;
							best_dp = dp;
							best_row = 2;
						}
					}
            }	
				
			}

			switch (best_row){
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

			f[i].pair = j;
			f[j].pair = i;
			f[i].type = P_PLmloop;
			f[j].type = P_PLmloop;

			energy_t branch1 = PLmloop10.get(i+1,j-1,k,l)+ ap_penalty + beta2P(j,i);
            energy_t branch2 = PLmloop01.get(i+1,j-1,k,l)+ ap_penalty + beta2P(j,i);

            cand_pos_t best_row = -1;
            if (branch1 < branch2) { best_row = 1;} else {best_row = 2;}

			switch (best_row){
				case 1:
					insert_node(i+1,l,j-1,k,P_PLmloop10);
					break;
				case 2:
					insert_node(i+1,l,j-1,k,P_PLmloop01);
					break;
			}

		}
		break;
		case P_PLmloop00:
		{
			cand_pos_t i = cur_interval->i;
			cand_pos_t l = cur_interval->j;
			cand_pos_t j = cur_interval->k;
			cand_pos_t k = cur_interval->l;
			if(debug) printf("At P_PLmloop00 with i=%d,l=%d,j=%d,k=%d\n",i,l,j,k);

			if (!(i <= j && j < k-1 && k <= l)){
				std::cerr << "border cases: This should not have happened!, P_PLmloop00" << std::endl;
				exit(EXIT_FAILURE);
			}

			if (i<=0 || j<=0 || k<=0 || l<=0 || i>n || j>n || k> n || l> n){
				std::cerr << "impossible cases: This should not have happened!, P_PLmloop00" << std::endl;
				exit(EXIT_FAILURE);
			}


			energy_t min_energy = PL.get(i,j,k,l)+beta2P(j,i),temp=INF;
			cand_pos_t best_row = 1, best_d=-1;
            for(cand_pos_t d = i; d<=j; ++d){
                if (d>i){
                    temp=get_WB(i,d-1)+PLmloop00.get(d,j,k,l);
                    if (temp < min_energy){
                        min_energy = temp;
                        best_row = 2;
                        best_d = d;
                    }

                }
                if(d<j){
                    temp=PLmloop00.get(i,d,k,l)+get_WB(d+1,j);
                    if (temp < min_energy){
                        min_energy = temp;
                        best_row = 3;
                        best_d = d;
                    }
                }
            }
			switch (best_row){
                case 1:
                    insert_node(i,l,j,k,P_PL);
                    break;
				case 2:
					insert_node(best_d,l,j,k,P_PLmloop00);
					insert_node(i,best_d-1,P_WB);
					break;
				case 3:
					insert_node(i,l,best_d,k,P_PLmloop00);
					insert_node(best_d+1,j,P_WB);
					break;
			}

		}
			break;
		case P_PLmloop01:
		{
			cand_pos_t i = cur_interval->i;
			cand_pos_t l = cur_interval->j;
			cand_pos_t j = cur_interval->k;
			cand_pos_t k = cur_interval->l;
			if(debug) printf("At P_PLmloop01 with i=%d,l=%d,j=%d,k=%d\n",i,l,j,k);

			if (!(i <= j && j < k-1 && k <= l)){
				std::cerr << "border cases: This should not have happened!, P_PLmloop01" << std::endl;
				exit(EXIT_FAILURE);
			}

			if (i<=0 || j<=0 || k<=0 || l<=0 || i>n || j>n || k> n || l> n){
				std::cerr << "impossible cases: This should not have happened!, P_PLmloop01" << std::endl;
				exit(EXIT_FAILURE);
			}

			energy_t min_energy = INF,temp=INF;
			cand_pos_t best_d=-1;
			for(cand_pos_t d = i; d < j; ++d){
                temp = PLmloop00.get(i,d,k,l) + WBP.get(d+1,j);
                if (temp < min_energy){
                    min_energy = temp;
                    best_d = d;
                }
            }
			insert_node(i,l,best_d,k,P_PLmloop00);
			insert_node(best_d+1,j,P_WBP);

		}
			break;
		case P_PLmloop10:
		{
			cand_pos_t i = cur_interval->i;
			cand_pos_t l = cur_interval->j;
			cand_pos_t j = cur_interval->k;
			cand_pos_t k = cur_interval->l;
			if(debug) printf("At P_PLmloop10 with i=%d,l=%d,j=%d,k=%d\n",i,l,j,k);

			if (!(i <= j && j < k-1 && k <= l)){
				std::cerr<< "border cases: This should not have happened!, P_PLmloop10" << std::endl;
				exit(EXIT_FAILURE);
			}

			if (i<=0 || j<=0 || k<=0 || l<=0 || i>n || j>n || k> n || l> n){
				std::cerr<< "impossible cases: This should not have happened!, P_PLmloop10" << std::endl;
				exit(EXIT_FAILURE);
			}

			energy_t min_energy = INF,temp=INF;
			cand_pos_t best_d=-1,best_row=-1;
			for(cand_pos_t d = i+1; d <= j; ++d){
                temp = WBP.get(i,d-1) + PLmloop00.get(d,j,k,l);
                if (temp < min_energy){
                    min_energy = temp;
                    best_row = 1;
                    best_d = d;
                }
                if(d<j){
                    temp = PLmloop10.get(i,d,k,l) + get_WB(d+1,j);
                    if (temp < min_energy){
                        min_energy = temp;
                        best_row = 2;
                        best_d = d;
                    }
                }
            }
			switch (best_row){
				case 1:
					insert_node(i, best_d-1, P_WBP);
                    insert_node(best_d,l,j,k, P_PLmloop00);
					break;
				case 2:
					insert_node(i,l,best_d,k, P_PLmloop10);
					insert_node(best_d+1,j,P_WB);
					break;
			}
		}
			break;

		case P_PRiloop:
		{
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

			f[k].pair = l;
			f[l].pair = k;
			f[k].type = P_PRiloop;
			f[l].type = P_PRiloop;

			energy_t min_energy = INF,temp=INF;
			cand_pos_t best_row=-1,best_d=-1,best_dp=-1;
			if (pair[S_[k]][S_[l]]>0){
				//branch 1
				temp = PR.get(i,j,k+1,l-1) + get_e_stP(k,l);
				if(temp < min_energy){
					min_energy = temp;
					best_row = 1;
				}
				//branch 2
				for(cand_pos_t d= k+1; d<std::min(l,k+MAXLOOP); ++d){
					for(cand_pos_t dp=l-1; dp > std::max(d+TURN,l-MAXLOOP); --dp){
						temp = get_e_intP(k,d,dp,l) + PR.get(i,j,d,dp);
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

			f[k].pair = l;
			f[l].pair = k;
			f[k].type = P_PRmloop;
			f[l].type = P_PRmloop;

			cand_pos_t best_row=-1;

			energy_t branch1 = PRmloop10.get(i,j,k+1,l-1) + ap_penalty + beta2P(l,k);
            energy_t branch2 = PRmloop01.get(i,j,k+1,l-1) + ap_penalty + beta2P(l,k);

            if (branch1 < branch2) {best_row = 1; } else {best_row = 2;}

			switch (best_row){
				case 1:
					insert_node(i,l-1,j,k+1,P_PRmloop10);
					break;
				case 2:
					insert_node(i,l-1,j,k+1,P_PRmloop01);
					break;
			}
		}
		break;

		case P_PRmloop00: // Mateo - This was omitted in the sparse code; not sure why.
		{
			cand_pos_t i = cur_interval->i;
			cand_pos_t l = cur_interval->j;
			cand_pos_t j = cur_interval->k;
			cand_pos_t k = cur_interval->l;
			if(debug) printf("At P_PRmloop00 with i=%d,l=%d,j=%d,k=%d\n",i,l,j,k);

			if (!(i <= j && j < k-1 && k <= l)){
				std::cerr << "border cases: This should not have happened!, P_PRmloop00" << std::endl;
				exit(EXIT_FAILURE);
			}

			if (i<=0 || j<=0 || k<=0 || l<=0 || i>n || j>n || k> n || l> n){
				std::cerr << "impossible cases: This should not have happened!, P_PRmloop00" << std::endl;
				exit(EXIT_FAILURE);
			}

			energy_t min_energy = INF,temp=INF;
			cand_pos_t best_d=-1;
            min_energy = PR.get(i,j,k,l)+beta2P(l,k);
			cand_pos_t best_row=1;
			for(cand_pos_t d=k; d<=l; ++d){
				if(d>k){
					temp=get_WB(k,d-1)+PRmloop00.get(i,j,d,l);
					if (temp < min_energy){
						min_energy = temp;
						best_row = 2;
						best_d = d;
                	}
				}
				if (d<l){
					temp = PRmloop00.get(i,j,k,d)+get_WB(d+1,l);
					if (temp < min_energy){
						min_energy = temp;
						best_row = 3;
						best_d = d;
                	}
				}
			}
			switch (best_row){
                case 1:
                    insert_node(i,j,k,l,P_PR);
                    break;
				case 2:
					insert_node(i,j,best_d,l,P_PRmloop00);
					insert_node(k,best_d-1,P_WB);
					break;
				case 3:
					insert_node(i,j,k,best_d,P_PRmloop00);
					insert_node(best_d+1,l,P_WB);
					break;
			}

		}
		break;

		case P_PRmloop01:
		{
			cand_pos_t i = cur_interval->i;
			cand_pos_t l = cur_interval->j;
			cand_pos_t j = cur_interval->k;
			cand_pos_t k = cur_interval->l;
			if(debug) printf("At P_PRmloop01 with i=%d,l=%d,j=%d,k=%d\n",i,l,j,k);

			if (!(i <= j && j < k-1 && k <= l)){
				std::cerr << "border cases: This should not have happened!, P_PRmloop01" << std::endl;
				exit(EXIT_FAILURE);
			}

			if (i<=0 || j<=0 || k<=0 || l<=0 || i>n || j>n || k> n || l> n){
				std::cerr << "impossible cases: This should not have happened!, P_PRmloop01" << std::endl;
				exit(EXIT_FAILURE);
			}

			energy_t min_energy = INF,temp=INF;
			cand_pos_t best_d=-1;
            min_energy = PRmloop01.get(i,j,k,l-1) + cp_penalty;
            cand_pos_t best_row = 1;

            for(cand_pos_t d = k; d < l; ++d){
                temp = PRmloop00.get(i,j,k,d) + WBP.get(d+1,l);
                if (temp < min_energy){
                    min_energy = temp;
                    best_row = 2;
                    best_d = d;
                }
            }
			switch (best_row) {
                case 1:
                    insert_node(i,l-1,j,k,P_PRmloop01);
                    break;
                case 2:
                    insert_node(best_d+1,l,P_WBP);
                    insert_node(i,best_d,j,k,P_PRmloop00);
                    break;
            }

		}
		break;

		case P_PRmloop10:
		{
			cand_pos_t i = cur_interval->i;
			cand_pos_t l = cur_interval->j;
			cand_pos_t j = cur_interval->k;
			cand_pos_t k = cur_interval->l;
			if(debug) printf("At P_PRmloop10 with i=%d,l=%d,j=%d,k=%d\n",i,l,j,k);

			if (!(i <= j && j < k-1 && k <= l)){
				std::cerr << "border cases: This should not have happened!, P_PRmloop10" << std::endl;
				exit(EXIT_FAILURE);
			}

			if (i<=0 || j<=0 || k<=0 || l<=0 || i>n || j>n || k> n || l> n){
				std::cerr << "impossible cases: This should not have happened!, P_PRmloop10" << std::endl;
				exit(EXIT_FAILURE);
			}

			energy_t min_energy = INF,temp=INF; 
			cand_pos_t best_d=-1;
			min_energy = PRmloop10.get(i,j,k+1,l)+ cp_penalty;
            cand_pos_t best_row = 1;

            for(cand_pos_t d = k+1; d <= l; ++d){
                temp = WBP.get(k,d-1) + PRmloop00.get(i,j,d,l);
                if (temp < min_energy){
                    min_energy = temp;
                    best_row = 2;
                    best_d = d;
                }
            }
			switch (best_row) {
                case 1:
                    insert_node(i,l,j,k+1,P_PRmloop10);
                    break;
                case 2:
                    insert_node(k,best_d-1,P_WBP);
                    insert_node(i,l,j,best_d,P_PRmloop00);
                    break;
            }

		}
		break;

		case P_PMiloop:
		{
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

			f[j].pair = k;
			f[k].pair = j;
			f[j].type = P_PMiloop;
			f[k].type = P_PMiloop;
			
			energy_t min_energy = INF,temp=INF;
			cand_pos_t best_d=-1,best_dp=-1,best_row=-1;

			if (pair[S_[j]][S_[k]]>0){
					// branch 1
				temp = PM.get(i,j-1,k+1,l) + get_e_stP(j-1,k+1);
				if (temp < min_energy){
					min_energy = temp;
					best_row=1;
				}

				for(cand_pos_t d= j-1; d>std::max(i,j-MAXLOOP); --d){
					for (cand_pos_t dp=k+1; dp <std::min(l,k+MAXLOOP); ++dp) {
						temp = get_e_intP(d,j,k,dp) + PM.get(i,d,dp,l);
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

			f[j].pair = k;
			f[k].pair = j;
			f[j].type = P_PMmloop;
			f[k].type = P_PMmloop;

			energy_t branch1 = PMmloop10.get(i,j-1,k+1,l)+ap_penalty+beta2P(j,k);
            energy_t branch2 = PMmloop01.get(i,j-1,k+1,l)+ap_penalty+beta2P(j,k);

            cand_pos_t best_row = -1;
            if (branch1 < branch2) { best_row = 1;}  else { best_row = 2;}

			switch (best_row){
				case 1:
					insert_node(i,l,j-1,k+1,P_PMmloop10);
					break;
				case 2:
					insert_node(i,l,j-1,k+1,P_PMmloop01);
					break;
			}
		}
		break;
		case P_PMmloop00:
		{
			cand_pos_t i = cur_interval->i;
			cand_pos_t l = cur_interval->j;
			cand_pos_t j = cur_interval->k;
			cand_pos_t k = cur_interval->l;
			if(debug) printf("At P_PMmloop00 with i=%d,l=%d,j=%d,k=%d\n",i,l,j,k);

			if (!(i <= j && j < k-1 && k <= l)){
				std::cerr << "border cases: This should not have happened!, P_PMmloop00" << std::endl;
				exit(EXIT_FAILURE);
			}

			if (i<=0 || j<=0 || k<=0 || l<=0 || i>n || j>n || k> n || l> n){
				std::cerr << "impossible cases: This should not have happened!, P_PMmloop00" << std::endl;
				exit(EXIT_FAILURE);
			}

			f[j].pair = k;
			f[k].pair = j;
			f[j].type = P_PMmloop;
			f[k].type = P_PMmloop;

			energy_t temp=INF;
			energy_t min_energy = PM.get(i,j,k,l)+beta2P(j,k);
			cand_pos_t best_row = 1, best_d = -1;

            for(cand_pos_t d=i; d<j; ++d){
                temp=get_WB(d+1,j)+PMmloop00.get(i,d,k,l);
                if (temp < min_energy){
                    min_energy = temp;
                    best_row = 2;
                    best_d = d;
                }
            }
            for(cand_pos_t d=k+1; d<=l; d++){
                temp=PMmloop00.get(i,j,d,l)+get_WB(k,d-1);
                if (temp < min_energy){
                    min_energy = temp;
                    best_row = 3;
                    best_d = d;
                }
            }
			

			switch (best_row){
                case 1:
                	insert_node(i,l,j,k,P_PM);
                	break;
                case 2:
                	insert_node(i,l,best_d,k,P_PMmloop00);
                	insert_node(best_d+1,j,P_WB);
                	break;
                case 3:
                	insert_node(i,l,j,best_d,P_PMmloop00);
                	insert_node(k,best_d-1,P_WB);
                	break;
            }

		}
			break;

		case P_PMmloop01:
		{
			cand_pos_t i = cur_interval->i;
			cand_pos_t l = cur_interval->j;
			cand_pos_t j = cur_interval->k;
			cand_pos_t k = cur_interval->l;
			if(debug) printf("At P_PMmloop01 with i=%d,l=%d,j=%d,k=%d\n",i,l,j,k);

			if (!(i <= j && j < k-1 && k <= l)){
				std::cerr << "border cases: This should not have happened!, P_PMmloop01" << std::endl;
				exit(EXIT_FAILURE);
			}

			if (i<=0 || j<=0 || k<=0 || l<=0 || i>n || j>n || k> n || l> n){
				std::cerr << "impossible cases: This should not have happened!, P_PMmloop01" << std::endl;
				exit(EXIT_FAILURE);
			}

			energy_t temp=INF;
			energy_t min_energy = PMmloop01.get(i,j,k+1,l)+cp_penalty;
			cand_pos_t best_row = 1, best_d = -1;

            for(cand_pos_t d = k+1; d <= l; ++d){
                temp = PMmloop00.get(i,j,d,l) + WBP.get(k,d-1);
                if (temp < min_energy){
                    min_energy = temp;
                    best_row = 2;
                    best_d = d;
                }
            }
			switch (best_row){
                case 1:
                	insert_node(i,l,j,k+1,P_PMmloop01);
                	break;
                case 2:
                	insert_node(i,l,j,best_d,P_PMmloop00);
                	insert_node(k,best_d-1,P_WBP);
                	break;
            }

		}
			break;
		case P_PMmloop10:
		{
			cand_pos_t i = cur_interval->i;
			cand_pos_t l = cur_interval->j;
			cand_pos_t j = cur_interval->k;
			cand_pos_t k = cur_interval->l;
			if(debug) printf("At P_PMmloop10 with i=%d,l=%d,j=%d,k=%d\n",i,l,j,k);

			if (!(i <= j && j < k-1 && k <= l)){
				std::cerr << "border cases: This should not have happened!, P_PMmloop10" << std::endl;
				exit(EXIT_FAILURE);
			}

			if (i<=0 || j<=0 || k<=0 || l<=0 || i>n || j>n || k> n || l> n){
				std::cerr << "impossible cases: This should not have happened!, P_PMmloop10" << std::endl;
				exit(EXIT_FAILURE);
			}

			energy_t temp=INF;
			energy_t min_energy = PMmloop10.get(i,j-1,k,l)+cp_penalty;
			cand_pos_t best_row = 1, best_d = -1;

            for(cand_pos_t d = i+1; d < j; ++d){
                temp = WBP.get(d,j) + PMmloop00.get(i,d-1,k,l);
                if (temp < min_energy){
                    min_energy = temp;
                    best_row = 2;
                    best_d = d;
                }
            }

            switch (best_row){
                case 1:
                	insert_node(i,l,j-1,k,P_PMmloop10);
                	break;
                case 2:
                	insert_node(i,l,best_d-1,k,P_PMmloop00);
                	insert_node(best_d,j,P_WBP);
                	break;
            }

		}
			break;
		case P_POiloop:
		{
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

			f[i].pair = l;
			f[l].pair = i;
			f[i].type = P_POiloop;
			f[l].type = P_POiloop;


			energy_t min_energy = INF;
			cand_pos_t best_d=-1,best_dp=-1,best_row=-1;
			if (pair[S_[i]][S_[l]]>0){
				//branch 1
				energy_t temp = PO.get(i+1,j,k,l-1) + get_e_stP(i,l);
				if(temp < min_energy){
					min_energy = temp;
					best_row=1;
				}
				for(cand_pos_t d= i+1; d<std::min(j,i+MAXLOOP); ++d){
					for (cand_pos_t dp=l-1; dp >std::max(l-MAXLOOP,k); --dp) {
						energy_t branch2 = get_e_intP(i,d,dp,l) + PO.get(d,j,dp,k);
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
			f[i].pair = l;
			f[l].pair = i;
			f[i].type = P_POmloop;
			f[l].type = P_POmloop;

            energy_t branch1 = POmloop10.get(i+1,j,k,l-1) + ap_penalty + beta2P(l,i);
            energy_t branch2 = POmloop01.get(i+1,j,k,l-1) + ap_penalty + beta2P(l,i);

            cand_pos_t best_row = -1;
            if (branch1 < branch2) { best_row = 1;} else { best_row = 2;}

			switch (best_row){
				case 1:
					insert_node(i+1,l-1,j,k,P_POmloop10);
					break;
				case 2:
					insert_node(i+1,l-1,j,k,P_POmloop01);
					break;
			}
		}
		break;

		case P_POmloop00:
		{
			cand_pos_t i = cur_interval->i;
			cand_pos_t l = cur_interval->j;
			cand_pos_t j = cur_interval->k;
			cand_pos_t k = cur_interval->l;
			if(debug) printf("At P_POmloop00 with i=%d,l=%d,j=%d,k=%d\n",i,l,j,k);

			if (!(i <= j && j < k-1 && k <= l)){
				std::cerr << "border cases: This should not have happened!, P_POmloop00" << std::endl;
				exit(EXIT_FAILURE);
			}

			if (i<=0 || j<=0 || k<=0 || l<=0 || i>n || j>n || k> n || l> n){
				std::cerr << "impossible cases: This should not have happened!, P_POmloop00" << std::endl;
				exit(EXIT_FAILURE);
			}

			energy_t min_energy = PO.get(i,j,k,l)+beta2P(l,i),temp=INF;
			cand_pos_t best_row=1,best_d=-1;

			for(cand_pos_t d=i+1; d<=j; ++d){
				temp=get_WB(i,d-1)+POmloop00.get(d,j,k,l);
				if (temp < min_energy){
					min_energy = temp;
					best_row = 2;
					best_d = d;
				}
			}
			for(cand_pos_t d=k; d<l; ++d){
				temp=POmloop00.get(i,j,k,d)+get_WB(d+1,l);
				if (temp < min_energy){
					min_energy = temp;
					best_row = 3;
					best_d = d;
				}
			}

			switch  (best_row){
				case 1:
					insert_node(i,l,j,k,P_PO);
					break;
				case 2:
					insert_node(best_d,l,j,k,P_POmloop00);
					insert_node(i,best_d-1,P_WBP);
					break;
				case 3:
					insert_node(i,best_d,j,k,P_POmloop00);
					insert_node(best_d+1,l,P_WB);
					break;
			}
		}
			break;

		case P_POmloop01:
		{
			cand_pos_t i = cur_interval->i;
			cand_pos_t l = cur_interval->j;
			cand_pos_t j = cur_interval->k;
			cand_pos_t k = cur_interval->l;
			if(debug) printf("At P_POmloop01 with i=%d,l=%d,j=%d,k=%d\n",i,l,j,k);

			if (!(i <= j && j < k-1 && k <= l)){
				std::cerr << "border cases: This should not have happened!, P_POmloop01" << std::endl;
				exit(EXIT_FAILURE);
			}

			if (i<=0 || j<=0 || k<=0 || l<=0 || i>n || j>n || k> n || l> n){
				std::cerr << "impossible cases: This should not have happened!, P_POmloop01" << std::endl;
				exit(EXIT_FAILURE);
			}

			energy_t min_energy = INF,temp=INF;
			cand_pos_t best_d=-1;
            for(cand_pos_t d = k; d < l; d++){
                temp = POmloop00.get(i,j,k,d) + WBP.get(d+1,l);
                if (temp < min_energy){
                    min_energy = temp;
					best_d = d;
                }
            }
			insert_node(i,best_d,j,k,P_POmloop00);
            insert_node(best_d+1,l,P_WBP);
		}
			break;
		case P_POmloop10:
		{
			cand_pos_t i = cur_interval->i;
			cand_pos_t l = cur_interval->j;
			cand_pos_t j = cur_interval->k;
			cand_pos_t k = cur_interval->l;
			if(debug) printf("At P_POmloop10 with i=%d,l=%d,j=%d,k=%d\n",i,l,j,k);

			if (!(i <= j && j < k-1 && k <= l)){
				std::cerr << "border cases: This should not have happened!, P_POmloop10" << std::endl;
				exit(EXIT_FAILURE);
			}

			if (i<=0 || j<=0 || k<=0 || l<=0 || i>n || j>n || k> n || l> n){
				std::cerr << "impossible cases: This should not have happened!, P_POmloop10" << std::endl;
				exit(EXIT_FAILURE);
			}

			energy_t min_energy = INF,temp=INF;
			cand_pos_t best_d=-1,best_row=-1;
			for(cand_pos_t d=i+1; d<=j;++d){
                temp=WBP.get(i,d-1)+POmloop00.get(d,j,k,l);
                if (temp < min_energy){
                    min_energy = temp;
                    best_row = 1;
                    best_d = d;
                }

            }
            for(cand_pos_t d = k+1; d < l; ++d){
                temp = POmloop10.get(i,j,k,d) + get_WB(d+1,l);
                if (temp < min_energy){
                    min_energy = temp;
                    best_row = 2;
                    best_d = d;
                }
            }

			switch  (best_row){
                case 1:
                    insert_node(best_d,l,j,k,P_POmloop00);
                    insert_node(i,best_d-1,P_WBP);
                    break;
                case 2:
                    insert_node(i,best_d,j,k,P_POmloop10);
                    insert_node(best_d+1,l,P_WB);
                    break;
			}
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
    tmp->type = type;
    tmp->next = stack_interval;
    stack_interval = tmp;
}

void pseudo_loop::set_stack_interval(seq_interval *stack_interval){
	this->stack_interval = stack_interval;
}
