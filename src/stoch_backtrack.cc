#include "part_func.hh"
#include "dot_plot.hh"
#include "h_globals.hh"
#include "pf_globals.hh"

#include <algorithm>
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <string>
#include <cassert>

#define debug 0

/*                                BPP                                            */

inline cand_pos_t boustrophedon_at(cand_pos_t start, cand_pos_t end, cand_pos_t pos) {
    cand_pos_t count = pos - 1;
    cand_pos_t advance = (cand_pos_t)(count / 2);

    return start + (end - start) * (count % 2) + advance - (2 * (count % 2)) * advance;
}

std::vector<cand_pos_t> boustrophedon(cand_pos_t start, cand_pos_t end) {
    std::vector<cand_pos_t> seq;

    if (end >= start) {
        seq.emplace_back(end - start + 1);
        for (cand_pos_t pos = 1; pos <= end - start + 1; pos++)
            seq.emplace_back(boustrophedon_at(start, end, pos));
    }

    return seq;
}

void W_final_pf::Sample_W(cand_pos_t start, cand_pos_t end, std::vector<int> &fres) {
    if (debug) printf("W at %d and %d with W[j]=%f,%f\n", start, end, W[end], to_Energy(W[end], end));
    cand_pos_t j = end;
    cand_pos_t m = end;

    pf_t W_temp = 0;
    if (end > start) {
        for (; j > start; --j) {         // Moving through the unpaired bases in j
            W_temp = W[j - 1] * scale[1];
            pf_t r = vrna_urn() * W[j];
            if (r > W_temp) { // Checking if our random sample means j is paired or unpaired
                break;        // j is paired
            }
        }
        if (j <= start + TURN) return; // No more base pairs can occur, but still successful
        pf_t r = vrna_urn() * (W[j] - W_temp);
        std::vector<cand_pos_t> is = boustrophedon(start, j - 1); // applies an alternating list so that the base pairing isn't biased to the right side
        cand_pos_t bous_n = is.size();
        pf_t qt = 0;
        cand_pos_t k = start;
        bool pseudoknot = false;
        for (m = 1; m < bous_n; ++m) {
            k = is[m];
            pf_t acc = (k > 1) ? W[k - 1] : 1;
            pf_t Wkl = acc * V.get(k,j) * exp_Extloop(k, j);
            qt += Wkl;
            if (qt > r) {
                break; // k pairs with j
            }

            Wkl = acc * P.get(k,j) * expPS_penalty;
            qt += Wkl;
            if (qt > r) {
                pseudoknot = true;
                break; // k pairs with j as a pseudoknot
            }
        }
        if (k + start > j) {
            printf("backtracking failed in ext loop at %d and %d with W[j] = %f, qt:%f < r:%f\n", start, end, W[j], qt, r);
            exit(0); /* error */
        }
        Sample_W(start, k-1, fres);
        if (!pseudoknot) {
            Sample_V(k, j, fres);
        } else {
            Sample_P(k, j, fres);
        }
    }
}

void W_final_pf::Sample_V(cand_pos_t i, cand_pos_t j, std::vector<int> &fres) {
    if (debug) printf("V at %d and %d\n", i, j);

    cand_pos_t k = i;
    cand_pos_t l = j;
    fres[i] = j;
    fres[j] = i; 

    pf_t qbr = V.get(i,j);
    pf_t V_temp = 0;

    std::pair<cand_pos_tu, cand_pos_tu> base_pair(i, j);
    std::pair<cand_pos_tu, cand_pos_tu> base_pair_reversed(j, i);
    ++samples[base_pair]; // Increments the base pair found in V
    ++samples[base_pair_reversed];

    pf_t r = vrna_urn() * qbr;
    pf_t qbt1 = 0;

    V_temp = HairpinE(i, j);
    qbt1 += V_temp;
    if (qbt1 >= r) return;

    cand_pos_t max_k = std::min(j - TURN - 2, i + MAXLOOP + 1); // i+1+tree.up[i+1]?
    const pair_type ptype_closing = pair[S_[i]][S_[j]];
    for (k = i + 1; k <= max_k; k++) {
        cand_pos_t min_l = std::max(k + TURN + 1 + MAXLOOP + 2, k + j - i) - MAXLOOP - 2;
        for (l = j - 1; l >= min_l; --l) {
            cand_pos_t u1 = k - i - 1;
            cand_pos_t u2 = j - l - 1;
            V_temp = V.get(k,l)
                        * exp_E_IntLoop(u1, u2, ptype_closing, rtype[pair[S_[k]][S_[l]]], S1_[i + 1], S1_[j - 1], S1_[k - 1], S1_[l + 1],
                                        exp_params_);
            V_temp *= scale[u1 + u2 + 2];
            qbt1 += V_temp;
            if (qbt1 >= r) break;
        }
        if (qbt1 >= r) break;
    }
    if (qbt1 >= r) {
        Sample_V(k, l, fres); // Backtrack the internal loop
        return;
    }

    V_temp = VM.get(i,j); // VM includes everything since it includes the basepair (i.e. not like WM2 region), so is this fine?
    qbt1 += V_temp;
    if (qbt1 < r) {
        printf("Backtracking failed for pair (%d,%d) with qbt1=%f and r=%f and qbr=%f\n", i, j,qbt1,r,qbr);
        exit(0);
    }

    // Must be a multiloop
    Sample_VM(i, j, fres);
}

void W_final_pf::Sample_VM(cand_pos_t i, cand_pos_t j, std::vector<int> &fres) {
    if (debug) printf("VM at %d and %d\n", i, j);
    cand_pos_t k;
    pf_t qt = 0;
    if ((i + 1) + 2 * TURN + 2 >= (j - 1)) {
        printf("backtracking impossible for VM[%d, %d]\n", i, j);
        exit(0); /* error */
    }
    pf_t V_temp = 0.;
    pf_t VM_inside = VM.get(i, j) / scale[2]; // If I remove scale from VM's saved values, I may save time here.
    pf_t r = vrna_urn() * VM_inside;
    bool unpaired = false;
    bool pseudoknot = false;
    for (k = i + 1; k <= j - TURN - 1; ++k) {
        V_temp = WM.get(i + 1, k - 1) * WMv.get(k, j - 1) * exp_Mbloop(i, j) * exp_params_->expMLclosing;
        qt += V_temp;
        if (qt > r) {
            break;
        }

        V_temp = (WM.get(i + 1, k - 1) * WMp.get(k, j - 1) * exp_Mbloop(i, j) * exp_params_->expMLclosing);
        qt += V_temp;
        if (qt > r) {
            pseudoknot = true;
            break;
        }

        V_temp = (expMLbase[k - i - 1] * WMp.get(k, j - 1) * exp_Mbloop(i, j) * exp_params_->expMLclosing);
        qt += V_temp;
        if (qt > r) {
            unpaired = true;
            pseudoknot = true;
            break;
        }
    }
    if (k > j - TURN) {
        printf("backtracking failed for VM at i=%d and j =%d\n", i, j);
        exit(0);
    }

    if (!unpaired) {
        Sample_WM(i + 1, k - 1, fres);
    }
    if (!pseudoknot) { // Case 1
        Sample_WMv(k, j - 1, fres);
    } else { // Case 2 or 3
        Sample_WMp(k, j - 1, fres);
    }
}
void W_final_pf::Sample_WM(cand_pos_t i, cand_pos_t j, std::vector<int> &fres) {
    if (debug) printf("WM at %d and %d\n", i, j);
    cand_pos_t k;
    pf_t qt = 0;
    pf_t qbt1 = 0;
    pf_t qbt2 = 0;
    bool unpaired = false;
    bool pseudoknot = false;

    pf_t V_temp = 0.;

    if (i + TURN >= j) {
        // return;
        printf("backtracking impossible for WM[%u, %u]\n", i, j);
        exit(0); /* error */
    }

    for (; j > i + TURN; --j) {
        pf_t r = vrna_urn() * (WM.get(i, j));

        V_temp = WM.get(i, j - 1) * expMLbase[1];
        qt = V_temp;
        if (r > qt) {
            break;
        }
    }

    if (i + TURN == j) {
        printf("backtracking failed for WM\n");
        exit(0); /* error */
    }

    qt = 0.;
    pf_t qm_rem = WM.get(i, j) - V_temp;
    pf_t r = vrna_urn() * qm_rem;
    for (k = i; k < j - TURN; ++k) {
        qbt1 = V.get(k, j) * exp_MLstem(k, j);
        qbt2 = P.get(k, j) * expPSM_penalty * expb_penalty;

        V_temp = static_cast<pf_t>(expMLbase[k - i]) * qbt1;
        qt += V_temp;
        if (qt >= r) {
            unpaired = true;
            break;
        }

        V_temp = static_cast<pf_t>(expMLbase[k - i]) * qbt2;
        qt += V_temp;
        if (qt >= r) {
            unpaired = true;
            pseudoknot = true;
            break;
        }

		V_temp = WM.get(i, k - 1) * qbt1;
		qt += V_temp;
		if (qt >= r) break;

		V_temp = WM.get(i, k - 1) * qbt2;
		qt += V_temp;
		if (qt >= r) {
			pseudoknot = true;
			break;
		}
        
    }
    if (k > j - TURN || qt < r) {
        printf("backtracking failed for WM at i=%d and j =%d with k=%d, qt=%f and r =%f and qt<r=%d\n", i, j,k,qt,r,qt<r);
        exit(0);
    }
    if (!unpaired) {
        Sample_WM(i, k - 1, fres);
    }
    if (!pseudoknot) {
        Sample_V(k, j, fres);
    } else {
        Sample_P(k, j, fres);
    }
}
void W_final_pf::Sample_WMv(cand_pos_t i, cand_pos_t j, std::vector<int> &fres) {
    if (debug) printf("WMv at %d and %d\n", i, j);
    pf_t qt = 0;

    pf_t V_temp = 0.;

    for (; j > i + TURN; --j) {
        pf_t r = vrna_urn() * WMv.get(i, j);

        V_temp = WMv.get(i, j - 1) * expMLbase[1];
        qt = V_temp;
        if (r > qt) {
            break;
        }
    }

    if (i + TURN == j) {
        printf("backtracking failed for WMV\n");
        exit(0); /* error */
    }

    Sample_V(i, j, fres);
}

void W_final_pf::Sample_WMp(cand_pos_t i, cand_pos_t j, std::vector<int> &fres) {
    if (debug) printf("WMp at %d and %d\n", i, j);
    pf_t qt = 0;

    pf_t V_temp = 0.;

    for (; j > i + TURN; --j) {
        pf_t r = vrna_urn() * WMp.get(i, j);

        V_temp = WMp.get(i, j - 1) * expMLbase[1];
        qt = V_temp;
        if (r > qt) {
            break;
        }
    }

    if (i + TURN == j) { // I'm kinda assuming something like ([..)] for this
        printf("backtracking failed for WMP\n");
        exit(0); /* error */
    }

    Sample_P(i, j, fres);
}
///////////////////// Util /////////////////////////////
void W_final_pf::Sample_WP(cand_pos_t i, cand_pos_t l,std::vector<int> &fres) {
    if (debug) printf("WP at %d and %d\n", i, l);
    assert(!(i<=0 || l<=0 || i>n || l>n));
    if(i>l) return;
    pf_t qt = 0;
    pf_t r = vrna_urn() * calc_WP(i,l);
    qt+=expPUP_pen[l-i+1];
    if(qt>r) return; // If just unpaired bases, stop
    Sample_WPP(i,l,fres);
}
void W_final_pf::Sample_WPP(cand_pos_t i, cand_pos_t l, std::vector<int> &fres) {
    if (debug) printf("WPP at %d and %d\n", i, l);
    pf_t WPP_temp = 0.;
    pf_t qt = 0;

    for (; l > i + TURN; --l) {
        pf_t r = vrna_urn() * WPP.get(i,l);

        WPP_temp = WPP.get(i,l-1) * expMLbase[1];
        qt = WPP_temp;
        if (r > qt) {
            break;
        }
    }
    pf_t r = vrna_urn() * (WPP.get(i,l) - WPP_temp);
    qt = 0;
    for(cand_pos_t d=i; d<l; ++d){
		qt += calc_WP(i,d-1)*V.get(d,l)*gamma2(l,d)*expPPS_penalty;
        if(qt>r){
            Sample_WP(i,d-1,fres);
            Sample_V(d,l,fres);
            return;
        }
		qt += calc_WP(i,d-1)*P.get(d,l)*expPSP_penalty*expPPS_penalty;
        if(qt>r){
            Sample_WP(i,d-1,fres);
            Sample_P(d,l,fres);
            return;
        }
	}
    UNREACHABLE();
}
void W_final_pf::Sample_WB(cand_pos_t i, cand_pos_t l, std::vector<int> &fres) {
    if (debug) printf("WB at %d and %d\n", i, l);
    assert(!(i<=0 || l<=0 || i>n || l>n));
    if(i>l) return;

    pf_t qt = 0;
    pf_t r = vrna_urn() * calc_WB(i,l);
    qt+=expcp_pen[l-i+1];
    if(qt>=r) return; // If just unpaired bases, stop. > or >=?
    Sample_WBP(i,l,fres);
}
void W_final_pf::Sample_WBP(cand_pos_t i, cand_pos_t l, std::vector<int> &fres) {
    if (debug) printf("WBP at %d and %d\n", i, l);
    pf_t WBP_temp = 0.;
    pf_t qt = 0;
    
    for (; l > i + TURN; --l) {
        pf_t r = vrna_urn() * WBP.get(i,l);

        WBP_temp = WBP.get(i,l-1) * expMLbase[1];
        qt = WBP_temp;
        if (r > qt) {
            break;
        }
    }
    pf_t r = vrna_urn() * (WBP.get(i,l) - WBP_temp);
    qt = 0;
    for(cand_pos_t d=i; d< l; ++d){
		qt += calc_WB(i,d-1)*V.get(d,l)*expbp_penalty*expPPS_penalty;
        if(qt>r){
            Sample_WB(i,d-1,fres);
            Sample_V(d,l,fres);
            return;
        }
		qt += calc_WB(i,d-1)*P.get(d,l)*expPSM_penalty*expPPS_penalty;
        if(qt>r){
            Sample_WB(i,d-1,fres);
            Sample_P(d,l,fres);
            return;
        }
	}
    UNREACHABLE();
}
/////////////////// Main ///////////////////////
void W_final_pf::Sample_P(cand_pos_t i, cand_pos_t l, std::vector<int> &fres) {
    if (debug) printf("P at %d and %d\n", i, l);
    pf_t qt = 0;
    pf_t r = vrna_urn() * P.get(i,l);
    for(cand_pos_t j=i; j<l; ++j){
		for (cand_pos_t d=j+1; d<l; ++d){
			for (cand_pos_t k=d+1; k<l; ++k){
				qt += PK1Om.get(i,j,d+1,k)*PK2Om.get(j+1,d,k+1,l);
                if(qt>r){
                    Sample_PK1X(i,j,d+1,k,MType::Om,fres);
                    Sample_PK2X(j+1,d,k+1,l,MType::Om,fres);
                    return;
                }
				qt += PK1Om.get(i,j,d+1,k)*PK2Os.get(j+1,d,k+1,l);
                if(qt>r){
                    Sample_PK1X(i,j,d+1,k,MType::Om,fres);
                    Sample_PK2X(j+1,d,k+1,l,MType::Os,fres);
                    return;
                }
				qt += PK1Om.get(i,j,d+1,k)*PK2LreO.get(j+1,d,k+1,l);
                if(qt>r){
                    Sample_PK1X(i,j,d+1,k,MType::Om,fres);
                    Sample_PK2X(j+1,d,k+1,l,MType::LreO,fres);
                    return;
                }
				qt += PK1Om.get(i,j,d+1,k)*PK2MreO.get(j+1,d,k+1,l);
                if(qt>r){
                    Sample_PK1X(i,j,d+1,k,MType::Om,fres);
                    Sample_PK2X(j+1,d,k+1,l,MType::MreO,fres);
                    return;
                }
				qt += PK1Om.get(i,j,d+1,k)*PK2LreR.get(j+1,d,k+1,l);
                if(qt>r){
                    Sample_PK1X(i,j,d+1,k,MType::Om,fres);
                    Sample_PK2X(j+1,d,k+1,l,MType::LreR,fres);
                    return;
                }
				qt += PK1Om.get(i,j,d+1,k)*PK2MreR.get(j+1,d,k+1,l);
                if(qt>r){
                    Sample_PK1X(i,j,d+1,k,MType::Om,fres);
                    Sample_PK2X(j+1,d,k+1,l,MType::MreR,fres);
                    return;
                }
				qt += PK1Om.get(i,j,d+1,k)*PK2R.get(j+1,d,k+1,l);
                if(qt>r){
                    Sample_PK1X(i,j,d+1,k,MType::Om,fres);
                    Sample_PK2X(j+1,d,k+1,l,MType::R,fres);
                    return;
                }

				qt += PK1Os.get(i,j,d+1,k)*PK2Om.get(j+1,d,k+1,l);
                if(qt>r){
                    Sample_PK1X(i,j,d+1,k,MType::Os,fres);
                    Sample_PK2X(j+1,d,k+1,l,MType::Om,fres);
                    return;
                }
				qt += PK1Os.get(i,j,d+1,k)*PK2Os.get(j+1,d,k+1,l);
                if(qt>r){
                    Sample_PK1X(i,j,d+1,k,MType::Os,fres);
                    Sample_PK2X(j+1,d,k+1,l,MType::Os,fres);
                    return;
                }
				qt += PK1Os.get(i,j,d+1,k)*PK2MreO.get(j+1,d,k+1,l);
                if(qt>r){
                    Sample_PK1X(i,j,d+1,k,MType::Os,fres);
                    Sample_PK2X(j+1,d,k+1,l,MType::MreO,fres);
                    return;
                }
				qt += PK1Os.get(i,j,d+1,k)*PK2MreR.get(j+1,d,k+1,l);
                if(qt>r){
                    Sample_PK1X(i,j,d+1,k,MType::Os,fres);
                    Sample_PK2X(j+1,d,k+1,l,MType::MreR,fres);
                    return;
                }
				qt += PK1Os.get(i,j,d+1,k)*PK2R.get(j+1,d,k+1,l);
                if(qt>r){
                    Sample_PK1X(i,j,d+1,k,MType::Os,fres);
                    Sample_PK2X(j+1,d,k+1,l,MType::R,fres);
                    return;
                }

				qt += PK1LreO.get(i,j,d+1,k)*PK2Om.get(j+1,d,k+1,l);
                if(qt>r){
                    Sample_PK1X(i,j,d+1,k,MType::LreO,fres);
                    Sample_PK2X(j+1,d,k+1,l,MType::Om,fres);
                    return;
                }
				qt += PK1LreO.get(i,j,d+1,k)*PK2LreO.get(j+1,d,k+1,l);
                if(qt>r){
                    Sample_PK1X(i,j,d+1,k,MType::LreO,fres);
                    Sample_PK2X(j+1,d,k+1,l,MType::LreO,fres);
                    return;
                }
				qt += PK1LreO.get(i,j,d+1,k)*PK2MreO.get(j+1,d,k+1,l);
                if(qt>r){
                    Sample_PK1X(i,j,d+1,k,MType::LreO,fres);
                    Sample_PK2X(j+1,d,k+1,l,MType::MreO,fres);
                    return;
                }

				qt += PK1LMreR.get(i,j,d+1,k)*PK2Om.get(j+1,d,k+1,l);
                if(qt>r){
                    Sample_PK1X(i,j,d+1,k,MType::LMreR,fres);
                    Sample_PK2X(j+1,d,k+1,l,MType::Om,fres);
                    return;
                }
				qt += PK1LMreR.get(i,j,d+1,k)*PK2LreO.get(j+1,d,k+1,l);
                if(qt>r){
                    Sample_PK1X(i,j,d+1,k,MType::LMreR,fres);
                    Sample_PK2X(j+1,d,k+1,l,MType::LreO,fres);
                    return;
                }
				qt += PK1LMreR.get(i,j,d+1,k)*PK2MreO.get(j+1,d,k+1,l);
                if(qt>r){
                    Sample_PK1X(i,j,d+1,k,MType::LMreR,fres);
                    Sample_PK2X(j+1,d,k+1,l,MType::MreO,fres);
                    return;
                }

				qt += PK1LMorO.get(i,j,d+1,k)*PK2LreR.get(j+1,d,k+1,l);
                if(qt>r){
                    Sample_PK1X(i,j,d+1,k,MType::LMorO,fres);
                    Sample_PK2X(j+1,d,k+1,l,MType::LreR,fres);
                    return;
                }
				qt += PK1LMorO.get(i,j,d+1,k)*PK2MreR.get(j+1,d,k+1,l);
                if(qt>r){
                    Sample_PK1X(i,j,d+1,k,MType::LMorO,fres);
                    Sample_PK2X(j+1,d,k+1,l,MType::MreR,fres);
                    return;
                }
				qt += PK1LMorO.get(i,j,d+1,k)*PK2R.get(j+1,d,k+1,l);
                if(qt>r){
                    Sample_PK1X(i,j,d+1,k,MType::LMorO,fres);
                    Sample_PK2X(j+1,d,k+1,l,MType::R,fres);
                    return;
                }
			}
		}
	}
    UNREACHABLE();
}
void W_final_pf::Sample_PK1X(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l, MType type, std::vector<int> &fres){
    if (debug) std::cout << "PK1X at " << i << " and " << j << " and " << k << " and " << l << " with type: " << type << std::endl;
    pf_t qt = 0;
    Matrix4DPF &PK1X = PK1X_by_mtype(type);
	Matrix4DPF &PX = PX_by_mtype(type);
    pf_t r = vrna_urn() * PK1X.get(i,j,k,l);
    for(cand_pos_t d = i+1;d<=j;++d){
		qt += PX.get(i,d,k,l)*calc_WP(d+1,j)*expPB_penalty;
        if(qt>r){
            Sample_PX(i,d,k,l,type,fres);
            Sample_WP(d+1,j,fres);
            return;
        }
	}
    UNREACHABLE();
}
void W_final_pf::Sample_PK2X(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l, MType type, std::vector<int> &fres){
    if (debug) std::cout << "PK2X at " << i << " and " << j << " and " << k << " and " << l << " with type: " << type << std::endl;
    pf_t qt = 0;
	Matrix4DPF &PK2X = PK2X_by_mtype(type);
    Matrix4DPF &PK1X = PK1X_by_mtype(type);
    pf_t r = vrna_urn() * PK2X.get(i,j,k,l);
    for(cand_pos_t d = k;d<l;++d){
		qt += PK1X.get(i,j,d,l)*calc_WP(k,d-1);
        if(qt>r){
            Sample_PK1X(i,j,d,l,type,fres);
            Sample_WP(k,d-1,fres);
            return;
        }
	}
    UNREACHABLE();
}
void W_final_pf::Sample_PXmloop(const Index4D &x, MType type,std::vector<int> &fres){
    if (debug) std::cout << "PXmloop at " << x.i() << " and " << x.j() << " and " << x.k() << " and " << x.l() << " with type: " << type << std::endl;
    Index4D xp(x);
	xp.shrink(type);
    Sample_PXmloop00(xp,type,fres); // Should probably have a qt>r as just a sanity check
}

void W_final_pf::Sample_PXiloop(const Index4D &x, MType type,std::vector<int> &fres){
	if (debug) std::cout << "PXiloop at " << x.i() << " and " << x.j() << " and " << x.k() << " and " << x.l() << " with type: " << type << std::endl;
    switch(type) {
    case MType::L: return Sample_PLiloop(x,type,fres);
    case MType::M: return Sample_PMiloop(x,type,fres);
    case MType::R: return Sample_PRiloop(x,type,fres);
    case MType::Om: return Sample_POiloop(x,type,fres);
	case MType::Os: return Sample_POiloop(x,type,fres);
	case MType::LreO: return Sample_PLiloop(x,type,fres);
	case MType::LreR: return Sample_PLiloop(x,type,fres);
	case MType::MreO: return Sample_PMiloop(x,type,fres);
	case MType::MreR: return Sample_PMiloop(x,type,fres);
	case MType::LMreR: return Sample_PLiloop(x,type,fres);
	case MType::LMorO: return Sample_PLiloop(x,type,fres);
    }
    UNREACHABLE();
}

void W_final_pf::Sample_PfromX(const Index4D &x, MType type,std::vector<int> &fres){
    if (debug) std::cout << "PfromX at " << x.i() << " and " << x.j() << " and " << x.k() << " and " << x.l() << " with type: " << type << std::endl;
    switch(type) {
    case MType::L: return Sample_PfromL(x.i(),x.j(),x.k(),x.l(),type,fres);
    case MType::M: return Sample_PfromM(x.i(),x.j(),x.k(),x.l(),type,fres);
    case MType::R: return Sample_PfromR(x.i(),x.j(),x.k(),x.l(),type,fres);
    case MType::Om: return Sample_PfromO(x.i(),x.j(),x.k(),x.l(),type,fres);
	case MType::Os: return;
	case MType::LreO: return Sample_PfromL(x.i(),x.j(),x.k(),x.l(),type,fres);
	case MType::LreR: return Sample_PfromL(x.i(),x.j(),x.k(),x.l(),type,fres);
	case MType::MreO: return Sample_PfromM(x.i(),x.j(),x.k(),x.l(),type,fres);
	case MType::MreR: return Sample_PfromM(x.i(),x.j(),x.k(),x.l(),type,fres);
	case MType::LMreR: return Sample_PfromL(x.i(),x.j(),x.k(),x.l(),type,fres);
	case MType::LMorO: return Sample_PfromL(x.i(),x.j(),x.k(),x.l(),type,fres);
    }
    UNREACHABLE();
}

void W_final_pf::Sample_PfromXprime(cand_pos_t i,cand_pos_t j,cand_pos_t k, cand_pos_t l, MType type,std::vector<int> &fres){
    if (debug) std::cout << "PfromXprime at " << i << " and " << j << " and " << k << " and " << l << " with type: " << type << std::endl;
    switch(type) {
    case MType::L: return Sample_PfromLprime(i,j,k,l,type,fres);
    case MType::M: return Sample_PfromMprime(i,j,k,l,type,fres);
    case MType::R: return Sample_PfromRprime(i,j,k,l,type,fres);
    case MType::Om: return Sample_PfromOprime(i,j,k,l,type,fres);
	case MType::Os: return;
	case MType::LreO: return Sample_PfromLprime(i,j,k,l,type,fres);
	case MType::LreR: return Sample_PfromLprime(i,j,k,l,type,fres);
	case MType::MreO: return Sample_PfromMprime(i,j,k,l,type,fres);
	case MType::MreR: return Sample_PfromMprime(i,j,k,l,type,fres);
	case MType::LMreR: return Sample_PfromLprime(i,j,k,l,type,fres);
	case MType::LMorO: return Sample_PfromLprime(i,j,k,l,type,fres);
    }
    UNREACHABLE();
}

void W_final_pf::Sample_PfromXdoubleprime(cand_pos_t i,cand_pos_t j,cand_pos_t k, cand_pos_t l, MType type,std::vector<int> &fres){
    if (debug) std::cout << "PfromXdoubleprime at " << i << " and " << j << " and " << k << " and " << l << " with type: " << type << std::endl;
    switch(type) {
    case MType::L: return Sample_PfromLdoubleprime(i,j,k,l,fres);
    case MType::M: return Sample_PfromMdoubleprime(i,j,k,l,fres);
    case MType::R: return Sample_PfromRdoubleprime(i,j,k,l,fres);
    case MType::Om: return Sample_PfromOdoubleprime(i,j,k,l,fres);
	case MType::Os: return;
	case MType::LreO: return Sample_PfromLreOdoubleprime(i,j,k,l,fres);
	case MType::LreR: return Sample_PfromLreRdoubleprime(i,j,k,l,fres);
	case MType::MreO: return Sample_PfromMreOdoubleprime(i,j,k,l,fres);
	case MType::MreR: return Sample_PfromMreRdoubleprime(i,j,k,l,fres);
	case MType::LMreR: return Sample_PfromLMreRdoubleprime(i,j,k,l,fres);
	case MType::LMorO: return Sample_PfromLMorOdoubleprime(i,j,k,l,fres);
    }
    UNREACHABLE();
}

void W_final_pf::Sample_PXmloop00(const Index4D &x, MType type,std::vector<int> &fres){
	if (debug) std::cout << "PXmloop00 at " << x.i() << " and " << x.j() << " and " << x.k() << " and " << x.l() << " with type: " << type << std::endl;
    switch(type) {
    case MType::L: return Sample_PLmloop00(x,type,fres);
    case MType::M: return Sample_PMmloop00(x,type,fres);
    case MType::R: return Sample_PRmloop00(x,type,fres);
    case MType::Om: return Sample_POmloop00(x,type,fres);
	case MType::Os: return Sample_POmloop00(x,type,fres);
	case MType::LreO: return Sample_PLmloop00(x,type,fres);
	case MType::LreR: return Sample_PLmloop00(x,type,fres);
	case MType::MreO: return Sample_PMmloop00(x,type,fres);
	case MType::MreR: return Sample_PMmloop00(x,type,fres);
	case MType::LMreR: return Sample_PLmloop00(x,type,fres);
	case MType::LMorO: return Sample_PLmloop00(x,type,fres);
    }
    UNREACHABLE();
}
void W_final_pf::Sample_PXmloop10(const Index4D &x, MType type,std::vector<int> &fres){
	if (debug) std::cout << "PXmloop10 at " << x.i() << " and " << x.j() << " and " << x.k() << " and " << x.l() << " with type: " << type << std::endl;
    switch(type) {
    case MType::L: return Sample_PLmloop10(x,type,fres);
    case MType::M: return Sample_PMmloop10(x,type,fres);
    case MType::R: return Sample_PRmloop10(x,type,fres);
    case MType::Om: return Sample_POmloop10(x,type,fres);
	case MType::Os: return Sample_POmloop10(x,type,fres);
	case MType::LreO: return Sample_PLmloop10(x,type,fres);
	case MType::LreR: return Sample_PLmloop10(x,type,fres);
	case MType::MreO: return Sample_PMmloop10(x,type,fres);
	case MType::MreR: return Sample_PMmloop10(x,type,fres);
	case MType::LMreR: return Sample_PLmloop10(x,type,fres);
	case MType::LMorO: return Sample_PLmloop10(x,type,fres);
    }
    UNREACHABLE();
}
void W_final_pf::Sample_PXmloop01(const Index4D &x, MType type,std::vector<int> &fres){
	if (debug) std::cout << "PXmloop01 at " << x.i() << " and " << x.j() << " and " << x.k() << " and " << x.l() << " with type: " << type << std::endl;
    switch(type) {
    case MType::L: return Sample_PLmloop01(x,type,fres);
    case MType::M: return Sample_PMmloop01(x,type,fres);
    case MType::R: return Sample_PRmloop01(x,type,fres);
    case MType::Om: return Sample_POmloop01(x,type,fres);
	case MType::Os: return Sample_POmloop01(x,type,fres);
	case MType::LreO: return Sample_PLmloop01(x,type,fres);
	case MType::LreR: return Sample_PLmloop01(x,type,fres);
	case MType::MreO: return Sample_PMmloop01(x,type,fres);
	case MType::MreR: return Sample_PMmloop01(x,type,fres);
	case MType::LMreR: return Sample_PLmloop01(x,type,fres);
	case MType::LMorO: return Sample_PLmloop01(x,type,fres);
    }
    UNREACHABLE();
}

void W_final_pf::Sample_PX(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l, MType type, std::vector<int> &fres){
     Matrix4DPF &PX = PX_by_mtype(type);
    if (debug) std::cout << "PX at " << i << " and " << j << " and " << k << " and " << l << " with type: " << type << " and val: " << PX.get(i,j,k,l) << std::endl;
    const Index4D x(i,j,k,l);
    pf_t qt = 0;
    pf_t r = vrna_urn() * PX.get(i,j,k,l);
    const int ptype_closing = pair[S_[x.lend(type)]][S_[x.rend(type)]];
    fres[x.lend(type)] = x.rend(type);
    fres[x.rend(type)] = x.lend(type);
    std::pair<cand_pos_tu, cand_pos_tu> base_pair(x.lend(type), x.rend(type));
    std::pair<cand_pos_tu, cand_pos_tu> base_pair_reversed(x.rend(type), x.lend(type));
    ++samples[base_pair]; // Increments the base pair found in V
    ++samples[base_pair_reversed];
    if (ptype_closing>0){
        if(type == MType::Os){
			if(i==j && k==l){
				qt+=gamma2(x.l(),x.i());
                return; // I think if I get here, I'm just supposed to stop.
			}
		} else if (x.difference(type)>TURN){
			Index4D xp(x);
			xp.shrink(type);
			Matrix4DPF &PfromX = PfromX_by_mtype(type);
			qt += PfromX.get(xp)*penalty(xp, gamma2, type);
            if(qt>r){
                Sample_PfromX(xp,type,fres);
                return;
            }
		}
        qt += calc_PXiloop(x, type);
        if(qt>r){
            Sample_PXiloop(x,type,fres);
            return;
        }
		qt += calc_PXmloop(x,type);
        if(qt>r){
            Sample_PXmloop(x,type,fres);
            return;
        }
    }
    UNREACHABLE();
}
void W_final_pf::Sample_PLiloop(const Index4D &x, MType type,std::vector<int> &fres){
    if (debug) std::cout << "PLiloop at " << x.i() << " and " << x.j() << " and " << x.k() << " and " << x.l() << " with type: " << type << std::endl;
    const cand_pos_t i = x.i(), j = x.j(), k = x.k(), l = x.l();
    pf_t qt = 0;
    Matrix4DPF &PX = PX_by_mtype(type);
    pf_t r = vrna_urn()*calc_PLiloop(x,type); // I have to recalculate this every time I come here because I don't have the value!!!!
    if (i+TURN+2<j) { 
		qt += PX.get(i+1,j-1,k,l)*get_e_stP(i,j);
        if(qt>r){
            Sample_PX(i+1,j-1,k,l,type,fres);
            return;
        }
	}
	cand_pos_t max_d = std::min(j,i+MAXLOOP);
	for(cand_pos_t d= i+1; d<max_d; ++d){
		cand_pos_t min_dp = std::max(d+TURN,j-MAXLOOP);
		for(cand_pos_t dp = j-1; dp > min_dp; --dp){
			if (!(pair[S_[d]][S_[dp]]>0)) continue;
			qt += get_e_intP(i,d,dp,j)*PX.get(d,dp,k,l);
            if(qt>r){
                Sample_PX(d,dp,k,l,type,fres);
                return;
            }
		}
	}
    UNREACHABLE();
}
void W_final_pf::Sample_PMiloop(const Index4D &x, MType type,std::vector<int> &fres){
    if (debug) std::cout << "PMiloop at " << x.i() << " and " << x.j() << " and " << x.k() << " and " << x.l() << " with type: " << type << std::endl;
    const cand_pos_t i = x.i(), j = x.j(), k = x.k(), l = x.l();
    pf_t qt = 0;
    Matrix4DPF &PX = PX_by_mtype(type);
    pf_t r = vrna_urn()*calc_PMiloop(x,type); // I have to recalculate this every time I come here because I don't have the value!!!!
    if (i<j && k<l) {
		qt += PX.get(i,j-1,k+1,l)*get_e_stP(j-1,k+1);
        if(qt>r){
            Sample_PX(i,j-1,k+1,l,type,fres);
            return;
        }
	}
	cand_pos_t max_d = std::max(i,j-MAXLOOP);
	for(cand_pos_t d= j-1; d>max_d; --d){
		cand_pos_t min_dp = std::min(l,k+MAXLOOP); // could switch these here so that we are increasing in the first for like all the others
		for (cand_pos_t dp=k+1; dp <min_dp; ++dp) {
			if (!(pair[S_[d]][S_[dp]]>0)) continue;
			qt += get_e_intP(d,j,k,dp)*PX.get(i,d,dp,l);
            if(qt>r){
                Sample_PX(i,d,dp,l,type,fres);
                return;
            }
		}
	}
    UNREACHABLE();
}
void W_final_pf::Sample_PRiloop(const Index4D &x, MType type,std::vector<int> &fres){
    if (debug) std::cout << "PRiloop at " << x.i() << " and " << x.j() << " and " << x.k() << " and " << x.l() << " with type: " << type << std::endl;
    const cand_pos_t i = x.i(), j = x.j(), k = x.k(), l = x.l();
    pf_t qt = 0;
    Matrix4DPF &PX = PX_by_mtype(type);
    pf_t r = vrna_urn()*calc_PRiloop(x,type); // I have to recalculate this every time I come here because I don't have the value!!!!
    if (k+TURN+2<l) { 
		qt += PX.get(i,j,k+1,l-1)*get_e_stP(k,l);
        if(qt>r){
            Sample_PX(i,j,k+1,l-1,type,fres);
            return;
        }
	}
	cand_pos_t max_d = std::min(l,k+MAXLOOP);
	for(cand_pos_t d= k+1; d<max_d; ++d){
		cand_pos_t min_dp = std::max(d+TURN,l-MAXLOOP);
		for(cand_pos_t dp=l-1; dp > min_dp; --dp){
			if (!(pair[S_[d]][S_[dp]]>0)) continue;
			qt += get_e_intP(k,d,dp,l)*PX.get(i,j,d,dp);
            if(qt>r){
                Sample_PX(i,j,d,dp,type,fres);
                return;
            }
		}
	}
    UNREACHABLE();
}
void W_final_pf::Sample_POiloop(const Index4D &x, MType type,std::vector<int> &fres){
    if (debug) std::cout << "POiloop at " << x.i() << " and " << x.j() << " and " << x.k() << " and " << x.l() << " with type: " << type << std::endl;
    const cand_pos_t i = x.i(), j = x.j(), k = x.k(), l = x.l();
    pf_t qt = 0;
    Matrix4DPF &PX = PX_by_mtype(type);
    pf_t r = vrna_urn()*calc_POiloop(x,type); // I have to recalculate this every time I come here because I don't have the value!!!!
    if (i<j && k<l ) { 
		qt += PX.get(i+1,j,k,l-1)*get_e_stP(i,l);
        if(qt>r){
            Sample_PX(i+1,j,k,l-1,type,fres);
            return;
        }
	}
	cand_pos_t max_d = std::min(j,i+MAXLOOP);
	for(cand_pos_t d= i+1; d<max_d; ++d){
		cand_pos_t min_dp = std::max(l-MAXLOOP,k);
		for (cand_pos_t dp=l-1; dp >min_dp; --dp) {
			if (!(pair[S_[d]][S_[dp]]>0)) continue;
			qt += get_e_intP(i,d,dp,l)*PX.get(d,j,k,dp);
            if(qt>r){
                Sample_PX(d,j,k,dp,type,fres);
                return;
            }
		}
	}
    UNREACHABLE();
}
/**
 * 
 * 
 * 
 */
void W_final_pf::Sample_PLmloop00(const Index4D &x, MType type,std::vector<int> &fres){
    if (debug) std::cout << "PLmloop00 at " << x.i() << " and " << x.j() << " and " << x.k() << " and " << x.l() << " with type: " << type << std::endl;
    const cand_pos_t i = x.i(), j = x.j(), k = x.k(), l = x.l();
    pf_t qt = 0;
    Matrix4DPF &PXmloop00 = PXmloop00_by_mtype(type);
    Matrix4DPF &PXmloop10 = PXmloop10_by_mtype(type);
    Matrix4DPF &PXmloop01 = PXmloop01_by_mtype(type);
    pf_t r = vrna_urn()*PXmloop00.get(i,j,k,l);
    for(cand_pos_t d = i+1; d<=j; ++d){
		qt += calc_WB(i,d-1)*PXmloop10.get(d,j,k,l);
        if(qt>r){
            Index4D xp(d,j,k,l);
            Sample_PXmloop10(xp,type,fres);
            Sample_WB(i,d-1,fres);
            return;
        }
	}
    for(cand_pos_t d = i; d<j; ++d){
		qt += PXmloop01.get(i,d,k,l)*calc_WB(d+1,j);
        if(qt>r){
            Index4D xp(i,d,k,l);
            Sample_PXmloop01(xp,type,fres);
            Sample_WB(d+1,j,fres);
            return;
        }
    }
    UNREACHABLE();
}
void W_final_pf::Sample_PMmloop00(const Index4D &x, MType type,std::vector<int> &fres){
    if (debug) std::cout << "PMmloop00 at " << x.i() << " and " << x.j() << " and " << x.k() << " and " << x.l() << " with type: " << type << std::endl;
    const cand_pos_t i = x.i(), j = x.j(), k = x.k(), l = x.l();
    pf_t qt = 0;
    Matrix4DPF &PXmloop00 = PXmloop00_by_mtype(type);
    Matrix4DPF &PXmloop10 = PXmloop10_by_mtype(type);
    Matrix4DPF &PXmloop01 = PXmloop01_by_mtype(type);
    pf_t r = vrna_urn()*PXmloop00.get(i,j,k,l);
    for(cand_pos_t d=i; d<j; ++d){
        qt += PXmloop10.get(i,d,k,l)*calc_WB(d+1,j);
        if(qt>r){
            Index4D xp(i,d,k,l);
            Sample_PXmloop10(xp,type,fres);
            Sample_WB(d+1,j,fres);
            return;
        }
    }
    for(cand_pos_t d=k+1; d<=l; ++d){
        qt += PXmloop01.get(i,j,d,l)*calc_WB(k,d-1);
        if(qt>r){
            Index4D xp(i,j,d,l);
            Sample_WB(k,d-1,fres);
            Sample_PXmloop01(xp,type,fres);
            return;
        }
    }
    UNREACHABLE();
}
void W_final_pf::Sample_PRmloop00(const Index4D &x, MType type,std::vector<int> &fres){
    if (debug) std::cout << "PRmloop00 at " << x.i() << " and " << x.j() << " and " << x.k() << " and " << x.l() << " with type: " << type << std::endl;
    const cand_pos_t i = x.i(), j = x.j(), k = x.k(), l = x.l();
    pf_t qt = 0;
    Matrix4DPF &PXmloop00 = PXmloop00_by_mtype(type);
    Matrix4DPF &PXmloop10 = PXmloop10_by_mtype(type);
    Matrix4DPF &PXmloop01 = PXmloop01_by_mtype(type);
    pf_t r = vrna_urn()*PXmloop00.get(i,j,k,l);
    for(cand_pos_t d=k+1; d<=l; ++d){
		qt += calc_WB(k,d-1)*PXmloop10.get(i,j,d,l);
        if(qt>r){
            Index4D xp(i,j,d,l);
            Sample_WB(k,d-1,fres);
            Sample_PXmloop10(xp,type,fres);
            return;
        }
	}
	for(cand_pos_t d=k; d<l; ++d){
		qt += PXmloop01.get(i,j,k,d)*calc_WB(d+1,l);
        if(qt>r){
            Index4D xp(i,j,k,d);
            Sample_PXmloop01(xp,type,fres);
            Sample_WB(d+1,l,fres);
            return;
        }
    }
    UNREACHABLE();
}
void W_final_pf::Sample_POmloop00(const Index4D &x, MType type,std::vector<int> &fres){
    if (debug) std::cout << "POmloop00 at " << x.i() << " and " << x.j() << " and " << x.k() << " and " << x.l() << " with type: " << type << std::endl;
    const cand_pos_t i = x.i(), j = x.j(), k = x.k(), l = x.l();
    pf_t qt = 0;
    Matrix4DPF &PXmloop00 = PXmloop00_by_mtype(type);
    Matrix4DPF &PXmloop10 = PXmloop10_by_mtype(type);
    Matrix4DPF &PXmloop01 = PXmloop01_by_mtype(type);
    pf_t r = vrna_urn()*PXmloop00.get(i,j,k,l);
    for(cand_pos_t d=i+1; d<=j; ++d){
        qt += calc_WB(i,d-1)*PXmloop10.get(d,j,k,l);
        if(qt>r){
            Index4D xp(d,j,k,l);
            Sample_WB(i,d-1,fres);
            Sample_PXmloop10(xp,type,fres);
            return;
        }
    }
    for(cand_pos_t d=k; d<l; ++d){
        qt += PXmloop01.get(i,j,k,d)*calc_WB(d+1,l);
        if(qt>r){
            Index4D xp(i,j,k,d);
            Sample_PXmloop01(xp,type,fres);
            Sample_WB(d+1,l,fres);
            return;
        }
    }
    UNREACHABLE();
}
/**
 * 
 * 
 * 
 */
void W_final_pf::Sample_PLmloop10(const Index4D &x, MType type,std::vector<int> &fres){
    if (debug) std::cout << "PLmloop10 at " << x.i() << " and " << x.j() << " and " << x.k() << " and " << x.l() << " with type: " << type << std::endl;
    const cand_pos_t i = x.i(), j = x.j(), k = x.k(), l = x.l();
    pf_t qt = 0;
    Matrix4DPF &PXmloop10 = PXmloop10_by_mtype(type);
    Matrix4DPF &PX = PX_by_mtype(type);
    pf_t r = vrna_urn()*PXmloop10.get(i,j,k,l);
    for(cand_pos_t d = i+1; d <= j; ++d){
        qt += PX.get(i,d,k,l)*calc_WB(d+1,j);
        if(qt>r){
            Sample_PX(i,d,k,l,type,fres);
            Sample_WB(d+1,j,fres);
            return;
        }
    }
    UNREACHABLE();
}
void W_final_pf::Sample_PMmloop10(const Index4D &x, MType type,std::vector<int> &fres){
    if (debug) std::cout << "PMmloop10 at " << x.i() << " and " << x.j() << " and " << x.k() << " and " << x.l() << " with type: " << type << std::endl;
    const cand_pos_t i = x.i(), j = x.j(), k = x.k(), l = x.l();
    pf_t qt = 0;
    Matrix4DPF &PXmloop10 = PXmloop10_by_mtype(type);
    Matrix4DPF &PX = PX_by_mtype(type);
    pf_t r = vrna_urn()*PXmloop10.get(i,j,k,l);
    for(cand_pos_t d = k+1; d <= l; ++d){
        qt += PX.get(i,j,d,l)*calc_WB(k,d-1);
        if(qt>r){
            Sample_PX(i,j,d,l,type,fres);
            Sample_WB(k,d-1,fres);
            return;
        }
    }
    UNREACHABLE();
}
void W_final_pf::Sample_PRmloop10(const Index4D &x, MType type,std::vector<int> &fres){
    if (debug) std::cout << "PRmloop10 at " << x.i() << " and " << x.j() << " and " << x.k() << " and " << x.l() << " with type: " << type << std::endl;
    const cand_pos_t i = x.i(), j = x.j(), k = x.k(), l = x.l();
    pf_t qt = 0;
    Matrix4DPF &PXmloop10 = PXmloop10_by_mtype(type);
    Matrix4DPF &PX = PX_by_mtype(type);
    pf_t r = vrna_urn()*PXmloop10.get(i,j,k,l);
    for(cand_pos_t d = k+1; d <= l; ++d){
        qt += PX.get(i,j,k,d)*calc_WB(d+1,l);
        if(qt>r){
            Sample_PX(i,j,k,d,type,fres);
            Sample_WB(d+1,l,fres);
            return;
        }
    }
    UNREACHABLE();
}
void W_final_pf::Sample_POmloop10(const Index4D &x, MType type,std::vector<int> &fres){
    if (debug) std::cout << "POmloop10 at " << x.i() << " and " << x.j() << " and " << x.k() << " and " << x.l() << " with type: " << type << std::endl;
    const cand_pos_t i = x.i(), j = x.j(), k = x.k(), l = x.l();
    pf_t qt = 0;
    Matrix4DPF &PXmloop10 = PXmloop10_by_mtype(type);
    Matrix4DPF &PX = PX_by_mtype(type);
    pf_t r = vrna_urn()*PXmloop10.get(i,j,k,l);
    for(cand_pos_t d=k+1; d<=l;++d){
        qt += PX.get(i,j,k,d)*calc_WB(d+1,l);
        if(qt>r){
            Sample_PX(i,j,k,d,type,fres);
            Sample_WB(d+1,l,fres);
            return;
        }
    }
    UNREACHABLE();
}
/**
 * 
 * 
 * 
 */
void W_final_pf::Sample_PLmloop01(const Index4D &x, MType type,std::vector<int> &fres){
    if (debug) std::cout << "PLmloop01 at " << x.i() << " and " << x.j() << " and " << x.k() << " and " << x.l() << " with type: " << type << std::endl;
    const cand_pos_t i = x.i(), j = x.j(), k = x.k(), l = x.l();
    pf_t qt = 0;
    Matrix4DPF &PXmloop01 = PXmloop01_by_mtype(type);
    Matrix4DPF &PX = PX_by_mtype(type);
    pf_t r = vrna_urn()*PXmloop01.get(i,j,k,l);
    for(cand_pos_t d = i; d < j; ++d){
        qt += expcp_pen[d-i]*PX.get(d,j,k,l);
        if(qt>r){
            Sample_PX(d,j,k,l,type,fres);
            return;
        }
    }
    UNREACHABLE();
}
void W_final_pf::Sample_PMmloop01(const Index4D &x, MType type,std::vector<int> &fres){
    if (debug) std::cout << "PMmloop01 at " << x.i() << " and " << x.j() << " and " << x.k() << " and " << x.l() << " with type: " << type << std::endl;
    const cand_pos_t i = x.i(), j = x.j(), k = x.k(), l = x.l();
    pf_t qt = 0;
    Matrix4DPF &PXmloop01 = PXmloop01_by_mtype(type);
    Matrix4DPF &PX = PX_by_mtype(type);
    pf_t r = vrna_urn()*PXmloop01.get(i,j,k,l);
    for(cand_pos_t d = i; d < j; ++d){
        qt += expcp_pen[j-d]*PX.get(i,d,k,l);
        if(qt>r){
            Sample_PX(i,d,k,l,type,fres);
            return;
        }
    }
    UNREACHABLE();
}
void W_final_pf::Sample_PRmloop01(const Index4D &x, MType type,std::vector<int> &fres){
    if (debug) std::cout << "PRmloop01 at " << x.i() << " and " << x.j() << " and " << x.k() << " and " << x.l() << " with type: " << type << std::endl;
    const cand_pos_t i = x.i(), j = x.j(), k = x.k(), l = x.l();
    pf_t qt = 0;
    Matrix4DPF &PXmloop01 = PXmloop01_by_mtype(type);
    Matrix4DPF &PX = PX_by_mtype(type);
    pf_t r = vrna_urn()*PXmloop01.get(i,j,k,l);
    for(cand_pos_t d = k; d < l; d++){
        qt += expcp_pen[d-k]*PX.get(i,j,d,l);
        if(qt>r){
            Sample_PX(i,j,d,l,type,fres);
            return;
        }
    }
    UNREACHABLE();
}
void W_final_pf::Sample_POmloop01(const Index4D &x, MType type,std::vector<int> &fres){
    if (debug) std::cout << "POmloop01 at " << x.i() << " and " << x.j() << " and " << x.k() << " and " << x.l() << " with type: " << type << std::endl;
    const cand_pos_t i = x.i(), j = x.j(), k = x.k(), l = x.l();
    pf_t qt = 0;
    Matrix4DPF &PXmloop01 = PXmloop01_by_mtype(type);
    Matrix4DPF &PX = PX_by_mtype(type);
    pf_t r = vrna_urn()*PXmloop01.get(i,j,k,l);
    for(cand_pos_t d = i; d < j; ++d){
        qt += expcp_pen[d-i]*PX.get(d,j,k,l);
        if(qt>r){
            Sample_PX(d,j,k,l,type,fres);
            return;
        }
    }
    UNREACHABLE();
}
/**
 * 
 * 
 * 
 */
void W_final_pf::Sample_PfromL(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l, MType type,std::vector<int> &fres){
	if (debug) std::cout << "PfromL at " << i << " and " << j << " and " << k << " and " << l << " with type: " << type << std::endl;
    pf_t qt = 0;
    Matrix4DPF &PfromX = PfromX_by_mtype(type);
    Matrix4DPF &PfromXprime = PfromXprime_by_mtype(type);
    pf_t r = vrna_urn()*PfromX.get(i,j,k,l);
    for(cand_pos_t d=i+1; d<=j; ++d){
		qt += calc_WP(i,d-1)*PfromXprime.get(d,j,k,l);
        if(qt>r){
            Sample_WP(i,d-1,fres);
            Sample_PfromXprime(d,j,k,l,type,fres);
            return;
        }
	}
    UNREACHABLE();
}
void W_final_pf::Sample_PfromM(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l, MType type,std::vector<int> &fres){
    if (debug) std::cout << "PfromM at " << i << " and " << j << " and " << k << " and " << l << " with type: " << type << std::endl;
    pf_t qt = 0;
    Matrix4DPF &PfromX = PfromX_by_mtype(type);
    Matrix4DPF &PfromXprime = PfromXprime_by_mtype(type);
    pf_t r = vrna_urn()*PfromX.get(i,j,k,l);
    for(cand_pos_t d=i; d<j; ++d){
		qt += PfromXprime.get(i,d,k,l)*calc_WP(d+1,j);
        if(qt>r){
            Sample_PfromXprime(i,d,k,l,type,fres);
            Sample_WP(d+1,j,fres);
            return;
        }
	}
    UNREACHABLE();
}
void W_final_pf::Sample_PfromR(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l, MType type,std::vector<int> &fres){
    if (debug) std::cout << "PfromR at " << i << " and " << j << " and " << k << " and " << l << " with type: " << type << std::endl;
    pf_t qt = 0;
    Matrix4DPF &PfromX = PfromX_by_mtype(type);
    Matrix4DPF &PfromXprime = PfromXprime_by_mtype(type);
    pf_t r = vrna_urn()*PfromX.get(i,j,k,l);
    for(cand_pos_t d=k+1; d<=l; ++d){
		qt += calc_WP(k,d-1)*PfromXprime.get(i,j,d,l);
        if(qt>r){
            Sample_WP(k,d-1,fres);
            Sample_PfromXprime(i,j,d,l,type,fres);
            return;
        }
	}
    UNREACHABLE();
}
void W_final_pf::Sample_PfromO(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l, MType type,std::vector<int> &fres){
    if (debug) std::cout << "PfromO at " << i << " and " << j << " and " << k << " and " << l << " with type: " << type << std::endl;
    pf_t qt = 0;
    Matrix4DPF &PfromX = PfromX_by_mtype(type);
    Matrix4DPF &PfromXprime = PfromXprime_by_mtype(type);
    pf_t r = vrna_urn()*PfromX.get(i,j,k,l);
    for(cand_pos_t d=i+1; d<=j; ++d){
		qt += calc_WP(i,d-1)*PfromXprime.get(d,j,k,l);
        if(qt>r){
            Sample_WP(i,d-1,fres);
            Sample_PfromXprime(d,j,k,l,type,fres);
            return;
        }
	}
    UNREACHABLE();
}
/**
 * 
 * 
 * 
 */
void W_final_pf::Sample_PfromLprime(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l, MType type,std::vector<int> &fres){
    if (debug) std::cout << "PfromLprime at " << i << " and " << j << " and " << k << " and " << l << " with type: " << type << std::endl;
    pf_t qt = 0;
    Matrix4DPF &PfromXprime = PfromXprime_by_mtype(type);
    pf_t r = vrna_urn()*PfromXprime.get(i,j,k,l);
    for(cand_pos_t d=i; d<j; ++d){
		qt += calc_PfromXdoubleprime(i,d,k,l,type)*calc_WP(d+1,j);
        if(qt>r){
            Sample_PfromXdoubleprime(i,d,k,l,type,fres);
            Sample_WP(d+1,j,fres);
            return;
        }
	}
    UNREACHABLE();
}
void W_final_pf::Sample_PfromMprime(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l, MType type,std::vector<int> &fres){
    if (debug) std::cout << "PfromMprime at " << i << " and " << j << " and " << k << " and " << l << " with type: " << type << std::endl;
    pf_t qt = 0;
    Matrix4DPF &PfromXprime = PfromXprime_by_mtype(type);
    pf_t r = vrna_urn()*PfromXprime.get(i,j,k,l);
    for(cand_pos_t d=k+1; d<=l; ++d){
		qt += calc_WP(k,d-1)*calc_PfromXdoubleprime(i,j,d,l,type);
        if(qt>r){
            Sample_WP(k,d-1,fres);
            Sample_PfromXdoubleprime(i,j,d,l,type,fres);
            return;
        }
	}
    UNREACHABLE();
}
void W_final_pf::Sample_PfromRprime(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l, MType type,std::vector<int> &fres){
    if (debug) std::cout << "PfromRprime at " << i << " and " << j << " and " << k << " and " << l << " with type: " << type << std::endl;
    pf_t qt = 0;
    Matrix4DPF &PfromXprime = PfromXprime_by_mtype(type);
    pf_t r = vrna_urn()*PfromXprime.get(i,j,k,l);
    for(cand_pos_t d=k; d<l; ++d){
		qt += calc_PfromXdoubleprime(i,j,k,d,type)*calc_WP(d+1,l);
        if(qt>r){
            Sample_PfromXdoubleprime(i,j,k,d,type,fres);
            Sample_WP(d+1,l,fres);
            return;
        }
	}
    UNREACHABLE();
}
void W_final_pf::Sample_PfromOprime(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l, MType type,std::vector<int> &fres){
    if (debug) std::cout << "PfromOprime at " << i << " and " << j << " and " << k << " and " << l << " with type: " << type << std::endl;
    pf_t qt = 0;
    Matrix4DPF &PfromXprime = PfromXprime_by_mtype(type);
    pf_t r = vrna_urn()*PfromXprime.get(i,j,k,l);
    for(cand_pos_t d=k; d<l; ++d){
		qt += calc_PfromXdoubleprime(i,j,k,d,type)*calc_WP(d+1,l);
        if(qt>r){
            Sample_PfromXdoubleprime(i,j,k,d,type,fres);
            Sample_WP(d+1,l,fres);
            return;
        }
	}
    UNREACHABLE();
}
/**
 * 
 * 
 * 
 */
void W_final_pf::Sample_PfromLdoubleprime(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l,std::vector<int> &fres){
    if (debug) std::cout << "PfromLdoubleprime at " << i << " and " << j << " and " << k << " and " << l << std::endl;
    pf_t qt = 0;
    pf_t r = vrna_urn()*calc_PfromLdoubleprime(i,j,k,l);
    qt += PR.get(i,j,k,l)*gamma2(l,k)*expPB_penalty;
    if(qt>r){
        Sample_PX(i,j,k,l,MType::R,fres);
        return;
    }
	qt += PM.get(i,j,k,l)*gamma2(j,k)*expPB_penalty;
    if(qt>r){
        Sample_PX(i,j,k,l,MType::M,fres);
        return;
    }
	qt += POm.get(i,j,k,l)*gamma2(l,i)*expPB_penalty;
    if(qt>r){
        Sample_PX(i,j,k,l,MType::Om,fres);
        return;
    }
	qt += POs.get(i,j,k,l)*gamma2(l,i)*expPB_penalty;
    if(qt>r){
        Sample_PX(i,j,k,l,MType::Os,fres);
        return;
    }
    UNREACHABLE();
}
void W_final_pf::Sample_PfromMdoubleprime(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l,std::vector<int> &fres){
    if (debug) std::cout << "PfromMdoubleprime at " << i << " and " << j << " and " << k << " and " << l << std::endl;
    pf_t qt = 0;
    pf_t r = vrna_urn()*calc_PfromMdoubleprime(i,j,k,l);
    qt += PL.get(i,j,k,l)*gamma2(j,i)*expPB_penalty;
    if(qt>r){
        Sample_PX(i,j,k,l,MType::L,fres);
        return;
    }
	qt += PR.get(i,j,k,l)*gamma2(l,k)*expPB_penalty;
    if(qt>r){
        Sample_PX(i,j,k,l,MType::R,fres);
        return;
    }
	qt += POm.get(i,j,k,l)*gamma2(l,i)*expPB_penalty;
    if(qt>r){
        Sample_PX(i,j,k,l,MType::Om,fres);
        return;
    }
    UNREACHABLE();
}
void W_final_pf::Sample_PfromRdoubleprime(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l,std::vector<int> &fres){
    if (debug) std::cout << "PfromRdoubleprime at " << i << " and " << j << " and " << k << " and " << l << std::endl;
    pf_t qt = 0;
    pf_t r = vrna_urn()*calc_PfromRdoubleprime(i,j,k,l);
    qt += PM.get(i,j,k,l)*gamma2(j,k)*expPB_penalty;
    if(qt>r){
        Sample_PX(i,j,k,l,MType::M,fres);
        return;
    }
	qt += POm.get(i,j,k,l)*gamma2(l,i)*expPB_penalty;
    if(qt>r){
        Sample_PX(i,j,k,l,MType::Om,fres);
        return;
    }
	qt += POs.get(i,j,k,l)*gamma2(l,i)*expPB_penalty;
    if(qt>r){
        Sample_PX(i,j,k,l,MType::Os,fres);
        return;
    }
    UNREACHABLE();
}
void W_final_pf::Sample_PfromOdoubleprime(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l,std::vector<int> &fres){
    if (debug) std::cout << "PfromOdoubleprime at " << i << " and " << j << " and " << k << " and " << l << std::endl;
    pf_t qt = 0;
    pf_t r = vrna_urn()*calc_PfromOdoubleprime(i,j,k,l);
    qt += PL.get(i,j,k,l)*gamma2(j,i)*expPB_penalty;
    if(qt>r){
        Sample_PX(i,j,k,l,MType::L,fres);
        return;
    }
	qt += PR.get(i,j,k,l)*gamma2(l,k)*expPB_penalty;
    if(qt>r){
        Sample_PX(i,j,k,l,MType::R,fres);
        return;
    }
    UNREACHABLE();
}
void W_final_pf::Sample_PfromLreOdoubleprime(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l,std::vector<int> &fres){
    if (debug) std::cout << "PfromLreOdoubleprime at " << i << " and " << j << " and " << k << " and " << l << std::endl;
    pf_t qt = 0;
    pf_t r = vrna_urn()*calc_PfromLreOdoubleprime(i,j,k,l);
    qt += PMreO.get(i,j,k,l)*gamma2(j,k)*expPB_penalty;
    if(qt>r){
        Sample_PX(i,j,k,l,MType::MreO,fres);
        return;
    }
	qt += POm.get(i,j,k,l)*gamma2(l,i)*expPB_penalty;
    if(qt>r){
        Sample_PX(i,j,k,l,MType::Om,fres);
        return;
    }
	qt += POs.get(i,j,k,l)*gamma2(l,i)*expPB_penalty;
    if(qt>r){
        Sample_PX(i,j,k,l,MType::Os,fres);
        return;
    }
    UNREACHABLE();
}
void W_final_pf::Sample_PfromLreRdoubleprime(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l,std::vector<int> &fres){
    if (debug) std::cout << "PfromLreRdoubleprime at " << i << " and " << j << " and " << k << " and " << l << std::endl;
    pf_t qt = 0;
    pf_t r = vrna_urn()*calc_PfromLreRdoubleprime(i,j,k,l);
    qt += PMreR.get(i,j,k,l)*gamma2(j,k)*expPB_penalty;
    if(qt>r){
        Sample_PX(i,j,k,l,MType::MreR,fres);
        return;
    }
	qt += PR.get(i,j,k,l)*gamma2(l,k)*expPB_penalty;
    if(qt>r){
        Sample_PX(i,j,k,l,MType::R,fres);
        return;
    }
    UNREACHABLE();
}
void W_final_pf::Sample_PfromMreOdoubleprime(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l,std::vector<int> &fres){
    if (debug) std::cout << "PfromMreOdoubleprime at " << i << " and " << j << " and " << k << " and " << l << std::endl;
    pf_t qt = 0;
    pf_t r = vrna_urn()*calc_PfromMreOdoubleprime(i,j,k,l);
    qt += PLreO.get(i,j,k,l)*gamma2(j,i)*expPB_penalty;
    if(qt>r){
        Sample_PX(i,j,k,l,MType::LreO,fres);
        return;
    }
	qt += POm.get(i,j,k,l)*gamma2(l,i)*expPB_penalty;
    if(qt>r){
        Sample_PX(i,j,k,l,MType::Om,fres);
        return;
    }
	qt += POs.get(i,j,k,l)*gamma2(l,i)*expPB_penalty;
    if(qt>r){
        Sample_PX(i,j,k,l,MType::Os,fres);
        return;
    }
    UNREACHABLE();
}
void W_final_pf::Sample_PfromMreRdoubleprime(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l,std::vector<int> &fres){
    if (debug) std::cout << "PfromMreRdoubleprime at " << i << " and " << j << " and " << k << " and " << l << std::endl;
    pf_t qt = 0;
    pf_t r = vrna_urn()*calc_PfromMreRdoubleprime(i,j,k,l);
    qt += PLreR.get(i,j,k,l)*gamma2(j,i)*expPB_penalty;
    if(qt>r){
        Sample_PX(i,j,k,l,MType::LreR,fres);
        return;
    }
	qt += PR.get(i,j,k,l)*gamma2(l,k)*expPB_penalty;
    if(qt>r){
        Sample_PX(i,j,k,l,MType::R,fres);
        return;
    }
    UNREACHABLE();
}
void W_final_pf::Sample_PfromLMreRdoubleprime(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l,std::vector<int> &fres){
    if (debug) std::cout << "PfromLMreRdoubleprime at " << i << " and " << j << " and " << k << " and " << l << std::endl;
    pf_t qt = 0;
    pf_t r = vrna_urn()*calc_PfromLMreRdoubleprime(i,j,k,l);
    qt += PMreR.get(i,j,k,l)*gamma2(j,k)*expPB_penalty;
    if(qt>r){
        Sample_PX(i,j,k,l,MType::MreR,fres);
        return;
    }
    UNREACHABLE();
}
void W_final_pf::Sample_PfromLMorOdoubleprime(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l,std::vector<int> &fres){
    if (debug) std::cout << "PfromLMorOdoubleprime at " << i << " and " << j << " and " << k << " and " << l << std::endl;
    pf_t qt = 0;
    pf_t r = vrna_urn()*calc_PfromLMorOdoubleprime(i,j,k,l);
    qt += PM.get(i,j,k,l)*gamma2(j,k)*expPB_penalty;
    if(qt>r){
        Sample_PX(i,j,k,l,MType::M,fres);
        return;
    }
	qt += POm.get(i,j,k,l)*gamma2(l,i)*expPB_penalty;
    if(qt>r){
        Sample_PX(i,j,k,l,MType::Om,fres);
        return;
    }
	qt += POs.get(i,j,k,l)*gamma2(l,i)*expPB_penalty;
    if(qt>r){
        Sample_PX(i,j,k,l,MType::Os,fres);
        return;
    }
    UNREACHABLE();
}
