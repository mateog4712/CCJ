#include "part_func.hh"
// #include "dot_plot.hh"
#include "h_externs.hh"
#include "pf_externs.hh"

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
        seq.push_back(end - start + 1);
        for (cand_pos_t pos = 1; pos <= end - start + 1; pos++)
            seq.push_back(boustrophedon_at(start, end, pos));
    }

    return seq;
}

void W_final_pf::Sample_W(cand_pos_t start, cand_pos_t end, std::string &structure,
                          std::unordered_map<std::pair<cand_pos_t, cand_pos_t>, cand_pos_t, SzudzikHash> &samples) {
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
        Sample_W(start, k-1, structure, samples);
        if (!pseudoknot) {
            Sample_V(k, j, structure, samples);
        } else {
            Sample_P(k, j, structure, samples);
        }
    }
}

void W_final_pf::Sample_V(cand_pos_t i, cand_pos_t j, std::string &structure,
                          std::unordered_map<std::pair<cand_pos_t, cand_pos_t>, cand_pos_t, SzudzikHash> &samples) {
    if (debug) printf("V at %d and %d\n", i, j);

    cand_pos_t k = i;
    cand_pos_t l = j;
    structure[i - 1] = '(';
    structure[j - 1] = ')';

    pf_t qbr = V.get(i, j);
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
        Sample_V(k, l, structure, samples); // Backtrack the internal loop
        return;
    }

    V_temp = VM.get(i,j); // VM includes everything since it includes the basepair (i.e. not like WM2 region), so is this fine?
    qbt1 += V_temp;
    if (qbt1 < r) {
        printf("Backtracking failed for pair (%d,%d)\n", i, j);
        exit(0);
    }

    // Must be a multiloop
    Sample_VM(i, j, structure, samples);
}

void W_final_pf::Sample_VM(cand_pos_t i, cand_pos_t j, std::string &structure,
                           std::unordered_map<std::pair<cand_pos_t, cand_pos_t>, cand_pos_t, SzudzikHash> &samples) {
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
        Sample_WM(i + 1, k - 1, structure, samples);
    }
    if (!pseudoknot) { // Case 1
        Sample_WMv(k, j - 1, structure, samples);
    } else { // Case 2 or 3
        Sample_WMp(k, j - 1, structure, samples);
    }
}
void W_final_pf::Sample_WM(cand_pos_t i, cand_pos_t j, std::string &structure,
                           std::unordered_map<std::pair<cand_pos_t, cand_pos_t>, cand_pos_t, SzudzikHash> &samples) {
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
        Sample_WM(i, k - 1, structure, samples);
    }
    if (!pseudoknot) {
        Sample_V(k, j, structure, samples);
    } else {
        Sample_P(k, j, structure, samples);
    }
}
void W_final_pf::Sample_WMv(cand_pos_t i, cand_pos_t j, std::string &structure,
                            std::unordered_map<std::pair<cand_pos_t, cand_pos_t>, cand_pos_t, SzudzikHash> &samples) {
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

    Sample_V(i, j, structure, samples);
}

void W_final_pf::Sample_WMp(cand_pos_t i, cand_pos_t j, std::string &structure,
                            std::unordered_map<std::pair<cand_pos_t, cand_pos_t>, cand_pos_t, SzudzikHash> &samples) {
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

    Sample_P(i, j, structure, samples);
}

void W_final_pf::Sample_P(cand_pos_t i, cand_pos_t j, std::string &structure,
                            std::unordered_map<std::pair<cand_pos_t, cand_pos_t>, cand_pos_t, SzudzikHash> &samples) {

}