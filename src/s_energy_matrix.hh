#ifndef ENERGY_MATRIX_H
#define ENERGY_MATRIX_H

#include "base_types.hh"
#include "sparse_tree.hh"
#include <string>
#include <vector>

extern "C" {
#include "ViennaRNA/pair_mat.h"
#include "ViennaRNA/loops/all.h"
#include "ViennaRNA/params/io.h"
}


class s_energy_matrix
{
    public:

        friend class s_multi_loop;

        s_energy_matrix (std::string seq, cand_pos_t length, short *S, short *S1, vrna_param_t *params);
        // The constructor

        ~s_energy_matrix ();
        // The destructor

        vrna_param_t *params_;

        short *S_;
        short *S1_;

        void compute_energy (cand_pos_t i, cand_pos_t j);


        free_energy_node* get_node (cand_pos_t i, cand_pos_t j) { cand_pos_t ij = index[i]+j-i; return &nodes[ij]; }
        // return the node at (i,j)

        energy_t get_energy (cand_pos_t i, cand_pos_t j) { if (i>=j) return INF; cand_pos_t ij = index[i]+j-i; return nodes[ij].energy; }

        energy_t get_energy_WM (cand_pos_t i, cand_pos_t j) { if (i>=j) return INF; cand_pos_t ij = index[i]+j-i; return WM[ij]; }
        energy_t get_energy_WMv (cand_pos_t i, cand_pos_t j) { if (i>=j) return INF; cand_pos_t ij = index[i]+j-i; return WMv[ij]; }
        energy_t get_energy_WMp (cand_pos_t i, cand_pos_t j) { if (i>=j) return INF; cand_pos_t ij = index[i]+j-i; return WMp[ij]; }

        char get_type (cand_pos_t i, cand_pos_t j) { cand_pos_t ij = index[i]+j-i; return nodes[ij].type; }

        energy_t HairpinE(const std::string& seq, const short* S, const short* S1,  const paramT* params, cand_pos_t i, cand_pos_t j);
        energy_t compute_stack(cand_pos_t i, cand_pos_t j, const paramT *params);
        energy_t compute_internal(cand_pos_t i, cand_pos_t j, const paramT *params);
        energy_t compute_int(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l, const paramT *params);

        void compute_energy_WM (cand_pos_t i, cand_pos_t j,std::vector<energy_t> &WMB);
        energy_t compute_energy_VM (cand_pos_t i, cand_pos_t j);
        energy_t E_MLStem(const energy_t& vij,const energy_t& vi1j,const energy_t& vij1,const energy_t& vi1j1,const short* S, paramT* params,cand_pos_t i, cand_pos_t j, cand_pos_t n);
        energy_t E_MbLoop(const energy_t WM2ij, const energy_t WM2ip1j, const energy_t WM2ijm1, const energy_t WM2ip1jm1, const short* S, paramT* params, cand_pos_t i, cand_pos_t j);
        void compute_WMv_WMp(cand_pos_t i, cand_pos_t j, energy_t WMB);

    protected:

        std::vector<energy_t> WM;
        std::vector<energy_t> WMv;
        std::vector<energy_t> WMp;

       
        std::string seq_;
        cand_pos_t n;              // sequence length
        std::vector<cand_pos_t> index;
        // int *index;                // an array with indexes, such that we don't work with a 2D array, but with a 1D array of length (n*(n+1))/2
        std::vector<free_energy_node> nodes;   // the free energy and type (i.e. base pair closing a hairpin loops, stacked pair etc), for each i and j
};



#endif
