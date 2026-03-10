#ifndef MATRICES_H
#define MATRICES_H

#include "base_types.hh"
#include <vector>
#include <iostream>
#include <limits>
#include <cassert>
#define INF 10000000 /* (INT_MAX/10) */

typedef std::vector<std::vector<std::vector<size_t>>> index_offset_t;

class TriangleMatrix {
public:
    TriangleMatrix() {}

    void init(int n, std::vector<cand_pos_t> &index, int return_val=INF){
        n_ = n;
        index_ = index;
        return_val_ = return_val;

        size_t tl=total_length(n);
        m_.clear();
        m_.resize(tl,INF+1);
    }

    int ij(int i, int j) const {
        return index_[i]+j-i;
    }

    //! unchecked get
    int get_uc(int i, int j) const {
        assert (i<=j);
        return m_[ij(i,j)];
    }

    int get (int i, int j) const {
        if (i>j) return return_val_;
        return get_uc(i,j);
    }


    int& operator [] (int ij) {
        return m_[ij];
    }

    int& set (int i, int j) {
        return m_[ij(i,j)];
    }

    void print() {
        for (int i=0; i<n_; i++) {
            std::cout << i << ": ";
            for (int j=i; j<n_; j++) {
                std::cout << get(i,j) << " ";
            }
            std::cout << std::endl;
        }
    }

    static int total_length(int n) {
        return (n *(n+1))/2;
    }

    static void new_index(std::vector<cand_pos_t> &index, int n) {
        index.resize(n);
        index[1] = 0;
        for (int i=2; i < n; i++) {
            index[i] = index[i-1]+n-i+1;
        }
    }

private:
    cand_pos_t n_;
    energy_t return_val_;
    std::vector<cand_pos_t> index_;
    std::vector<energy_t> m_;
};

class TriangleMatrix_PF {
public:
    TriangleMatrix_PF() {}

    void init(int n, std::vector<cand_pos_t> &index, pf_t return_val=0.0){
        n_ = n;
        index_ = index;
        return_val_ = return_val;

        size_t tl=total_length(n);
        m_.clear();
        m_.resize(tl,0);
    }

    int ij(int i, int j) const {
        return index_[i]+j-i;
    }

    //! unchecked get
    pf_t get_uc(int i, int j) const {
        assert (i<=j);
        return m_[ij(i,j)];
    }

    pf_t get (int i, int j) const {
        if (i>j) return return_val_;
        return get_uc(i,j);
    }


    pf_t& operator [] (int ij) {
        return m_[ij];
    }

    pf_t& set (int i, int j) {
        return m_[ij(i,j)];
    }

    void print() {
        for (int i=0; i<n_; i++) {
            std::cout << i << ": ";
            for (int j=i; j<n_; j++) {
                std::cout << get(i,j) << " ";
            }
            std::cout << std::endl;
        }
    }

    static int total_length(int n) {
        return (n *(n+1))/2;
    }

    static void new_index(std::vector<cand_pos_t> &index, int n) {
        index.resize(n);
        index[1] = 0;
        for (int i=2; i < n; i++) {
            index[i] = index[i-1]+n-i+1;
        }
    }

private:
    cand_pos_t n_;
    pf_t return_val_;
    std::vector<cand_pos_t> index_;
    std::vector<pf_t> m_;
};

class Matrix4D {
public:
    const energy_16t INTERN_INF = std::numeric_limits<energy_16t>::max();

    //! construct empty
    Matrix4D() {}

    void
    init(cand_pos_t n, index_offset_t &index3D) {
        n_=n;
        index3D_ = index3D;
        int slice_size_ = index3D_[n-1][n-1][n-1] + 1;
        assert( slice_size_ == n*(n+1)*(n+2)*(n+3)/24);
        m_.clear();
        m_.resize(slice_size_,INTERN_INF);

    }

    int get_uc(int i, int j, int k, int l) const {
        assert(!(i<=0 || l> n_));

        int val = m_[index(i,j,k,l)];
        return val;
    }

    int get(int i, int j, int k, int l) const {
        if (!(i <= j && j < k-1 && k <= l)){
            return INF;
        }
        return get_uc(i,j,k,l);
    }

    void set(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l, energy_t e) {
        if (e >= INTERN_INF) e=INTERN_INF; // If I remove this, I need to set the m type to be 32 bit int instead of 16 bit int
        m_[index(i,j,k,l)] = e;
    }

    static void construct_index(index_offset_t &index, cand_pos_t n){
        // Construct for i<=j<=k<=l (even if j<k-1)
        index.resize(n);
        size_t idx = 0;
        for (cand_pos_t i = 0; i < n; ++i) {
            index[i].resize(n);
            for (cand_pos_t j = i; j < n; ++j) {
                index[i][j].resize(n);
                for (cand_pos_t k = j; k < n; ++k) {
                    index[i][j][k] = idx;
                    idx += (n - k);  // remaining values from k to n
                }
            }
        }
    }

    private:
    cand_pos_t n_;
    std::vector<energy_16t> m_; // I had m as a 16 bit int. This seemed to result in some kind of of overflow or underflow, most likely because I allowed pushing INF(32 bit int/10) despite being 16 bit
    index_offset_t index3D_;

    size_t index(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l) const {
        return index3D_[i-1][j-1][k-1] + (l-k);
    }
};

class Matrix4DPF {
public:

    //! construct empty
    Matrix4DPF() {}

    void
    init(cand_pos_t n, index_offset_t &index3D) {
        n_=n;
        index3D_ = index3D;
        int slice_size_ = index3D_[n-1][n-1][n-1] + 1;
        assert( slice_size_ == n*(n+1)*(n+2)*(n+3)/24);
        m_.clear();
        m_.resize(slice_size_,0.0);

    }

    int get_uc(int i, int j, int k, int l) const {
        assert(!(i<=0 || l> n_));

        int val = m_[index(i,j,k,l)];
        return val;
    }

    int get(int i, int j, int k, int l) const {
        if (!(i <= j && j < k-1 && k <= l)){
            return 0.0;
        }
        return get_uc(i,j,k,l);
    }

    void set(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l, energy_t e) {
        m_[index(i,j,k,l)] = e;
    }

    static void construct_index(index_offset_t &index, cand_pos_t n){
        // Construct for i<=j<=k<=l (even if j<k-1)
        index.resize(n);
        size_t idx = 0;
        for (cand_pos_t i = 0; i < n; ++i) {
            index[i].resize(n);
            for (cand_pos_t j = i; j < n; ++j) {
                index[i][j].resize(n);
                for (cand_pos_t k = j; k < n; ++k) {
                    index[i][j][k] = idx;
                    idx += (n - k);  // remaining values from k to n
                }
            }
        }
    }

    private:
    cand_pos_t n_;
    std::vector<pf_t> m_; // I had m as a 16 bit int. This seemed to result in some kind of of overflow or underflow, most likely because I allowed pushing INF(32 bit int/10) despite being 16 bit
    index_offset_t index3D_;

    size_t index(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l) const {
        return index3D_[i-1][j-1][k-1] + (l-k);
    }
};

#endif