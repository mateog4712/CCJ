#ifndef MATRICES_H
#define MATRICES_H

#include "base_types.hh"
#include "index4D.hh"
#include <vector>
#include <iostream>
#include <limits>
#include <cassert>
#define INF 10000000 /* (INT_MAX/10) */

typedef std::vector<size_t> index_offset_t;

class TriangleMatrix {
public:
    TriangleMatrix() {}

    void init(cand_pos_t n, std::vector<cand_pos_t> &index, energy_t return_val=INF){
        n_ = n;
        index_ = &index;
        return_val_ = return_val;

        size_t tl=total_length(n);
        m_.clear();
        m_.resize(tl,INF+1);
    }

    cand_pos_t ij(cand_pos_t i, cand_pos_t j) const {
        return (*index_)[i]+j-i;
    }

    //! unchecked get
    energy_t get_uc(cand_pos_t i, cand_pos_t j) const {
        assert (i<=j);
        return m_[ij(i,j)];
    }

    energy_t get (cand_pos_t i, cand_pos_t j) const {
        if (i>j) return return_val_;
        return get_uc(i,j);
    }


    energy_t& operator [] (cand_pos_t ij) {
        return m_[ij];
    }

    energy_t& set (cand_pos_t i, cand_pos_t j) {
        return m_[ij(i,j)];
    }

    void print() {
        for (cand_pos_t i=0; i<n_; i++) {
            std::cout << i << ": ";
            for (cand_pos_t j=i; j<n_; j++) {
                std::cout << get(i,j) << " ";
            }
            std::cout << std::endl;
        }
    }

    static cand_pos_t total_length(cand_pos_t n) {
        return (n *(n+1))/2;
    }

    static void new_index(std::vector<cand_pos_t> &index, cand_pos_t n) {
        index.resize(n);
        index[1] = 0;
        for (cand_pos_t i=2; i < n; i++) {
            index[i] = index[i-1]+n-i+1;
        }
    }

private:
    cand_pos_t n_;
    energy_t return_val_;
    std::vector<cand_pos_t> *index_;
    std::vector<energy_t> m_;
};

class TriangleMatrix_PF {
public:
    TriangleMatrix_PF() {}

    void init(cand_pos_t n, std::vector<cand_pos_t> &index, pf_t return_val=0.0){
        n_ = n;
        index_ = &index;
        return_val_ = return_val;

        size_t tl=total_length(n);
        m_.clear();
        m_.resize(tl,0);
    }

    cand_pos_t ij(cand_pos_t i, cand_pos_t j) const {
        return (*index_)[i]+j-i;
    }

    //! unchecked get
    pf_t get_uc(cand_pos_t i, cand_pos_t j) const {
        assert (i<=j);
        return m_[ij(i,j)];
    }

    pf_t get (cand_pos_t i, cand_pos_t j) const {
        if (i>j) return return_val_;
        return get_uc(i,j);
    }


    pf_t& operator [] (cand_pos_t ij) {
        return m_[ij];
    }

    pf_t& set (cand_pos_t i, cand_pos_t j) {
        return m_[ij(i,j)];
    }

    void print() {
        for (cand_pos_t i=0; i<n_; i++) {
            std::cout << i << ": ";
            for (cand_pos_t j=i; j<n_; j++) {
                std::cout << get(i,j) << " ";
            }
            std::cout << std::endl;
        }
    }

    static cand_pos_t total_length(int n) {
        return (n *(n+1))/2;
    }

    static void new_index(std::vector<cand_pos_t> &index, cand_pos_t n) {
        index.resize(n);
        index[1] = 0;
        for (cand_pos_t i=2; i < n; i++) {
            index[i] = index[i-1]+n-i+1;
        }
    }

private:
    cand_pos_t n_;
    pf_t return_val_;
    std::vector<cand_pos_t> *index_;
    std::vector<pf_t> m_;
};

class Matrix4D {
public:
    const energy_16t INTERN_INF = std::numeric_limits<energy_16t>::max();

    //! construct empty
    Matrix4D() {}

    void init(cand_pos_t n, index_offset_t &index3D) {
        n_=n;
        nn_ = n*n;
        index3D_ = &index3D;
        int slice_size_ = (*index3D_)[(n-1)*nn_ + (n-1)*n + n-1] + 1;
        assert( slice_size_ == n*(n+1)*(n+2)*(n+3)/24);
        m_.clear();
        m_.resize(slice_size_,INTERN_INF);

    }

    int get_uc(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l) const {
        assert(!(i<=0 || l> n_));

        int val = m_[index(i,j,k,l)];
        return val;
    }

    int get_uc(const Index4D &x) const {
        return get_uc(x.i(),x.j(),x.k(),x.l());
    }

    int get(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l) const {
        if (!(i <= j && j < k-1 && k <= l)){
            return INF;
        }
        return get_uc(i,j,k,l);
    }

    int get(const Index4D &x) const {
        return get(x.i(),x.j(),x.k(),x.l());
    }

    void set(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l, energy_t e) {
        if (e >= INTERN_INF) e=INTERN_INF; // If I remove this, I need to set the m type to be 32 bit int instead of 16 bit int
        m_[index(i,j,k,l)] = e;
    }

    void set(const Index4D &x, int e) {
        set(x.i(),x.j(),x.k(),x.l(), e);
    }

    void setI(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l, int e) {
        if (e >= INTERN_INF) e=INTERN_INF;
        assert( e >= std::numeric_limits<energy_t>::min() );

        m_[index(i,j,k,l)] = e;
    }

    void setI(const Index4D &x, int e) {
        setI(x.i(),x.j(),x.k(),x.l(),e);
    }

    static void construct_index(index_offset_t &index, cand_pos_t n){
        // Construct for i<=j<=k<=l (even if j<k-1)
        index.resize(n*n*n);
        cand_pos_t nn = n*n;
        size_t idx = 0;
        for (cand_pos_t i = 0; i < n; ++i) {
            for (cand_pos_t j = i; j < n; ++j) {
                for (cand_pos_t k = j; k < n; ++k) {
                    index[i*nn+j*n+k] = idx;
                    idx += (n - k);  // remaining values from k to n
                }
            }
        }
    }

    private:
    cand_pos_t n_;
    cand_pos_t nn_;
    std::vector<energy_16t> m_; // I had m as a 16 bit int. This seemed to result in some kind of of overflow or underflow, most likely because I allowed pushing INF(32 bit int/10) despite being 16 bit
    index_offset_t *index3D_;

    size_t index(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l) const {
        return (*index3D_)[(i-1)*nn_+(j-1)*n_+k-1] + (l-k);
    }
};

class Matrix4DPF {
public:

    //! construct empty
    Matrix4DPF() {}

    void init(cand_pos_t n, index_offset_t &index3D) {
        n_=n;
        nn_ = n*n;
        index3D_ = &index3D;
        int slice_size_ = (*index3D_)[(n-1)*nn_ + (n-1)*n + n-1] + 1;
        assert( slice_size_ == n*(n+1)*(n+2)*(n+3)/24);
        m_.clear();
        m_.resize(slice_size_,0.0);

    }

    int get_uc(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l) const {
        assert(!(i<=0 || l> n_));

        int val = m_[index(i,j,k,l)];
        return val;
    }

    int get(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l) const {
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
        index.resize(n*n*n);
        cand_pos_t stride1 = n*n;
        size_t idx = 0;
        for (cand_pos_t i = 0; i < n; ++i) {
            for (cand_pos_t j = i; j < n; ++j) {
                for (cand_pos_t k = j; k < n; ++k) {
                    index[i*stride1+j*n+k] = idx;
                    idx += (n - k);  // remaining values from k to n
                }
            }
        }
    }

    private:
    cand_pos_t n_;
    cand_pos_t nn_;
    std::vector<pf_t> m_; // I had m as a 16 bit int. This seemed to result in some kind of of overflow or underflow, most likely because I allowed pushing INF(32 bit int/10) despite being 16 bit
    index_offset_t *index3D_;

    size_t index(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l) const {
        return (*index3D_)[(i-1)*nn_+(j-1)*n_+k-1] + (l-k);
    }
};

#endif