#ifndef H_STRUCT_H_
#define H_STRUCT_H_

#include "constants.hh"
#include "base_types.hh"
#include <vector>

// the data structure stored in the V array
typedef struct minimum_fold
{
    short int pair;
    char type;                  // type can be 'H', 'S', 'I', 'M'
    minimum_fold()
    {
        pair = -1;
        type = NONE;
    }
} minimum_fold;

struct seq_interval
{
  cand_pos_t i;
  cand_pos_t j;
  energy_t energy;
  char type;
  seq_interval* next = nullptr;
	// Hosna, Feb 15, 2014
	// adding the following so that I can use stack_interval for backtracking the gapped region in CCJ also.
	// I am defining the gapped region as [i,k]U[l,j] instead of [i,j]U[k,l] in the recurrences for compatibility with simfold's seq_interval
	cand_pos_t k;
	cand_pos_t l;
	cand_pos_t asym;

	void copy (seq_interval *other)
	{
		other->i = i;
		other->j = j;
		other->energy = energy;
		other->type = type;
		other ->k = k;
		other ->l = l;
		other -> asym = asym;
	}

};

struct free_energy_node
{
    energy_t energy;
    char type;          // type may be: N (NONE), H (HAIRPIN), S (STACKED), I (INTERNAL), M (MULTI)
    free_energy_node()
    {
        energy = 10000; // INF
        type = NONE;
    }
};

#endif /*H_STRUCT_H_*/
