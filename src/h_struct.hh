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

// This node is used to keep the intervals that need to be further backtracked
// struct seq_interval
// {
//   int i;
//   int j;
//   int energy;                        // it is used
//   char type;
//   seq_interval* next;

//   void copy (seq_interval *other)
//   {
//     other->i = i;
//     other->j = j;
//     other->energy = energy;
//     other->type = type;
//   }
// };

typedef struct brack_type
	{
		char open;
		char close;
		void copy(brack_type *other){
			open = other->open;
			close = other->close;
		}
        brack_type(char _open, char _close):
        open(_open),
        close(_close)
        {
            
        }
	}brack_type;

typedef struct band_elem
{	
	band_elem *next;
	char open;
	char close;
	int outer_start;
	int outer_end;
	int inner_start;
	int inner_end;
	void copy(band_elem *other){
		other->outer_start = outer_start;
		other->outer_end = outer_end;
		other->inner_start = inner_start;
		other->inner_end = inner_end;
		other->open = open;
		other->close = close;
	}
    band_elem(char _open,char _close,int _outer_start, int _outer_end, int _inner_start, int _inner_end):
    open(_open),
    close(_close),
    outer_start(_outer_start),
    outer_end(_outer_end),
    inner_start(_inner_start),
    inner_end(_inner_end)
    {
        
    }
}band_elem;

struct seq_interval
{
  cand_pos_t i;
  cand_pos_t j;
  energy_t energy;                        // it is used
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
		//Hosna Feb 15, 2014
		// this part was added with similar reasoning as follows
		other ->k = k;
		other ->l = l;
		other -> asym = asym;
	}

};

struct free_energy_node
{
    int energy;
    char type;          // type may be: N (NONE), H (HAIRPIN), S (STACKED), I (INTERNAL), M (MULTI)
    free_energy_node()
    {
        energy = 10000; // INF
        type = NONE;
    }
};

#endif /*H_STRUCT_H_*/
