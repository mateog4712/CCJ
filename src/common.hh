#ifndef COMMON_H
#define COMMON_H
#include "base_types.hh"

energy_t alpha1P(cand_pos_t z);					//penalty for z unpaired bases in an internal loop that spans a band
//int alpha2P(int i, int l);				// penalty for closing pair i.l of an internal loop that spans a band
energy_t alpha2P( cand_pos_t int_sequence_l, cand_pos_t int_sequence_i, cand_pos_t int_sequence_lminus1, cand_pos_t int_sequence_iplus1);
energy_t alpha3P(cand_pos_t z);	

// I feel like the above does not need to be their own functions. They are essentially just getting and interior loop or stack value for a pseudoknot. We have that with intP

energy_t beta2(cand_pos_t i, cand_pos_t l);		// penalty for closing pair i.l or l.i of an ordinary multiloop
energy_t beta2P(cand_pos_t i, cand_pos_t l);		// penalty for closing pair i.l or l.i of a multiloop that spans a band

energy_t gamma2(cand_pos_t i, cand_pos_t l);		// penalty for closing pair i.l or l.i of a pseudoloop

#endif