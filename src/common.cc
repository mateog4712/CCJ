#include "common.hh"
#include "h_externs.hh"


// penalty for closing pair i.l or l.i of an ordinary multiloop
energy_t beta2(cand_pos_t i, cand_pos_t l)
{
	// Hosna, April 2, 2014
	// I don't think this is the complete value, but since HFold's WM and CCJ's Vmloop recurrences are the same I am not changing this value here, unless I find out it is needed
	// the correct value should be: Non-GC-penalty(i,l)+b_penalty
	return b_penalty;
}

// penalty for closing pair i.l or l.i of a multiloop that spans a band
energy_t beta2P(cand_pos_t i, cand_pos_t l)
{
	// Hosna, April 2, 2014
	// I don't think this is the complete value, but since HFold's WM and CCJ's Vmloop recurrences are the same I am not changing this value here, unless I find out it is needed
	// the correct value should be: Non-GC-penalty(i,l) *0.74 + bp_penalty
	return bp_penalty;
}

// penalty for closing pair i.l or l.i of a pseudoloop
energy_t gamma2(cand_pos_t i, cand_pos_t l)
{
	// Hosna, April 2, 2014
	// I changed this value to be 0 as I can't find its correct value
	// the correct value for this penalty should be similar to what we have in case of an internal loop or a multiloop, but the value is missing here
	//return 0;
	// Hosna July 17, 2014
	// To avoid addition of single base pair bands I am giving a very small non-zero value to gamma2
	return 1;
}