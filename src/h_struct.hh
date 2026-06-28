#ifndef H_STRUCT_H_
#define H_STRUCT_H_

#include "constants.hh"
#include "base_types.hh"
#include <vector>
#include <array>
#include <string>

static constexpr std::array<std::pair<char,char>, 4> brackets = {{
    {'(', ')'},
    {'[', ']'},
    {'{', '}'},
    {'<', '>'},
}};

struct Band {
    cand_pos_t i, j;
    int  color;
	Band(cand_pos_t i, cand_pos_t j, int color): i(i), j(j), color(color){
	}
};

inline bool crosses(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l) {
    return (i < k && k < j && j < l) || (k < i && i < l && l < j);
}

// This is a graph coloring problem in essence
inline void fill_structure(std::vector<int> &fres,std::string &structure) {
	cand_pos_t n = structure.length();
    std::vector<Band> bands;
	for (cand_pos_t i = 0; i < n; ++i) {
        if (fres[i] != -2 && i < fres[i]){
			bands.emplace_back(i, fres[i], -1);
		}
    }
	for(cand_pos_t i = 0; i<(cand_pos_t) bands.size();++i){
		std::vector<int> cross = {0,0,0};
		for(cand_pos_t j=0; j<i;++j){
			if(crosses(bands[i].i,bands[i].j,bands[j].i,bands[j].j)){
				cross[bands[j].color] = 1; // Can just reassign 1 if it crosses. Better than an if else as branches can be expensive
			}
		}
		// The order here should ensure that it picks correctly. cross[color] ensures we stop at the first 0 in the cross array
		int color = 0;
        while (color < (int)cross.size() && cross[color]) ++color;
        bands[i].color = color;
	}
    for (cand_pos_t i = 0; i < (cand_pos_t) bands.size(); ++i) {
		auto [open, close] = brackets[bands[i].color];
		structure[bands[i].i] = open;
		structure[bands[i].j] = close;
    }
}

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
