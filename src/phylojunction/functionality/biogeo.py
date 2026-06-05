import sys
import os
import enum
import math
from itertools import product
import numpy as np
import typing as ty

# pj imports
import phylojunction.utility.helper_functions as pjh
import phylojunction.functionality.feature_io as pjf

__author__ = "Fabio K. Mendes"
__email__ = "fmendes@lsu.edu"


class State2BitLookup:
    """
    Stash and machinery for checking and converting character
    compound-states into bit patterns, and vice-versa

    For example, if we are thinking about regions as characters,
    then compound-state means different ranges (e.g., "A", "AB"),
    while "state" is the number of different values a single
    character can take. In the case of 1-character-1-region, then
    the number of states is 2, because a species is either present
    or not in a region.

    Parameters:
        n_char (int): Number of places in bit pattern, which could be
            the number of atomic regions in a biogeographic system.
        n_states_per_char (int): Number of states a bit pattern can
            take. Normally this number is 2, since a bit can be either
            on (1) or off (0).
        n_states (int):
        int2bit_dict (dict):
        bit2int_dict (dict):
    """

    # e.g., for 2 states
    # Region | compound-state  | Bit pattern
    # --------------------------------------
    # Dead   |        0        | 000
    # A      |        1        | 100
    # B      |        2        | 010
    # C      |        3        | 001
    # AB     |        4        | 110
    # AC     |        5        | 101
    # BC     |        6        | 011
    # ABC    |        7        | 111

    n_char: int
    n_states_per_char: int
    n_states: int
    int2bit_dict: ty.Dict[str, str] = dict()
    bit2int_dict: ty.Dict[str, str] = dict()

    def __init__(self,
                 n_characters: int,
                 n_states_per_char: int,
                 geosse: bool = True) -> None:
        
        self.n_char = n_characters
        self.n_states_per_char = n_states_per_char
        self.n_states = int(n_states_per_char ** n_characters)

        # side-effect
        self._populate_dicts(geosse)

        if geosse:
            self.n_states -= 1  # remove null range


    # internal
    def _populate_dicts(self, geosse: bool = True) -> None:
        """
        Non-recursively generate all bit patterns for a given number
        of characters and states per characters

        Note that if we are thinking about regions as characters, then
        "compound state" means the range (e.g., "AB"), and state is the
        number of options per character (e.g., binary if present/absent
        in a region)
        """

        # non-recursive method
        list_for_sorting: ty.List[ty.List[str]] = \
            [[] for i in range(self.n_char + 1)]

        for compound_state_idx in range(self.n_states):
            bit_pattern = [0 for i in range(self.n_char)]
            z = compound_state_idx

            for bit_idx in range(self.n_char-1, -1, -1):
                # n_states_per_char here is present/absent for regions
                # (i.e., binary state)
                v = self.n_states_per_char ** bit_idx

                if z >= v:
                    bit_pattern[bit_idx] = 1
                    z -= v

                else:
                    bit_pattern[bit_idx] = 0

            bit_pattern_str = "".join(str(b) for b in bit_pattern)
            n_bits_on = sum(bit_pattern)

            # we do not want null range
            if geosse and n_bits_on == 0:
                continue

            list_for_sorting[n_bits_on].append(bit_pattern_str)

        # will store all bit patterns after sorting first by number of bits
        # that are on, and then by the decimal value underlying the bit pattern
        sorted_list_patterns: ty.List[str] = []
        for list_of_patterns in list_for_sorting:
            sorted_list_patterns.extend(list_of_patterns)

        # getting the other dict by reversing this one
        self.int2bit_dict = \
            dict((idx, bit_pattern)
                 for idx, bit_pattern
                 in enumerate(sorted_list_patterns))
        for k, v in self.int2bit_dict.items():
            self.bit2int_dict[v] = k

    # getters
    def get_bit(self, state: int, geosse: bool = True) -> str:
        return self.int2bit_dict[state]
    
    def get_bit_patts(self) -> ty.List[str]:
        return self.int2bit_dict.values()
    
    def get_state(self, bit_patt: str) -> int:
        return self.bit2int_dict[bit_patt]
    
    def get_widespread_bits(self) -> ty.Tuple[str]:
        widespread_bit_patt_list: ty.List[str] = list()
        
        for bit_patt in self.int2bit_dict.values():
            if bit_patt.count("1") > 1:
                widespread_bit_patt_list.append(bit_patt)

        return tuple(widespread_bit_patt_list)


class StatesGeoCondChangeLookup():

    n_regions: int  # i.e., number of characters
    state_bit_lookup_tup: ty.Tuple[State2BitLookup]
    geo_query: pjf.GeoFeatureQuery
    widespread_bit_patt_lookup_dict: \
        ty.Dict[str,
                ty.Dict[ty.Tuple[str],
                        ty.Dict[ty.Tuple[int],
                                ty.List[float]]]] \
                                = pjh.autovivify(3)
    # { widespread_bit1: {
    #                      (disjoint_region_sets_pair1): { (region_pair1): [ barrier_ages ], (region_pair2): [barrier_ages] },
    #                      (disjoint_region_sets_pair2): { (region_pair1): [ barrier_ages ], (region_pair2): [barrier_ages] },
    #                      (disjoint_region_sets_pair3): { (region_pair1): [ barrier_ages ], (region_pair2): [barrier_ages] }, ...
    #                     },
    #   widespread_bit2: { ... }
    #   ...
    # }

    def __init__(self,
                 n_regions: int,
                 geo_query: pjf.GeoFeatureQuery,
                 requirement_fn: pjf.MyCallableType,
                 requirement_fn_name: str):

        # this is like the number of characters
        self.n_regions = n_regions

        state_bit_lookup_list = list()
        for i in range(n_regions):
            state_bit_lookup_list.append(State2BitLookup(i+1, 2))

        self.state_bit_lookup_tup = tuple(state_bit_lookup_list)
        # self.state_bit_lookup = State2BitLookup(n_regions, 2)

        self.requirement_fn = requirement_fn
        
        self.geo_query = geo_query
        self.geo_query.populate_geo_cond_member_dicts(
            requirement_fn_name,
            requirement_fn)

        # debugging
        # for k in self.geo_query.geo_cond_change_times_dict["distance"]:
        #     print(*k)
        
        # debugging
        # print(self.state_bit_lookup.int2bit_dict)
        # print(self.state_bit_lookup.bit2int_dict)

        self._init_widespread_bit_patt_lookup_dict()

        # generate all bit patterns here
        #    self.generate_disjoint_bit_patt(bit_patt)

        self._debug_print_dict()
    
    def _init_widespread_bit_patt_lookup_dict(self):

        def _pairwise_one_pos_between_patts(one_pos_set1,
                                            negative_one_pos_set2):
            for i in product(one_pos_set1, negative_one_pos_set2):
                yield i
        
        state_bit_lookup = self.state_bit_lookup_tup[self.n_regions-1]

        # query from state2bit object for this number of 1's
        for widespread_bit_patt \
                in state_bit_lookup.get_widespread_bits():
            for disjoint_bit_patt_str_one_pos_tup, \
                negative_disjoint_bit_patt_str_one_pos_tup \
                    in self.disjoint_bit_patt_str_iter(widespread_bit_patt):
                
                disjoint_bit_patt_str, one_pos_set \
                    = disjoint_bit_patt_str_one_pos_tup
                
                # "negative" pattern is one in which the 0's are 1's and vice-versa
                neg_disjoint_bit_patt_str, neg_one_pos_set \
                    = negative_disjoint_bit_patt_str_one_pos_tup
                
                for pair_one_pos \
                        in _pairwise_one_pos_between_patts(one_pos_set, neg_one_pos_set):

                    from_pos, to_pos = pair_one_pos

                    # populating with geographic condition change times
                    self.widespread_bit_patt_lookup_dict[widespread_bit_patt] \
                        [disjoint_bit_patt_str + " " + neg_disjoint_bit_patt_str] \
                        [pair_one_pos] \
                            = self.geo_query.geo_cond_change_times_dict["distance"][from_pos][to_pos]
                    
                    # still populating! not assuming symmetry!
                    self.widespread_bit_patt_lookup_dict[widespread_bit_patt] \
                        [neg_disjoint_bit_patt_str + " " + disjoint_bit_patt_str] \
                        [(to_pos, from_pos)] \
                            = self.geo_query.geo_cond_change_times_dict["distance"][to_pos][from_pos]
                
        # print(self.widespread_bit_patt_lookup_dict)
                
    def _debug_print_dict(self):
        print("parent ch1 ch2 (reg pair partial barrier) [event times]")
        for widespread_bit_patt, disjoint_bit_patts_str_dict in self.widespread_bit_patt_lookup_dict.items():
            for disjoint_bit_patts_str, region_pair_tup_dict in disjoint_bit_patts_str_dict.items():
                for region_pair_tip, age_list in region_pair_tup_dict.items():
                        print(widespread_bit_patt,
                              disjoint_bit_patts_str,
                              region_pair_tip,
                              age_list)

    def disjoint_bit_patt_str_iter(self, widespread_bit_patt: str):
        idx_for_ones = set(idx for idx, el in enumerate(widespread_bit_patt) if el == "1")
        widespread_range_size = len(idx_for_ones)
        forbidden_patts = set(["0" * widespread_range_size,
                               "1" * widespread_range_size])

        state_bit_lookup = self.state_bit_lookup_tup[widespread_range_size-1]

        # the number of elements in bit_patt will be the
        # number of 1's in widespread_bit_patt
        for bit_patt_one_sized in state_bit_lookup.get_bit_patts():
            if bit_patt_one_sized not in forbidden_patts:
                bit_patt_one_sized_pos = set([idx for idx, one \
                                              in enumerate(bit_patt_one_sized) if one == "1"])

                lhs_bit_patt = list(widespread_bit_patt)
                rhs_bit_patt = list(widespread_bit_patt)

                negative_bit_patt = ["0" if i == "1" else "1" for i in list(bit_patt_one_sized)]
                negative_bit_patt_pos = idx_for_ones.difference(bit_patt_one_sized_pos)

                # debugging
                # print(bit_patt_one_sized_pos)
                # print(negative_bit_patt_pos)

                for i, bit_idx in enumerate(idx_for_ones):
                    lhs_bit_patt[bit_idx] = bit_patt_one_sized[i]
                    rhs_bit_patt[bit_idx] = negative_bit_patt[i]
                
                yield (("".join(lhs_bit_patt), bit_patt_one_sized_pos),
                       ("".join(rhs_bit_patt), negative_bit_patt_pos))

        # for bit_idx, el in enumerate(widespread_bit_patt):
        #     lhs_bit_patt = list(widespread_bit_patt)
        #     rhs_bit_patt = list(widespread_bit_patt)

        #     if bit_idx in idx_for_ones[:-1]:
        #         pos = idx_for_ones.index(bit_idx) + 1
        #         idx_for_ones_left_of_current = idx_for_ones[:pos]
        #         idx_for_ones_right_of_current = idx_for_ones[pos:]

        #         # debugging
        #         # print("pos = " + str(pos))
        #         # print("left non inclusive")
        #         # print(idx_for_ones_left_of_current)
        #         # print("right")
        #         # print(idx_for_ones_right_of_current)

        #         for i in idx_for_ones_left_of_current:
        #             lhs_bit_patt[i] = "0"
                
        #         for i in idx_for_ones_right_of_current:
        #             rhs_bit_patt[i] = "0"

        #         yield ("".join(lhs_bit_patt), "".join(rhs_bit_patt))

if __name__ == "__main__":
    
    fc = pjf.GeoFeatureCollection(
        sys.argv[1],
        age_summary_fp=sys.argv[2])
    
    fq = pjf.GeoFeatureQuery(fc)

    requirement_fn = \
        pjf.GeoFeatureQuery.qb_feature_threshold(fc, 15.0, True, feat_name="qb_1")

    # first argument here, number of regions, must <= the
    # number of features provided in feature_summary.csv
    biogeo_lookup = StatesGeoCondChangeLookup(4, fq, requirement_fn, "distance")

    # testing methods inside #
    test_bit_pattern = False
    if test_bit_pattern:

        biogeo_lookup2 = StatesGeoCondChangeLookup(2, fq, requirement_fn, "distance")
        print("\nex1: ", end="")
        for i, (lhs, rhs) in enumerate(biogeo_lookup2.disjoint_bit_patt_str_iter("11")):
            if i > 0:
                print("    ", lhs, " and ", rhs)

            else:
                print(lhs, " and ", rhs)
        # ex1: ('10', {0})  and  ('01', {1})
        #      ('01', {1})  and  ('10', {0})

        biogeo_lookup4 = StatesGeoCondChangeLookup(4, fq, requirement_fn, "distance")
        print("\nex2: ", end="")
        for i, (lhs, rhs) in enumerate(biogeo_lookup4.disjoint_bit_patt_str_iter("1001")):
            if i > 0:
                print("    ", lhs, " and ", rhs)

            else:
                print(lhs, " and ", rhs)
        # ex2: ('1000', {0})  and  ('0001', {3})
        #      ('0001', {1})  and  ('1000', {0, 3})        

        print("\nex3: ", end="")
        for i, (lhs, rhs) in enumerate(biogeo_lookup4.disjoint_bit_patt_str_iter("1011")):
            if i > 0:
                print("    ", lhs, " and ", rhs)

            else:
                print(lhs, " and ", rhs)
        # ex3: ('1000', {0})  and  ('0011', {2, 3})
        #      ('0010', {1})  and  ('1001', {0, 2, 3})
        #      ('0001', {2})  and  ('1010', {0, 3})
        #      ('1010', {0, 1})  and  ('0001', {2, 3})
        #      ('1001', {0, 2})  and  ('0010', {3})
        #      ('0011', {1, 2})  and  ('1000', {0, 3})
    
        print("\nex4: ", end="")
        for i, (lhs, rhs) in enumerate(biogeo_lookup4.disjoint_bit_patt_str_iter("1111")):
            if i > 0:
                print("    ", lhs, " and ", rhs)

            else:
                print(lhs, " and ", rhs)
        # ex4: ('1000', {0})  and  ('0111', {1, 2, 3})
        #      ('0100', {1})  and  ('1011', {0, 2, 3})
        #      ('0010', {2})  and  ('1101', {0, 1, 3})
        #      ('0001', {3})  and  ('1110', {0, 1, 2})
        #      ('1100', {0, 1})  and  ('0011', {2, 3})
        #      ('1010', {0, 2})  and  ('0101', {1, 3})
        #      ('0110', {1, 2})  and  ('1001', {0, 3})
        #      ('1001', {0, 3})  and  ('0110', {1, 2})
        #      ('0101', {1, 3})  and  ('1010', {0, 2})
        #      ('0011', {2, 3})  and  ('1100', {0, 1})
        #      ('1110', {0, 1, 2})  and  ('0001', {3})
        #      ('1101', {0, 1, 3})  and  ('0010', {2})
        #      ('1011', {0, 2, 3})  and  ('0100', {1})
        #      ('0111', {1, 2, 3})  and  ('1000', {0})