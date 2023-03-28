"""
This script generate insertion and deletion list according to given parameters
"""
import numpy as np


def generate_indel_lens(alpha, beta, num, len_max=10):
    """
    generate array of insertion/deletion lengths according to beta distribution
    :param alpha: alpha in beta distribution, default 0.5
    :param beta: beta in beta distribution, default 3
    :param num: total number of length to generate
    :param len_max: the max length of insertion/deletion  default 10
    :return: array of lengths
    """
    # alpha, beta, num, len_max = 0.5, 3, 20, 10
    indel_len_array = np.random.beta(alpha, beta, num) * len_max
    indel_len_array = indel_len_array.astype(np.int) + 1
    return indel_len_array


def generate_insertion_array(insertion_len_array, base_prop_dict):
    """
    :param insertion_len_array:
    :param base_prop_dict: like {"A": 1/4, "G": 1/4, "C": 1/4, "T": 1/4}
    :return: ['A', 'TCTAT', 'C', 'C',]
    the dict of base proportion, like geno_id: {"A": 1/4, "G": 1/4, "C": 1/4, "T": 1/4}
    """
    # insertion_len_array = np.array([1, 5, 1, 1, 8, 1, 2, 6, 4, 4, 1, 5, 1, 1, 2, 9, 5, 2, 6, 1])
    base_list = []
    prop_list = []
    for i, j in base_prop_dict.items():
        base_list.append(i)
        prop_list.append(j)
    # randomly select base from base_list according to prop_list
    insertion_array = ["".join(np.random.choice(base_list, i, replace=True, p=prop_list)) for i in insertion_len_array]
    return insertion_array


def random_long_indel_from_blocks(chr_len, avoided_region_list, long_indel_length_list):
    """
    assign long indels to the chromosome, not avoided_region_list
    # random_long_indel_from_blocks(10000, [[0, 1000], [5000,8000], [8500, 9500]], [500, 700, 100, 2000])
    :param chr_len: the length of current chromosome
    :param avoided_region_list:
    :param long_indel_length_list:
    :return: start position of long_indel, if -1000, this indel fail
    """
    remaining_blocks = [[0, chr_len]]               # list
    remaining_block_lengths = np.array([chr_len])   # array
    ############################################
    # remove nonmutated region from chromosome #
    ############################################
    if len(avoided_region_list) > 0:
        # check whether the avoided_region exceed chr_len
        avoided_region_passed = []
        for avoided_region in avoided_region_list:
            if 0 <= avoided_region[0] <= chr_len and 0 <= avoided_region[1] <= chr_len:
                avoided_region_passed.append(avoided_region)
            else:
                print(f"Error! The region {avoided_region} exceed the length of chromosome! Skip this region.")
        if len(avoided_region_passed) > 0:
            # ## need to consider avoid_region ###
            bin_left = [i[1] for i in avoided_region_passed]
            bin_right = [i[0] for i in avoided_region_passed]    # now len(bin_right) - bin_left = 1
            # check whether avoided_region start from 0
            if bin_right[0] <= 0:
                bin_right = bin_right[1:]       # right -1
            else:
                bin_left = [0] + bin_left       # left +1
            # check whether avoided_region end with chr_len
            if bin_left[-1] >= chr_len:
                bin_left = bin_left[:-1]        # left -1
            else:
                bin_right.append(chr_len)       # right +1
            remaining_blocks = [[i, j] for i, j in zip(bin_left, bin_right)]
            remaining_block_lengths = np.array([j-i for i, j in zip(bin_left, bin_right)])

    long_indel_pos_list = []        # where these indels should start from
    # insert long_indel_length_list
    for long_indel_length in long_indel_length_list:
        # check whether indel_length > remaining_block_lengths
        # indel can only be put into the block longer than the indel
        feasible_block_index = np.where(remaining_block_lengths > long_indel_length)[0]
        if len(feasible_block_index) == 0:
            # no feasible blocks for indel, assign a negative value to the start position
            long_indel_pos_list.append(-1000)
        else:
            # feasible_block_index = np.array([2, 5, 7]), feasible_block_prob = np.array([0.2, 0.2, 0.6])
            feasible_block_prob = remaining_block_lengths[feasible_block_index]/sum(remaining_block_lengths[feasible_block_index])
            chosen_block_index = np.random.choice(feasible_block_index, p=feasible_block_prob)
            chosen_block_range = remaining_blocks[chosen_block_index]       # [1000, 500000]
            # determine the start point of this indel, can not exceed chosen_block_range
            indel_start = np.random.choice(range(chosen_block_range[0], chosen_block_range[1]-long_indel_length))
            indel_end = indel_start + long_indel_length
            long_indel_pos_list.append(indel_start)
            # split the chosen block by indel_start and indel_end, update the remaining blocks
            remaining_blocks = remaining_blocks[:chosen_block_index] + [[chosen_block_range[0], indel_start], [indel_end, chosen_block_range[1]]] + remaining_blocks[chosen_block_index + 1:]
            remaining_block_lengths = np.array([i[1] - i[0] for i in remaining_blocks])
    return long_indel_pos_list
