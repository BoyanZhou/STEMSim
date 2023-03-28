"""
This script aims to modify sam files to generate evolution signals
module add bowtie2/2.3.5.1
module add samtools/1.9
"""

import numpy as np
import os
import pysam
import stemsim.distance_to_end

"""
@@@ testing example @@@
cd /gpfs/data/lilab/home/zhoub03/long_short_term_evolution/simulation_study/single_species/Bifidobacterium_breve/Bifidobacterium_breve_2strains_simu1/rep2/2022.11.01_16.59.20_sample_2/bam
unaligned_bam_abs_path = "/gpfs/data/lilab/home/zhoub03/long_short_term_evolution/simulation_study/single_species/Bifidobacterium_breve/Bifidobacterium_breve_2strains_simu1/rep2/2022.11.01_16.59.20_sample_2/bam/Genome1.0.bam"
bowtie_reference = "/gpfs/data/lilab/home/zhoub03/teddy/ncbi/dbGaP-21556/systematic_simulation/species_reference_genomes/Bifidobacterium_breve/Bifidobacterium_breve_12L/Bifidobacterium_breve"
thread = 4
bowtie2_path = "bowtie2"
samtools_path = "samtools"
process_simulated_data(unaligned_bam_abs_path, bowtie_reference, thread, bowtie2_path, samtools_path)

ref_fas_path = "/gpfs/data/lilab/home/zhoub03/teddy/ncbi/dbGaP-21556/systematic_simulation/species_reference_genomes/Bifidobacterium_breve/Bifidobacterium_breve_12L/Bifidobacterium_breve.fas"
pysam_ref = pysam.FastaFile(ref_fas_path)
aligned_bam_abs_path = "/gpfs/data/lilab/home/zhoub03/long_short_term_evolution/simulation_study/single_species/Bifidobacterium_breve/Bifidobacterium_breve_2strains_simu1/rep2/2022.11.01_16.59.20_sample_2/bam/Genome1.0_aligned_sorted.bam"
pysam_bam = pysam.AlignmentFile(aligned_bam_abs_path, "rb", reference_filename=ref_fas_path)
pysam_ref.references[0]
sample_pileup = pysam_bam.pileup(pysam_ref.references[0], 1000, 2000, truncate=True, stepper="samtools", fastafile=pysam_ref)
for pc in sample_pileup:
    print(pc.get_query_names())
"""


complementary_base_dict = {"A": "T",
                           "T": "A",
                           "U": "A",
                           "G": "C",
                           "C": "G"}


def reverse_complement(sequence):
    """
    Get the reverse complement of the input sequence
    :param sequence: sequence = "AGCTACGT"
    :return:
    """
    reverse_complement_seq = "".join([complementary_base_dict[i] for i in sequence[::-1]])
    return reverse_complement_seq


def get_nodes(chr_len, avoided_region_list, bin_len):
    """
    Get the inter_nodes of chromosome by setting length of bin; divide the whole genome by avoided_region
    For the divided genome segments, split each into bins if length > bin_len
    @ avoided_region_list = [[76,250], [700,900], [3100, 3333]]
    @ bin_len = 155
    @ chr_len = 3333
    :param chr_len:
    :param avoided_region_list: like [[0, 500], [2500, 3100]], can also be [] if no avoid region
    :param bin_len:
    :return:
    """
    if len(avoided_region_list) == 0:
        # ## no avoid_region provided ###
        bin_left = [i*bin_len for i in range(chr_len//bin_len + 1)]
        bin_right = [i * bin_len for i in range(1, chr_len // bin_len + 1)] + [chr_len]
        boundary = [[i, j] for i, j in zip(bin_left, bin_right)]
    else:
        # check whether the avoided_region exceed chr_len
        avoided_region_list2 = []
        for avoided_region in avoided_region_list:
            if 0 <= avoided_region[0] <= chr_len and 0 <= avoided_region[1] <= chr_len:
                avoided_region_list2.append(avoided_region)
            else:
                print(f"Error! The region {avoided_region} exceed the length of chromosome! Skip this region.")
        if len(avoided_region_list2) == 0:
            # ## no avoid_region provided ###
            bin_left = [i * bin_len for i in range(chr_len // bin_len + 1)]
            bin_right = [i * bin_len for i in range(1, chr_len // bin_len + 1)] + [chr_len]
            boundary = [[i, j] for i, j in zip(bin_left, bin_right)]
        else:
            # ## need to consider avoid_region ###
            bin_left = [i[1] for i in avoided_region_list2]
            bin_right = [i[0] for i in avoided_region_list2]    # now len(bin_right) - bin_left = 1
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
            """ split the divided segments longer than bin_len """
            boundary = []
            for region_left, region_right in zip(bin_left, bin_right):
                region_len = region_right - region_left
                if region_len >= bin_len*2:
                    bin_in_region_left = [(i * bin_len + region_left) for i in range(region_len // bin_len)]
                    bin_in_region_right = [(i * bin_len + region_left) for i in range(1, region_len // bin_len)] + \
                                          [region_right]
                    boundary.extend([[i, j] for i, j in zip(bin_in_region_left, bin_in_region_right)])
                else:
                    boundary.append([region_left, region_right])
    return boundary


# main function
def summarize_reads_to_modify(input_bam_list, ref_fas_path, ref_base_prop_vec, mutation_prop_combination_dict,
                              nonmutated_region_dict, long_insertion_dict, long_deletion_dict, insertion_dict,
                              deletion_length_dict, mutation_number_dict, base_mutation_freq_dict, depth_threshold=5,
                              bin_len=10000, space_between_mutations=150):
    """
    For each chromosome, then for each bin, summarize all info in all_mutation_info = {}, merge to all_mutation_info_chr
    at last merge to all_mutation_info.

    :param input_bam_list: list of abs path of bam
    :param ref_fas_path:
    :param ref_base_prop_vec: the proportion of four bases in mutations
    @ like [1/6, 2/6, 1/6, 2/6], in the order of "A", "C", "G", "T"

    :param mutation_prop_combination_dict: {"combination1": {"prob": 0.3, "longitudinal_prop": [0.1, 0.7, 0.8]}}
    @ the scenario of longitudinal proportion is sampled from mutation_prop_combination_dict
    @ there are many combination generated according to some distribution or given by user

    :param nonmutated_region_dict: {"chr1":[[200, 500], [700, 1000]]}, no overlap, from small to large
    :param long_insertion_dict: {"chr1":[[200, 500], [700, 1000]]}, no overlap, from small to large
    :param long_deletion_dict: {"chr1":[[200, 500], [700, 1000]]}, no overlap, from small to large
    :param insertion_dict: {"chr1":["AA", "G", "T", "AGT"]}
    :param deletion_length_dict: {"chr1":[1, 1, 1, 3, 5]}
    :param mutation_number_dict: total number of short-term mutations assigned to each chr, {"chr1": 20, "chr2":4}
    :param base_mutation_freq_dict: {"A": {"alt_base_list": ["G", "C", "T"], "alt_freq_list": [0.25, 0.5, 0.25]}, ...}
    :param depth_threshold: pos with depth > threshold will be used to generate mutations
    :param bin_len:
    :param space_between_mutations: distance between mutations (including small indels)
    :return: reads_mutations_record, reads_to_remove_record, all_mutation_info
    "type" is insertion or deletion or SNV
    {chr_name: {{position: {"type": "deletion", "base": [ref_base, mutated_base], "ref_count": [0] * len(input_bam_list), "alt_count": [0] * len(input_bam_list)}}}}
    """

    """ all results recorded """
    # record all created mutations, {"read_name1": {distance_to_end1: mutation_base}}, distance is 0-based
    reads_mutations_record = [{} for i in range(len(input_bam_list))]
    # reads that need to be removed, list of set
    reads_to_remove_record = [set() for i in range(len(input_bam_list))]
    # mutation info record, including position, base type (SNV, or insertion or deletion), mutation number
    all_mutation_info = {}  # {"chr1": {pos: {"type": "SNV", "base":[ref, alt], "ref_count": [], "alt_count": []}}}
    # long indel record
    long_indel_info = {}    # {"chr1": {pos: {"range": [], "type": "long_insertion", "prop": [0.1, 0.2, 0.5, 0.9]}}}

    pysam_ref = pysam.FastaFile(ref_fas_path)
    pysam_bam_list = [pysam.AlignmentFile(i, "rb", reference_filename=ref_fas_path) for i in input_bam_list]

    longitudinal_prop_combination_name_list = list(mutation_prop_combination_dict.keys())
    combination_prob_list = [mutation_prop_combination_dict[i]["prob"] for i in longitudinal_prop_combination_name_list]

    # ------------------------------------------------------------------------------------------------------------------
    """ assign mutations to bins for each chromosome """
    for chr_name, chr_len in zip(pysam_ref.references, pysam_ref.lengths):
        # summarize all information in all_mutation_info_chr for this chr_name
        all_mutation_info_chr = {}  # {{pos: {"type": "SNV", "base":[ref, alt], "ref_count": [], "alt_count": []}}}
        long_indel_info_chr = {}    # {"chr1": {pos: {"range": [,], "type": "long_insertion", "prop": [0.1, 0.2, 0.5, 0.9]}}}

        # get the nonmutated_region/mutations in this chromosome
        nonmutated_region_list = []
        if chr_name in nonmutated_region_dict:
            nonmutated_region_list = nonmutated_region_dict[chr_name]   # [[200, 500], [700, 1000]]

        long_insertion_list = []
        if chr_name in long_insertion_dict:
            long_insertion_list = long_insertion_dict[chr_name]

        long_deletion_list = []
        if chr_name in long_deletion_dict:
            long_deletion_list = long_deletion_dict[chr_name]

        insertion_list = []
        if chr_name in insertion_dict:
            insertion_list = insertion_dict[chr_name]               # list of insertion in this chr
        insertion_start_index = 0                               # insertion index in insertion_list, used in each bin

        deletion_length_list = []
        if chr_name in deletion_length_dict:
            deletion_length_list = deletion_length_dict[chr_name]   # list of len of each deletion
        deletion_start_index = 0                                # deletion index in insertion_list, used in each bin
        mutation_number = mutation_number_dict[chr_name]

        ##########################################
        # 1. check whether need to avoid regions #
        ##########################################
        # first need to combine all regions in nonmutated_region_list, long_insertion_list, long_deletion_list
        # these need to be avoided when generating mutations
        avoided_region_list = nonmutated_region_list + long_insertion_list + long_deletion_list
        avoided_region_sorted_list = sorted(avoided_region_list, key=lambda x: x[0])

        """ coverage summary by bins in primary chromosome """
        # add the function of avoiding regions to reads_count_summary()
        # reads_count_in_bins; row: bins, col: start, end, reads of all sample (n+2 cols)
        reads_count_in_bins = reads_count_summary(pysam_bam_list, chr_name, chr_len, bin_len,
                                                  avoided_region_sorted_list)

        #########################################
        # 2. create long insertion and deletion #
        #########################################
        if len(long_deletion_list) > 0:
            # long_deletion_list is the list of ranges of deletions, [[100, 1000], [2500, 3600]]
            # longitudinal_prop_combination_name_list = ["combination1", "combination2", "combination3"]
            # combination_prob_list = [0.2, 0.5, 0.3]; combination_assigned_to_mutation_list = ["combination1", ""]
            combination_assigned_to_long_deletion_list = np.random.choice(longitudinal_prop_combination_name_list,
                                                                     size=len(long_deletion_list),
                                                                     replace=True, p=combination_prob_list)
            prop_assigned_to_long_deletion_list = [mutation_prop_combination_dict[i]["longitudinal_prop"] for i in combination_assigned_to_long_deletion_list]
            # {pos: {"range": [,], "type": "long_insertion", "prop": [0.1, 0.2, 0.5, 0.9]}
            long_indel_info_chr.update({i[0]: {"range": i, "type": "long_deletion", "prop": j} for i, j in zip(long_deletion_list, prop_assigned_to_long_deletion_list)})
            # mutation_prop_of_pos_in_bin_array.shape is m * n, n: number of longitudinal sample, m: mutation number
            mutation_prop_of_long_deletion_array = np.array(prop_assigned_to_long_deletion_list)

            for bam_index, pysam_bam in enumerate(pysam_bam_list):
                # read ids that need to be removed for each pysam_bam
                read_ids_to_remove_set = reads_id_to_remove(pysam_bam, chr_name, long_deletion_list,
                                                            mutation_prop_of_long_deletion_array[:, bam_index])
                reads_to_remove_record[bam_index].update(read_ids_to_remove_set)

        """ insertion is like the complement of deletion """
        if len(long_insertion_list) > 0:
            combination_assigned_to_long_insertion_list = np.random.choice(longitudinal_prop_combination_name_list,
                                                                           size=len(long_insertion_list), replace=True,
                                                                           p=combination_prob_list)
            prop_assigned_to_long_insertion_list = [mutation_prop_combination_dict[i]["longitudinal_prop"] for i in combination_assigned_to_long_insertion_list]
            long_indel_info_chr.update({i[0]: {"range": i, "type": "long_insertion", "prop": j} for i, j in zip(long_insertion_list, prop_assigned_to_long_insertion_list)})
            mutation_prop_of_long_insertion_array = np.array(prop_assigned_to_long_insertion_list)
            """ need to calculate the complement, which is different from deletion """
            mutation_prop_of_long_insertion_complement_array = 1 - mutation_prop_of_long_insertion_array

            for bam_index, pysam_bam in enumerate(pysam_bam_list):
                # read ids that need to be removed for each pysam_bam
                read_ids_to_remove_set = reads_id_to_remove(pysam_bam, chr_name, long_insertion_list,
                                                            mutation_prop_of_long_insertion_complement_array[:,
                                                            bam_index])
                reads_to_remove_record[bam_index].update(read_ids_to_remove_set)

        #############################################
        # 3. create short mutations in the each chr #
        #############################################
        # assign SNVs and InDels to each bin according to the length of each bin, multinomial distribution
        bin_len_array = reads_count_in_bins[:, 1] - reads_count_in_bins[:, 0]   # length of each bin
        bin_prob_list = list(bin_len_array/np.sum(bin_len_array))
        # mutations_each_bin_array = array([0, 1, 0, 0, 0, 0, 0, 1, 0, 0]), including snv and indel
        snv_each_bin_array = np.random.multinomial(n=mutation_number, pvals=bin_prob_list)
        insertion_each_bin_array = np.random.multinomial(n=len(insertion_list), pvals=bin_prob_list)
        deletion_each_bin_array = np.random.multinomial(n=len(deletion_length_list), pvals=bin_prob_list)
        mutations_each_bin_array = snv_each_bin_array + insertion_each_bin_array + deletion_each_bin_array
        # bin_boundary_list = get_nodes(pysam_ref.lengths[0], bin_len)     # [[0, 10000], [10000, 20000]]
        index_of_bin_with_mutations_array = np.where(mutations_each_bin_array > 0)[0]
        # ------------------------------------------------------------------------------------------------------------------
        """ Assign mutations in each bin """
        for bin_index in index_of_bin_with_mutations_array:
            """ Get the coverage at each pos in the bin for each BAM """
            bin_boundary = reads_count_in_bins[bin_index, :2]       # [0, 10000]
            # depth of each position in this bin
            bams_bin_depth_array = get_depth_in_bin(chr_name, bin_boundary[0], bin_boundary[1], pysam_bam_list)
            if len(bams_bin_depth_array) == 0:
                # means no coverage in this block
                continue
            # sum of depth in all bams
            bams_bin_depth_sum_array = np.sum(bams_bin_depth_array, axis=1)         # len of end-start

            # ----------------------------------------------------------------------------------------------------------
            """ Get passed pos and corresponding base """
            pos_index_passed_bool = bams_bin_depth_sum_array >= depth_threshold     # whether pass depth threshold
            if sum(pos_index_passed_bool) == 0:
                # no pos pass threshold
                continue
            pos_index_passed_blocks = get_true_blocks(pos_index_passed_bool)        # [[0, 30], [50, 70]]
            # remove two ends (default length is ) from blocks
            pos_index_passed_blocks_filtered = filter_blocks_by_boundary(pos_index_passed_blocks)

            if len(pos_index_passed_blocks_filtered) == 0:
                # all blocks too short
                continue

            # pos_index_passed = np.where(pos_index_passed_bool)[0]                           # index in the bin
            fasta_bin = pysam_ref.fetch(reference=chr_name, start=bin_boundary[0], end=bin_boundary[1], region=None)
            ref_bases_bin_array = np.array(list(fasta_bin.upper()))     # array(['A', "G", "C", "T"])
            # ref_bases_passed_array = ref_bases_bin_array[pos_index_passed]

            # ----------------------------------------------------------------------------------------------------------
            """ Assign mutations of each type to passed positions """
            # ## get the pos need to be modified and corresponding original genotype ###
            # ref_base_prop_vec = [1/6, 2/6, 1/6, 2/6], mutation_number=60  in the order of "A", "C", "G", "T"
            # get a vector of number of 4 bases, array([10, 19,  9, 22])
            ##############
            # Attention! #
            ##############
            # InDels should not be too close to the end of reads or border of pos_index_passed
            all_mutation_info_in_bin = {}  # {pos: {"type": "SNV", "base":[ref, alt], "ref_count": [[0] * len(input_bam_list)], "alt_count": [[0] * len(input_bam_list)]}}
            """ 3.1 If deletions exist """
            # ref: GTCT, alt: G;
            deletion_num_this_bin = deletion_each_bin_array[bin_index]
            if deletion_num_this_bin > 0:
                # update pos_index_passed_blocks_filtered, because selected pos and surrounding regions removed
                # e.g. [0, 500], remove 200 and surrounding, become [[0, 50], [350, 500]]
                deletion_positions, pos_index_passed_blocks_filtered = random_positions_from_blocks(
                    pos_index_passed_blocks_filtered, deletion_num_this_bin, space_between_mutations)
                for deletion_position in deletion_positions:
                    deletion_length_current = deletion_length_list[deletion_start_index]    # the len of deletion
                    ref_base = "".join(
                        ref_bases_bin_array[deletion_position: (deletion_position + deletion_length_current + 1)])      # "".join(np.array(['A', "G"]))
                    mutated_base = ref_bases_bin_array[deletion_position]
                    all_mutation_info_in_bin.update({deletion_position + + bin_boundary[0]: {"type": "deletion",
                                                                         "base": [ref_base, mutated_base],
                                                                         "ref_count": [0] * len(input_bam_list),
                                                                         "alt_count": [0] * len(input_bam_list)}})
                    deletion_start_index += 1
            # #----------------------------------------------------------------
            """ 3.2 If insertions exist and still region available """
            # ref: G, alt: GTCT;
            insertion_num_this_bin = insertion_each_bin_array[bin_index]
            if len(pos_index_passed_blocks_filtered) > 0 and insertion_num_this_bin > 0:
                insertion_positions, pos_index_passed_blocks_filtered = random_positions_from_blocks(
                    pos_index_passed_blocks_filtered, insertion_num_this_bin, space_between_mutations)
                for insertion_position in insertion_positions:
                    insertion_current = insertion_list[insertion_start_index]  # real insertion here, "AGT"
                    ref_base = ref_bases_bin_array[insertion_position]
                    mutated_base = ref_base + insertion_current
                    all_mutation_info_in_bin.update({insertion_position + bin_boundary[0]: {"type": "insertion",
                                                                         "base": [ref_base, mutated_base],
                                                                         "ref_count": [0] * len(input_bam_list),
                                                                         "alt_count": [0] * len(input_bam_list)}})
                    insertion_start_index += 1
            # #----------------------------------------------------------------
            """ 3.3 If snv exist and still region available """
            snv_chosen_index_base_dict = {}  # {1: "A", 17: "G"}
            snv_num_this_bin = snv_each_bin_array[bin_index]
            if len(pos_index_passed_blocks_filtered) > 0 and snv_num_this_bin > 0:
                # pos_index_passed_blocks_filtered = [[10, 25], [40, 60]]
                # to the form like pos_index_passed = array([10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23,
                # 24, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59])
                pos_index_passed = np.concatenate([np.arange(start, end) for start, end in pos_index_passed_blocks_filtered])
                ref_bases_passed_array = ref_bases_bin_array[pos_index_passed]          # array(["A", "G", ... ...])
                base_count_AGCT = np.random.multinomial(n=snv_num_this_bin, pvals=ref_base_prop_vec)

                for base_i, base_i_count in zip(["A", "C", "G", "T"], base_count_AGCT):
                    """ randomly assign base i to passed pos """
                    # base_i = "A"
                    base_i_index = pos_index_passed[ref_bases_passed_array == base_i]   # index of pos in the bin
                    if len(base_i_index) > base_i_count:
                        # more positions for mutations
                        snv_chosen_index_base_dict.update(
                            {i: base_i for i in list(np.random.choice(base_i_index, size=base_i_count, replace=False))})
                    else:
                        # candidate index is not enough
                        snv_chosen_index_base_dict.update({i: base_i for i in list(base_i_index)})

                """
                Get base that it mutates to 
                """
                for pos_i, ref_base in snv_chosen_index_base_dict.items():
                    alt_base_list = base_mutation_freq_dict[ref_base]["alt_base_list"]
                    alt_freq_list = base_mutation_freq_dict[ref_base]["alt_freq_list"]
                    mutated_base = np.random.choice(alt_base_list, 1, alt_freq_list)[0]  # randomly mutated to another base
                    all_mutation_info_in_bin.update({pos_i + bin_boundary[0]: {"type": "SNV",
                                                                               "base": [ref_base, mutated_base],
                                                                               "ref_count": [0] * len(input_bam_list),
                                                                               "alt_count": [0] * len(input_bam_list)}})
            # #------------------------------------------------------------------------------------------------------
            """ 3.4 Generate mutation_freq at each pos from mutation_prop_combination_dict """
            # longitudinal_prop_combination_name_list = ["combination1", "combination2", "combination3"]
            # combination_prob_list = [0.2, 0.5, 0.3]; combination_assigned_to_mutation_list = ["combination1", ""]
            combination_assigned_to_mutation_list = np.random.choice(longitudinal_prop_combination_name_list,
                                                                         size=len(all_mutation_info_in_bin),
                                                                         replace=True, p=combination_prob_list)

            # mutation_prop_of_pos_in_bin_array.shape is m * n, n: number of longitudinal sample, m: mutation number
            mutation_prop_of_pos_in_bin_array = np.array(
                    [mutation_prop_combination_dict[i]["longitudinal_prop"] for i in
                     combination_assigned_to_mutation_list])

            """ Record reads info to change for each bam """
            pos_index_in_chr_sorted = sorted(all_mutation_info_in_bin.keys())
            for bam_index, pysam_bam in enumerate(pysam_bam_list):
                # in current bin, for each longitudinal sample
                chosen_index_mutation_freq_dict = {chosen_index: mutation_freq for chosen_index, mutation_freq in
                                                   zip(pos_index_in_chr_sorted, mutation_prop_of_pos_in_bin_array[:, bam_index])}
                # all_mutation_info_in_bin = {}  # {pos: {"base":[ref, alt], "ref_count": [], "alt_count": []}}
                read_info_to_change_in_a_bam_dict = get_reads_info_to_modify(pysam_bam, pysam_ref,
                                                                             chr_name, bin_boundary[0], bin_boundary[1],
                                                                             chosen_index_mutation_freq_dict,
                                                                             all_mutation_info_in_bin,
                                                                             bam_index)
                # update mutations info on reads to total record list (no need to discriminate chr)
                reads_mutations_record[bam_index].update(read_info_to_change_in_a_bam_dict)
            # record the info of generated mutations in all_mutation_info_chr
            all_mutation_info_chr.update(all_mutation_info_in_bin)
        # for each chromosome, update all_mutation_info_chr to all_mutation_info
        all_mutation_info.update({chr_name: all_mutation_info_chr})  # pos and reads count
        long_indel_info.update({chr_name: long_indel_info_chr})
    return reads_mutations_record, reads_to_remove_record, all_mutation_info, long_indel_info


def get_depth_in_bin(chr_name, start, end, pysam_bam_list):
    """
    Get the depth of all bams in the region of chr_name
    :param chr_name:
    :param start:
    :param end:
    :param pysam_bam_list:
    :return: depth array, (end-start, pysam_bam_number)
    """
    bams_bin_depth = []
    for pysam_bam in pysam_bam_list:
        # bin_depth.shape=(4, end-start)
        bin_depth_array = np.array(pysam_bam.count_coverage(contig=chr_name, start=start, stop=end))
        bams_bin_depth.append(list(np.sum(np.array(bin_depth_array), axis=0)))
    return np.array(bams_bin_depth).T


def reads_count_summary(pysam_bam_list, chr_name, chr_len, bin_len, avoided_region_sorted_dict):
    """
    Summary reads count in each bin of each chr
    :param pysam_bam_list:
    :param chr_name:
    :param chr_len:
    :param bin_len:
    :param avoided_region_sorted_dict:{"chr1": [[0, 500], [2500, 3100]]}, can be {}
    mutations will not be generated in these regions
    :return: {chr1: array([1000, 2000, 50,])}, reads count of each bin
    """
    # reads_in_bin_by_chromosome = {}
    # for this chromosome of pysam_ref
    reads_bin_chromosome = []       # row: bins, col: start, end, reads of all sample (n+2 cols)
    # for each bin, 10000bp example; consider avoid region
    avoided_region_chr = []
    if chr_name in avoided_region_sorted_dict:
        avoided_region_chr = avoided_region_sorted_dict[chr_name]     # like [[0, 500], [2500, 3100]]
    for boundary in get_nodes(chr_len, avoided_region_chr, bin_len):
        # get reads of all samples at this bin
        boundary.extend([pysam_sample.count(contig=chr_name, start=boundary[0], stop=boundary[1]) for pysam_sample in
                         pysam_bam_list])
        reads_bin_chromosome.append(boundary)
    # update reads in bins
    # reads_in_bin_by_chromosome.update({chr_name: np.array(reads_bin_chromosome)})
    return np.array(reads_bin_chromosome)


def get_true_blocks(bool_arr):
    """
    Get the block_boundary where is True
    :param bool_arr: like [True, True, True, False, False, True, True]
    :return: [[0: 3], [5: 7]
    """
    if not isinstance(bool_arr, np.ndarray) or bool_arr.dtype != np.bool:
        raise ValueError("Input must be a NumPy array of booleans.")
    true_blocks = []
    start, end = None, None
    for i, val in enumerate(bool_arr):
        if val and start is None:
            start = i
        elif not val and start is not None:
            end = i
            true_blocks.append([start, end])
            start, end = None, None
    if start is not None:
        end = len(bool_arr)
        true_blocks.append([start, end])
    return true_blocks


def filter_blocks_by_boundary(input_blocks, min_distance_to_boundary=5):
    """
    :param input_blocks: [[0, 30], [50, 70], [90, 93]]
    :param min_distance_to_boundary: 5bp to boundar was excluded
    :return: [[5, 25], [55, 65]]
    """
    filter_blocks = []
    for block_boundary in input_blocks:
        block_boundary[0] += min_distance_to_boundary
        block_boundary[1] -= min_distance_to_boundary
        if block_boundary[1] > block_boundary[0]:
            # only report block not too close to boundary
            filter_blocks.append(block_boundary)
    return filter_blocks


def get_reads_info_to_modify(pysam_bam, pysam_ref, chr_name, start, end,
                             chosen_index_mutation_freq_dict, mutation_info_in_bin, index_in_bam_list):
    """
    Given the index (in bin) of pos to be modified, summarize the read names and position in reads need to be changed
    from the pysam_bam input
    :param pysam_bam:
    :param pysam_ref:
    :param chr_name:
    :param start: start pos of the bin
    :param end:
    # :param chosen_index_base_dict: {1: "A", 17: "G"}, the index within the bin and corresponding base to change
    :param chosen_index_mutation_freq_dict: {1: 0.1, 17: 0.7"}, key is the pos in chr, not in bin
    :param mutation_info_in_bin: {pos: {"type": "SNV", "base":[ref, alt], "ref_count": [0, 0, 0],
    "alt_count": [0, 0, 0]}}, key is the pos in genome
    :param index_in_bam_list: 0, 1, ...
    :return: {"read_name1": {"range_list":[[distance_to_end_start, distance_to_end_end]], "alt_list": [mutation_base]}},
    distance is 0-based change the base in region [distance_to_end_start, distance_to_end_end] of read to mutation_base
    there may be multiple range if multiple mutations
    For deletion, e.g. we change "ACGT" to "A"
    """
    read_info_to_change_dict = {}   # {"read_name1": {distance_to_end1: mutation_base}}, distance is 0-based

    bin_pileup = pysam_bam.pileup(chr_name, start, end, truncate=True, stepper="samtools", fastafile=pysam_ref)
    for pc in bin_pileup:
        # pc_index = pc.reference_pos - start
        if pc.reference_pos in mutation_info_in_bin:
            """ 
            If it is the position need to be mutated 
            """
            # {"type": "SNV", "base":[ref, alt], "ref_count": [0, 0, 0], "alt_count": [0, 0, 0]}
            mutation_info_this_pos = mutation_info_in_bin[pc.reference_pos]
            # two arrays
            pos_to_left_end, whether_reverse = stemsim.distance_to_end.summarize_distance(pc)
            read_name_array = np.array(pc.get_query_names())
            """ 
            check whether reads number is the same with pos_to_left_end, whether_reverse 
            """
            if len(read_name_array) != len(pos_to_left_end) or len(read_name_array) != len(whether_reverse):
                print(f"A warning: The mutation at {pc.reference_pos + 1} is not generated.")
                continue    # skip this pos

            """ record the number of ref and alt """
            # number_of_mutated_reads = reads_num * mutation_rate
            print(f"Start is {start}, end is {end}")
            print(f"chosen_index_mutation_freq_dict is {chosen_index_mutation_freq_dict}")
            print(f"mutation_info_in_bin is {mutation_info_in_bin}")
            number_of_mutated_reads = round(chosen_index_mutation_freq_dict[pc.reference_pos] * len(read_name_array))
            number_of_not_mutated_reads = len(read_name_array) - number_of_mutated_reads
            mutation_info_this_pos["ref_count"][index_in_bam_list] += number_of_not_mutated_reads
            mutation_info_this_pos["alt_count"][index_in_bam_list] += number_of_mutated_reads

            # the index of chosen read in the pc
            chosen_reads_index = np.random.choice(np.array(range(len(read_name_array))), number_of_mutated_reads,
                                                  replace=False)

            ref_base = mutation_info_this_pos["base"][0]        # the ref base
            mutated_base = mutation_info_this_pos["base"][1]    # the alt base
            for chosen_index in chosen_reads_index:
                read_name = read_name_array[chosen_index]
                if mutation_info_this_pos["type"] == "SNV":
                    # ref is like "A", alt is like "T"
                    if whether_reverse[chosen_index]:
                        """ read mapped to the reverse strand """
                        if read_name in read_info_to_change_dict:
                            read_info_to_change_dict[read_name]["range_list"].append(
                                [pos_to_left_end[chosen_index], pos_to_left_end[chosen_index] + 1])
                            read_info_to_change_dict[read_name]["alt_list"].append(
                                complementary_base_dict[mutated_base])
                        else:
                            read_info_to_change_dict.update({read_name: {
                                "range_list": [[pos_to_left_end[chosen_index], pos_to_left_end[chosen_index] + 1]],
                                "alt_list": [complementary_base_dict[mutated_base]]}})
                    else:
                        """ read mapped to the forward strand """
                        if read_name in read_info_to_change_dict:
                            read_info_to_change_dict[read_name]["range_list"].append([pos_to_left_end[chosen_index],
                                                                                      pos_to_left_end[chosen_index] + 1])
                            read_info_to_change_dict[read_name]["alt_list"].append(mutated_base)
                        else:
                            read_info_to_change_dict.update({read_name: {"range_list": [[pos_to_left_end[chosen_index],
                                                                                         pos_to_left_end[chosen_index] + 1]],
                                                                         "alt_list": [mutated_base]}})

                elif mutation_info_this_pos["type"] == "insertion":
                    # ref is like "A", alt is like "AGTT"
                    if whether_reverse[chosen_index]:
                        """ read mapped to the reverse strand """
                        if read_name in read_info_to_change_dict:
                            read_info_to_change_dict[read_name]["range_list"].append(
                                [pos_to_left_end[chosen_index], pos_to_left_end[chosen_index] + 1])
                            read_info_to_change_dict[read_name]["alt_list"].append(
                                reverse_complement(mutated_base))
                        else:
                            read_info_to_change_dict.update({read_name: {
                                "range_list": [[pos_to_left_end[chosen_index], pos_to_left_end[chosen_index] + 1]],
                                "alt_list": [reverse_complement(mutated_base)]}})
                    else:
                        """ read mapped to the forward strand """
                        if read_name in read_info_to_change_dict:
                            read_info_to_change_dict[read_name]["range_list"].append([pos_to_left_end[chosen_index],
                                                                                      pos_to_left_end[chosen_index] + 1])
                            read_info_to_change_dict[read_name]["alt_list"].append(mutated_base)
                        else:
                            read_info_to_change_dict.update({read_name: {"range_list": [[pos_to_left_end[chosen_index],
                                                                                         pos_to_left_end[chosen_index] + 1]],
                                                                         "alt_list": [mutated_base]}})
                elif mutation_info_this_pos["type"] == "deletion":
                    # ref is like "AGTT", alt is like "A"
                    # !!! be careful the deletion can't be out of the range of reads
                    if whether_reverse[chosen_index]:
                        """ read mapped to the reverse strand """
                        if read_name in read_info_to_change_dict:
                            read_info_to_change_dict[read_name]["range_list"].append(
                                [pos_to_left_end[chosen_index] - len(ref_base) + 1,
                                 pos_to_left_end[chosen_index] + 1])
                            read_info_to_change_dict[read_name]["alt_list"].append(
                                complementary_base_dict[mutated_base])
                        else:
                            read_info_to_change_dict.update({read_name: {
                                "range_list": [[pos_to_left_end[chosen_index] - len(ref_base) + 1,
                                                pos_to_left_end[chosen_index] + 1]],
                                "alt_list": [complementary_base_dict[mutated_base]]}})
                    else:
                        """ read mapped to the forward strand """
                        if read_name in read_info_to_change_dict:
                            read_info_to_change_dict[read_name]["range_list"].append([pos_to_left_end[chosen_index],
                                                                                      pos_to_left_end[chosen_index] +
                                                                                      len(ref_base)])
                            read_info_to_change_dict[read_name]["alt_list"].append(mutated_base)
                        else:
                            read_info_to_change_dict.update({read_name: {"range_list": [
                                [pos_to_left_end[chosen_index], pos_to_left_end[chosen_index] + len(ref_base)]],
                                                                         "alt_list": [mutated_base]}})
    return read_info_to_change_dict


def random_positions_from_blocks(blocks, num_positions, min_distance=150):
    """
    Select position from blocks
    blocks = [[0, 100], [150, 400], [500, 700]]
    num_positions = 5
    min_distance = 50
    :param blocks:
    :param num_positions:
    :param min_distance: minimal distance between selected positions, and remove region surrounding positions
    :return: positions, remaining_blocks
    """
    positions = []
    remaining_blocks = blocks.copy()                            # [[0, 100], [150, 400], [500, 700]]
    for i in range(num_positions):
        # Check if there are any remaining blocks left
        if len(remaining_blocks) == 0:
            break
        # Choose a random block from the list
        remaining_blocks_array = np.array(remaining_blocks)     # np.array([[0, 100], [150, 400], [500, 700]])
        remaining_blocks_len_array = remaining_blocks_array[:, 1] - remaining_blocks_array[:, 0]

        # choose block according to the length of each block
        chosen_block_index = np.random.choice(range(remaining_blocks_array.shape[0]),
                                              p=remaining_blocks_len_array/sum(remaining_blocks_len_array))

        # Generate a random position within the block
        chosen_block = remaining_blocks[chosen_block_index]     # [500, 700]
        position = np.random.randint(chosen_block[0], chosen_block[1])
        positions.append(position)

        # delete the chosen block
        del remaining_blocks[chosen_block_index]

        # check whether the chosen block ([500, 700]) remain after pick this position (589)
        if (position - min_distance) > chosen_block[0]:
            remaining_blocks.append([chosen_block[0], position - min_distance])
        if (position + min_distance) < chosen_block[1]:
            remaining_blocks.append([position + min_distance, chosen_block[1]])

    return positions, remaining_blocks


def reads_id_to_remove(pysam_bam, chr_name, deletion_ranges_list, mutation_prop_of_long_deletion_array):
    """
    Get the ids of reads that need to be removed in all long deletions
    :param pysam_bam:
    :param chr_name:
    :param deletion_ranges_list:
    :param mutation_prop_of_long_deletion_array: prop of each deletion, one-dim array
    :return: set of reads id need to remove
    """
    remove_reads_set = set()  # remove_reads_set for all deletion_ranges

    for deletion_range, mutation_prop in zip(deletion_ranges_list, mutation_prop_of_long_deletion_array):
        # for each deletion region, get all reads id in that region; mutation_prop is value (0-1) indicating the prop
        reads_id_list = [i.query_name for i in pysam_bam.fetch(chr_name, deletion_range[0], deletion_range[1])]
        reads_to_remove = np.random.choice(reads_id_list, size=int(len(reads_id_list)*mutation_prop), replace=False)
        remove_reads_set.update(set(reads_to_remove))

    return remove_reads_set



