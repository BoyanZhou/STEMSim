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


def get_nodes(chr_len, bin_len):
    """ get the inter_nodes of chromosome by setting length of bin """
    bin_left = [i*bin_len for i in range(chr_len//bin_len + 1)]
    bin_right = [i * bin_len for i in range(1, chr_len // bin_len + 1)] + [chr_len]
    boundary = [[i, j] for i, j in zip(bin_left, bin_right)]
    return boundary


# main function
def summarize_reads_to_modify(input_bam_list, ref_fas_path, ref_base_prop_vec, mutation_prop_combination_dict,
                              mutation_number, base_mutation_freq_dict, depth_threshold=5, bin_len=10000):
    """

    :param input_bam_list: list of abs path of bam
    :param ref_fas_path:
    :param ref_base_prop_vec: the proportion of four bases in mutations
    @ like [1/6, 2/6, 1/6, 2/6], in the order of "A", "G", "C", "T"

    :param mutation_prop_combination_dict: {"combination1": {"prob": 0.3, "longitudinal_prop": [0.1, 0.7, 0.8]}}
    @ the scenario of longitudinal proportion is sampled from mutation_prop_combination_dict
    @ there are many combination generated according to some distribution or given by user

    :param mutation_number: total number of short-term mutations
    :param base_mutation_freq_dict: {"A": {"alt_base_list": ["G", "C", "T"], "alt_freq_list": [0.25, 0.5, 0.25]}, ...}
    :param depth_threshold:
    :param bin_len:
    :return:
    """

    """ all results recorded """
    # record all created mutations, {"read_name1": {distance_to_end1: mutation_base}}, distance is 0-based
    reads_mutations_record = [{} for i in range(len(input_bam_list))]
    # mutation info record, including position, base type, mutation number
    all_mutation_info = {}  # {pos: {"base":[ref, alt], "ref_count": [], "alt_count": []}}

    pysam_ref = pysam.FastaFile(ref_fas_path)
    pysam_bam_list = [pysam.AlignmentFile(i, "rb", reference_filename=ref_fas_path) for i in input_bam_list]

    longitudinal_prop_combination_name_list = list(mutation_prop_combination_dict.keys())
    combination_prob_list = [mutation_prop_combination_dict[i]["prob"] for i in longitudinal_prop_combination_name_list]

    ############################
    # assign mutations to bins #
    ############################
    """ coverage summary by bins in primary chromosome """
    reads_count_each_chr_dict = reads_count_summary(pysam_bam_list, pysam_ref, bin_len)
    # create mutations in the primary chr: pysam_ref.references[0]
    reads_count_in_bins = reads_count_each_chr_dict[pysam_ref.references[0]]
    # row: bins, col: start, end, reads of all sample (n+2 cols)
    bin_num = reads_count_in_bins.shape[0]
    # mutations_each_bin_array = array([0, 1, 0, 0, 0, 0, 0, 1, 0, 0])
    mutations_each_bin_array = np.random.multinomial(n=mutation_number, pvals=[1/bin_num] * bin_num)
    bin_boundary_list = get_nodes(pysam_ref.lengths[0], bin_len)     # [[0, 10000], [10000, 20000]]
    index_of_bin_with_mutations_array = np.where(mutations_each_bin_array > 0)[0]
    # ------------------------------------------------------------------------------------------------------------------
    """ Assign mutations in each bin """
    for bin_index in index_of_bin_with_mutations_array:

        """ Get the coverage at each pos in the bin for each BAM """
        bin_boundary = bin_boundary_list[bin_index]     # [0, 10000]
        bams_bin_depth_array = get_depth_in_bin(pysam_ref.references[0], bin_boundary[0], bin_boundary[1], pysam_bam_list)
        bams_bin_depth_sum_array = np.sum(bams_bin_depth_array, axis=1)     # len of end-start

        # --------------------------------------------------------------------------------------------------------------

        """ Get passed pos and corresponding base """
        pos_index_passed = np.where(bams_bin_depth_sum_array >= depth_threshold)[0] # index in the bin
        fasta_bin = pysam_ref.fetch(reference=pysam_ref.references[0], start=bin_boundary[0], end=bin_boundary[1], region=None)
        ref_bases_bin_array = np.array(list(fasta_bin.upper()))     # array(['A', "G", "C", "T"])
        ref_bases_passed_array = ref_bases_bin_array[pos_index_passed]

        # --------------------------------------------------------------------------------------------------------------

        """ Assign mutations of each type to passed positions """
        # ## get the pos need to be modified and corresponding original genotype ###
        # ref_base_prop_vec = [1/6, 2/6, 1/6, 2/6], mutation_number=60  in the order of "A", "G", "C", "T"
        # get a vector of number of 4 bases, array([10, 19,  9, 22])
        base_count_AGCT = np.random.multinomial(n=mutations_each_bin_array[bin_index], pvals=ref_base_prop_vec)
        all_chosen_index_base_dict = {}     # {1: "A", 17: "G"}
        for base_i, base_i_count in zip(["A", "G", "C", "T"], base_count_AGCT):
            """ randomly assign base i to passed pos """
            # base_i = "A"
            base_i_index = pos_index_passed[ref_bases_passed_array == base_i]   # index of pos in the bin
            if len(base_i_index) > base_i_count:
                # more positions for mutations
                all_chosen_index_base_dict.update({i: base_i for i in list(np.random.choice(base_i_index, size=base_i_count, replace=False))})
            else:
                # candidate index is not enough
                all_chosen_index_base_dict.update({i: base_i for i in list(base_i_index)})
        """
        Get base that it mutates to 
        """
        all_mutation_info_in_bin = {}  # {pos: {"base":[ref, alt], "ref_count": [], "alt_count": []}}

        for pos_i, ref_base in all_chosen_index_base_dict.items():
            alt_base_list = base_mutation_freq_dict[ref_base]["alt_base_list"]
            alt_freq_list = base_mutation_freq_dict[ref_base]["alt_freq_list"]
            mutated_base = np.random.choice(alt_base_list, 1, alt_freq_list)[0]  # randomly mutated to another base
            all_mutation_info_in_bin.update({pos_i + bin_boundary[0]: {"base": [ref_base, mutated_base],
                                                                       "ref_count": [0] * len(input_bam_list),
                                                                       "alt_count": [0] * len(input_bam_list)}})

        # --------------------------------------------------------------------------------------------------------------

        """ Generate mutation_freq at each pos from mutation_prop_combination_dict """
        # longitudinal_prop_combination_name_list = ["combination1", "combination2", "combination3"]
        # combination_prob_list = [0.2, 0.5, 0.3]; combination_assigned_to_mutation_list = ["combination1", ""]
        combination_assigned_to_mutation_list = np.random.choice(longitudinal_prop_combination_name_list, size=len(all_chosen_index_base_dict), replace=True, p=combination_prob_list)
        # mutation_prop_of_pos_in_bin_array.shape is m * n, n: number of longitudinal sample, m: mutation number
        mutation_prop_of_pos_in_bin_array = np.array([mutation_prop_combination_dict[i]["longitudinal_prop"] for i in combination_assigned_to_mutation_list])

        """ Record reads info to change for each bam """
        for bam_index, pysam_bam in enumerate(pysam_bam_list):
            # in current bin, for each longitudinal sample
            chosen_index_mutation_freq_dict = {chosen_index: mutation_freq for chosen_index, mutation_freq in
                                               zip(sorted(all_chosen_index_base_dict.keys()),
                                                   mutation_prop_of_pos_in_bin_array[:, bam_index])}
            # all_mutation_info_in_bin = {}  # {pos: {"base":[ref, alt], "ref_count": [], "alt_count": []}}
            read_info_to_change_in_a_bam_dict = get_reads_info_to_modify(pysam_bam, pysam_ref, pysam_ref.references[0],
                                                                         bin_boundary[0], bin_boundary[1],
                                                                         chosen_index_mutation_freq_dict,
                                                                         all_mutation_info_in_bin,
                                                                         bam_index)
            # update to total record list
            reads_mutations_record[bam_index].update(read_info_to_change_in_a_bam_dict)

        all_mutation_info.update(all_mutation_info_in_bin)  # pos and reads count
    return reads_mutations_record, all_mutation_info


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


def reads_count_summary(pysam_bam_list, pysam_ref, bin_len):
    """
    Summary reads count in eacg bin of each chr
    :param pysam_bam_list:
    :param pysam_ref:
    :param bin_len:
    :return: {chr1: array([1000, 2000, 50,])}, reads count of each bin
    """
    reads_in_bin_by_chromosome = {}
    # for each chromosome of pysam_ref
    for chromosome, chromosome_len in zip(pysam_ref.references, pysam_ref.lengths):
        reads_bin_chromosome = []   # row: bins, col: start, end, reads of all sample (n+2 cols)
        # for each bin, 10000bp example
        for boundary in get_nodes(chromosome_len, bin_len):
            # get reads of all samples at this bin
            boundary.extend(
                [pysam_sample.count(contig=chromosome, start=boundary[0], stop=boundary[1]) for pysam_sample in
                 pysam_bam_list])
            reads_bin_chromosome.append(boundary)
        # update reads in bins
        reads_in_bin_by_chromosome.update({chromosome: np.array(reads_bin_chromosome)})
    return reads_in_bin_by_chromosome


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
    :param chosen_index_mutation_freq_dict: {1: 0.1, 17: 0.7"}, key is the pos in bin
    :param mutation_info_in_bin: {pos: {"base":[ref, alt], "ref_count": [], "alt_count": []}}, key is the pos in genome
    :return: {"read_name1": {distance_to_end1: mutation_base}}, distance is 0-based
    """
    read_info_to_change_dict = {}   # {"read_name1": {distance_to_end1: mutation_base}}, distance is 0-based

    bin_pileup = pysam_bam.pileup(chr_name, start, end, truncate=True, stepper="samtools", fastafile=pysam_ref)
    for pc in bin_pileup:
        pc_index = pc.reference_pos - start
        if pc.reference_pos in mutation_info_in_bin:
            """ 
            If it is the position need to be mutated 
            """
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
            number_of_mutated_reads = round(chosen_index_mutation_freq_dict[pc_index] * len(read_name_array))
            number_of_not_mutated_reads = len(read_name_array) - number_of_mutated_reads
            mutation_info_in_bin[pc.reference_pos]["ref_count"][index_in_bam_list] += number_of_not_mutated_reads
            mutation_info_in_bin[pc.reference_pos]["alt_count"][index_in_bam_list] += number_of_mutated_reads

            # the index of chosen read in the pc
            chosen_reads_index = np.random.choice(np.array(range(len(read_name_array))), number_of_mutated_reads,
                                                  replace=False)

            mutated_base = mutation_info_in_bin[pc.reference_pos]["base"][1]    # the alt base
            for chosen_index in chosen_reads_index:
                read_name = read_name_array[chosen_index]
                if whether_reverse[chosen_index]:
                    """ read mapped to the reverse strand """
                    if read_name in read_info_to_change_dict:
                        read_info_to_change_dict[read_name].update({pos_to_left_end[chosen_index]: complementary_base_dict[mutated_base]})
                    else:
                        read_info_to_change_dict.update({read_name: {pos_to_left_end[chosen_index]: complementary_base_dict[mutated_base]}})
                else:
                    """ read mapped to the forward strand """
                    if read_name in read_info_to_change_dict:
                        read_info_to_change_dict[read_name].update({pos_to_left_end[chosen_index]: mutated_base})
                    else:
                        read_info_to_change_dict.update({read_name: {pos_to_left_end[chosen_index]: mutated_base}})

    return read_info_to_change_dict


