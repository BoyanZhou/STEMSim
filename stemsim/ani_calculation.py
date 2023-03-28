"""
This script calculates ANI between pairwise longitudinal samples, for each given reference genome

ConANI substitutions are called whenever the consensus base differs,
and popANI substitutions are only called when there is no allelic overlap between samples.

Reference:
Olm, M. R., Crits-Christoph, A., Bouma-Gregson, K., Firek, B. A., Morowitz, M. J., & Banfield, J. F. (2021).
inStrain profiles population microdiversity from metagenomic data and sensitively detects shared microbial strains.
Nature Biotechnology, 39(6), 727-736.

"""
import numpy as np


def calculate_ani(mutation_info_dict, chr_name_list, chr_len_list, sample_name_list, abs_path_output):
    """
    In calculation of ANI, we only consider SNV
    :param mutation_info_dict: {chr_name: {{position: {"type": "deletion", "base": [ref_base, mutated_base],
                                    "ref_count": [0] * len(input_bam_list), "alt_count": [0] * len(input_bam_list)}}}}
    :param chr_name_list:
    :param chr_len_list:
    :param sample_name_list: [sample1, ...]
    :param abs_path_output: output ANI file
    :return:
    """
    conANI_count = np.zeros([len(sample_name_list), len(sample_name_list)], dtype=int)  # count qualified sites
    popANI_count = np.zeros([len(sample_name_list), len(sample_name_list)], dtype=int)  # count qualified sites
    for chr_name in chr_name_list:
        if chr_name in mutation_info_dict:
            mutation_info_chr_dict = mutation_info_dict[chr_name]
            for pos, mutation_info in mutation_info_chr_dict.items():
                # {"type": "deletion", "base": [ref_base, mutated_base],
                # "ref_count": [0] * len(input_bam_list), "alt_count": [0] * len(input_bam_list)}
                if mutation_info["type"] == "SNV":
                    ref_count_array = np.array(mutation_info["ref_count"])
                    alt_count_array = np.array(mutation_info["alt_count"])
                    # ref_count_array = np.array([10, 50, 60, 130, 30]),
                    # alt_count_array = np.array([70, 0, 300, 30, 230])
                    total_count_array = ref_count_array + alt_count_array + 0.0001  # avoid divide 0
                    ref_prop_array = ref_count_array/total_count_array
                    ref_if_dominate_array = ref_prop_array >= 0.5
                    ref_if_nondominate_array = ~ref_if_dominate_array
                    # conANI is where ref_if_dominate_array is different ref_if_dominate_array
                    # bool matrix, indicate whether conANI between samples
                    conANI_bool_matrix = np.outer(ref_if_dominate_array, ref_if_nondominate_array) + np.outer(
                        ref_if_nondominate_array, ref_if_dominate_array)
                    conANI_count += conANI_bool_matrix
                    # popANI is where base completely substitute
                    ref_if_0_array = ref_prop_array == 0
                    ref_if_1_array = ref_prop_array == 1
                    popANI_bool_matrix = np.outer(ref_if_0_array, ref_if_1_array) + np.outer(ref_if_1_array,
                                                                                             ref_if_0_array)
                    popANI_count += popANI_bool_matrix
    chr_len_sum = sum(chr_len_list)
    conANI = (1 - conANI_count / chr_len_sum).astype(str)
    popANI = (1 - popANI_count / chr_len_sum).astype(str)
    with open(abs_path_output, "w") as ani_f:
        ani_f.write("conANI\t" + "\t".join(sample_name_list) + "\n")
        for sample_index, sample_name in enumerate(sample_name_list):
            ani_f.write(sample_name + "\t" + "\t".join(conANI[sample_index]) + "\n")
        ani_f.write("\n")
        ani_f.write("popANI\t" + "\t".join(sample_name_list) + "\n")
        for sample_index, sample_name in enumerate(sample_name_list):
            ani_f.write(sample_name + "\t" + "\t".join(popANI[sample_index]) + "\n")
