"""
This script aims to modify the fq file according to the summarize information
"""

import numpy as np
import os


# main function1
def generate_modified_fq(input_fq_path, output_fq_path, read_info_to_change_dict, whether_gzip=True):
    """
    Modify the raw simulated reads according to read_info_to_change_dict, interrogate the name of each read, to check
    whether it is in the read_info_to_change_dict, then change it accordingly
    :param input_fq_path:
    :param output_fq_path:
    :param read_info_to_change_dict: {"read_name1": {distance_2_end1: mutation_base1, distance_2_end2: mutation_base2}}
    {"read_name1": {10: "T", 53: "C"}}. in most cases, only one distance_2_end in a read, 0-based
    :param whether_gzip:
    :return:
    """
    output_fq_f = open(output_fq_path, "w")
    input_fq_f = open(input_fq_path, "r")
    while True:
        line = input_fq_f.readline()
        if not line:
            break   # break at the end of input fastq

        output_fq_f.write(line)                             # 1st line of 4
        read_name = line[1:-1]
        line_of_sequence = input_fq_f.readline()            # 2nd line of 4
        if read_name in read_info_to_change_dict:
            # need to modify
            for distance_to_end, mutation_base in read_info_to_change_dict[read_name].items():
                # replace the original base with the mutated base
                line_of_sequence = line_of_sequence[:distance_to_end] + mutation_base + \
                                   line_of_sequence[distance_to_end + 1:]

        # no need to modify
        output_fq_f.write(line_of_sequence)
        output_fq_f.write(input_fq_f.readline())  # 3rd line of 4
        output_fq_f.write(input_fq_f.readline())  # 4th line of 4

    input_fq_f.close()
    output_fq_f.close()

    if whether_gzip:
        os.system(f"gzip {output_fq_path}")

