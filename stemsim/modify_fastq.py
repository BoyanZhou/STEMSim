"""
This script aims to modify the fq file according to the summarize information
"""

import os
import numpy as np
import sys


# main function1 for single fq
def generate_modified_fq(input_fq_path, output_fq_path, read_info_to_change_dict, reads_to_remove_set,
                         whether_gzip=True):
    """
    Modify the raw simulated reads according to read_info_to_change_dict, interrogate the name of each read, to check
    whether it is in the read_info_to_change_dict, then change it accordingly
    :param input_fq_path:
    :param output_fq_path:
    :param read_info_to_change_dict: {"read_name1":
    {"range_list":[[distance_to_end_start, distance_to_end_end]], "alt_list": [mutation_base]}}}
    {"read_name1": {"range_list":[[10, 11], [15, 16]], "alt_list": ["A", "T"]}}.
    in most cases, only one mutation in a read, 0-based
    :param reads_to_remove_set:     {"read_name1", "read_name2", "read_name3"}
    :param whether_gzip:
    :return:
    """
    output_fq_f = open(output_fq_path, "w")
    input_fq_f = open(input_fq_path, "r")
    while True:
        line = input_fq_f.readline()
        if not line:
            break   # break at the end of input fastq

        read_name = line[1:-1]
        line_of_sequence = input_fq_f.readline().strip()        # 2nd line of 4

        if len(reads_to_remove_set) > 0:
            if read_name in reads_to_remove_set:
                # this read should be removed from the raw data
                input_fq_f.readline()       # not output 3rd line of 4
                input_fq_f.readline()       # not output 4th line of 4
                continue

        output_fq_f.write(line)  # output 1st line of 4

        if read_name in read_info_to_change_dict:
            whether_contain_indel = False
            # get the region of bases that do not need to be changed, that is remaining_regions_list
            range_list = read_info_to_change_dict[read_name]["range_list"]
            alt_list = read_info_to_change_dict[read_name]["alt_list"]
            if len(range_list) > 1:
                # read located on reverse strand, mutation should be reversed
                if range_list[0][0] > range_list[-1][0]:
                    range_list = range_list[::-1]
                    alt_list = alt_list[::-1]

            remaining_regions_list, mutation_out_of_bound, whether_indel = get_remaining_region(len(line_of_sequence),
                                                                                                range_list)
            if whether_indel:
                whether_contain_indel = True

            # get the start of updated sequence
            new_sequence = line_of_sequence[remaining_regions_list[0][0]: remaining_regions_list[0][1]]
            # # get the region that do not need to be change and base that need to be replaced
            alt_list_bound = len(alt_list)-mutation_out_of_bound
            for remaining_region, alt in zip(remaining_regions_list[1:], alt_list[:alt_list_bound]):
                if len(alt) > 1:
                    whether_contain_indel = True
                new_sequence += alt + line_of_sequence[remaining_region[0]:remaining_region[1]]

            output_fq_f.write(new_sequence + "\n")          # output 2nd line of 4
            output_fq_f.write(input_fq_f.readline())        # output 3rd line of 4
            # we need to modify the base quality sequence
            if whether_contain_indel:
                base_quality_sequence = input_fq_f.readline().strip()
                base_quality_sequence_new = base_quality_sequence[remaining_regions_list[0][0]: remaining_regions_list[0][1]]
                for remaining_region, alt, ref_range in zip(remaining_regions_list[1:], alt_list[:alt_list_bound], range_list[:alt_list_bound]):
                    if len(alt) == 1:
                        # 1. it is SNV or deletion
                        base_quality_sequence_new += base_quality_sequence[ref_range[0]] + base_quality_sequence[remaining_region[0]:remaining_region[1]]
                    else:
                        # 2. it is an insertion, base quality is sampled from surrounding bases
                        base_quality_array = np.array(list(base_quality_sequence))  # need to remove base close to end
                        # random select qualities from quality of the sequence not close to end
                        added_base_quality = np.random.choice(base_quality_array[5:-5], size=len(alt)-1, replace=True)
                        base_quality_sequence_new += base_quality_sequence[ref_range[0]] + "".join(added_base_quality) + base_quality_sequence[remaining_region[0]: remaining_region[1]]
                output_fq_f.write(base_quality_sequence_new + "\n")     # output 4th line of 4

                if len(new_sequence) != len(base_quality_sequence_new):
                    print(f"Error! The length of read does not match with the length of quality.")
                    print(f"new_sequence is {new_sequence}")
                    print(f"base_quality_sequence_new is {base_quality_sequence_new}")
                    print(f"line_of_sequence is {line_of_sequence}")
                    print(f"range_list is {range_list}")
                    print(f"remaining_regions_list is {remaining_regions_list}")
                    print(f"mutation_out_of_bound is {mutation_out_of_bound}")
                    print(f"whether_indel is {whether_indel}")
                    print(read_info_to_change_dict[read_name])
                    sys.exit()

            else:
                # not contain indel, just output original base quality sequence
                line4 = input_fq_f.readline()
                output_fq_f.write(line4)

                if len(new_sequence) != len(line4.strip()):
                    print(f"Error! The length of read does not match with the length of quality.")
                    print(f"new_sequence is {new_sequence}")
                    print(f"line4 is {line4.strip()}")
                    print(f"line_of_sequence is {line_of_sequence}")
                    print(f"range_list is {range_list}")
                    print(f"remaining_regions_list is {remaining_regions_list}")
                    print(f"mutation_out_of_bound is {mutation_out_of_bound}")
                    print(f"whether_indel is {whether_indel}")
                    print(read_info_to_change_dict[read_name])
                    for remaining_region, alt in zip(remaining_regions_list[1:], alt_list[:-mutation_out_of_bound]):
                        print(remaining_region, alt)
                        if len(alt) > 1:
                            whether_contain_indel = True
                        new_sequence += alt + line_of_sequence[remaining_region[0]:remaining_region[1]]
                        print(new_sequence)
                    sys.exit()

        else:
            # no need to modify
            output_fq_f.write(line_of_sequence + "\n")
            output_fq_f.write(input_fq_f.readline())            # 3rd line of 4
            output_fq_f.write(input_fq_f.readline())            # 4th line of 4

    input_fq_f.close()
    output_fq_f.close()

    if whether_gzip:
        os.system(f"gzip {output_fq_path}")


# main function2 for paired fqs
def generate_modified_paired_fqs(input_fq1_path, input_fq2_path, output_fq1_path, output_fq2_path,
                                 read_info_to_change_dict, reads_to_remove_set, whether_gzip=True):
    """
    Modify the raw simulated reads according to read_info_to_change_dict, interrogate the name of each read, to check
    whether it is in the read_info_to_change_dict, then change it accordingly
    :param input_fq1_path:
    :param input_fq2_path:
    :param output_fq1_path:
    :param output_fq2_path:
    :param read_info_to_change_dict: {"read_name1":
    {"range_list":[[distance_to_end_start, distance_to_end_end]], "alt_list": [mutation_base]}}}
    {"read_name1": {"range_list":[[10, 11], [15, 16]], "alt_list": ["A", "T"]}}.
    in most cases, only one mutation in a read, 0-based
    :param reads_to_remove_set:     {"read_name1", "read_name2", "read_name3"}
    :param whether_gzip:
    :return:
    """
    output_fq1_f = open(output_fq1_path, "w")
    output_fq2_f = open(output_fq2_path, "w")
    output_fq_f_list = [output_fq1_f, output_fq2_f]
    input_fq1_f = open(input_fq1_path, "r")
    input_fq2_f = open(input_fq2_path, "r")
    input_fq_f_list = [input_fq1_f, input_fq2_f]

    while True:
        line1 = input_fq1_f.readline()
        line2 = input_fq2_f.readline()
        if not line1:
            break   # break at the end of input fastq

        read_name1 = line1[1:-1]
        read_name2 = line2[1:-1]
        read_name_list = [read_name1, read_name2]
        sequences_list = [input_fq1_f.readline().strip(), input_fq2_f.readline().strip()]   # 2nd line of 4

        if len(reads_to_remove_set) > 0:
            if read_name1 in reads_to_remove_set or read_name2 in reads_to_remove_set:
                # this read should be removed from the raw data, if one of them in the avoid region
                input_fq1_f.readline()       # not output 3rd line of 4
                input_fq1_f.readline()       # not output 4th line of 4
                input_fq2_f.readline()
                input_fq2_f.readline()
                continue

        output_fq1_f.write(line1)  # output 1st line of 4
        output_fq2_f.write(line2)  # output 1st line of 4

        for fq_index, read_name in enumerate(read_name_list):
            # fq_index = 0, read_name = read_name1
            if read_name in read_info_to_change_dict:
                whether_contain_indel = False
                # get the region of bases that do not need to be changed
                range_list = read_info_to_change_dict[read_name]["range_list"]
                alt_list = read_info_to_change_dict[read_name]["alt_list"]
                if len(range_list) > 1:
                    # read located on reverse strand, mutation should be reversed
                    if range_list[0][0] > range_list[-1][0]:
                        range_list = range_list[::-1]
                        alt_list = alt_list[::-1]

                remaining_regions_list, mutation_out_of_bound, whether_indel = get_remaining_region(
                    len(sequences_list[fq_index]), range_list)
                if whether_indel:
                    whether_contain_indel = True

                # get the updated sequence
                new_sequence = sequences_list[fq_index][remaining_regions_list[0][0]: remaining_regions_list[0][1]]
                # # get the region that do not need to be change and base that need to be replaced
                alt_list_bound = len(alt_list) - mutation_out_of_bound
                for remaining_region, alt in zip(remaining_regions_list[1:], alt_list[:alt_list_bound]):
                    if len(alt) > 1:
                        whether_contain_indel = True
                    new_sequence += alt + sequences_list[fq_index][remaining_region[0]:remaining_region[1]]

                output_fq_f_list[fq_index].write(new_sequence + "\n")                   # output 2nd line of 4
                output_fq_f_list[fq_index].write(input_fq_f_list[fq_index].readline())  # output 3rd line of 4
                # we need to modify the base quality sequence
                if whether_contain_indel:
                    base_quality_sequence = input_fq_f_list[fq_index].readline().strip()
                    base_quality_sequence_new = base_quality_sequence[remaining_regions_list[0][0]: remaining_regions_list[0][1]]
                    for remaining_region, alt, ref_range in zip(remaining_regions_list[1:], alt_list[:alt_list_bound], range_list[:alt_list_bound]):
                        if len(alt) == 1:
                            # 1. it is SNV or deletion
                            base_quality_sequence_new += base_quality_sequence[ref_range[0]] + base_quality_sequence[remaining_region[0]:remaining_region[1]]
                        else:
                            # 2. it is an insertion, base quality is sampled from surrounding bases
                            base_quality_array = np.array(list(base_quality_sequence))      # need to remove base close to end
                            # random select qualities from quality of the sequence not close to end
                            added_base_quality = np.random.choice(base_quality_array[5:-5], size=len(alt) - 1, replace=True)
                            base_quality_sequence_new += base_quality_sequence[ref_range[0]] + "".join(added_base_quality) + base_quality_sequence[remaining_region[0]: remaining_region[1]]
                    output_fq_f_list[fq_index].write(base_quality_sequence_new + "\n")  # output 4th line of 4
                else:
                    # not contain indel, just output original base quality sequence
                    output_fq_f_list[fq_index].write(input_fq_f_list[fq_index].readline())
            else:
                # no need to modify
                output_fq_f_list[fq_index].write(sequences_list[fq_index] + "\n")
                output_fq_f_list[fq_index].write(input_fq_f_list[fq_index].readline())  # 3rd line of 4
                output_fq_f_list[fq_index].write(input_fq_f_list[fq_index].readline())  # 4th line of 4

    output_fq1_f.close()
    output_fq2_f.close()
    input_fq1_f.close()
    input_fq2_f.close()

    if whether_gzip:
        os.system(f"gzip {output_fq1_path}")
        os.system(f"gzip {output_fq2_path}")


def get_remaining_region(sequence_len, avoid_regions_list):
    """
    :param sequence_len: the total length of the sequence, e.g. 150
    :param avoid_regions_list: regions of indel or bases that need to be changed, [[10, 12], [50, 51]]
    :return: list of regions removed avoid_regions_list, [[0, 10], [12, 50], [51, 150]]
    and number of mutations out of bound
    """
    left_boundary = [0]
    right_boundary = []
    mutation_out_of_bound = 0       # record the number of mutations out of bound
    whether_indel = False

    for avoid_region in avoid_regions_list:
        # avoid_region = [51, 150]
        if avoid_region[1] > sequence_len:
            # out of bound
            mutation_out_of_bound += 1
        else:
            if (avoid_region[1] - avoid_region[0]) > 1:
                whether_indel = True
            left_boundary.append(avoid_region[1])
            right_boundary.append(avoid_region[0])
    right_boundary.append(sequence_len)

    remaining_regions_list = [[i, j] for i, j in zip(left_boundary, right_boundary)]
    return remaining_regions_list, mutation_out_of_bound, whether_indel






