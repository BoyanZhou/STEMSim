"""
For pc, get distance to end of reads, 0-based
Read mapped to the forward strand, distance to the left end; mapped to the reverse strand, distance to the right end
Which means they are all distance to the 5' in raw fq sequence
"""

import numpy as np


# main function
def summarize_distance(pc):
    """
    Accept a pysam pileup column
    :param pc:
    :return: the distance array to the left end (mapped to forward strand), to the right end (mapped to reverse strand);
    0 based (soft clipped bases counted); and the array indicate whether mapped to reverse strand
    """
    # this left end is the left end shown in tview, not the real left end in fq
    pos_to_left_end = np.array(pc.get_query_positions())  # 0-based, array, to left end, soft clipped bases counted
    whether_reverse = []            # whether each read is mapped to the reverse strand of the genome
    reads_length = []               # the real length of each reads
    for read in pc.pileups:
        whether_reverse.append(read.alignment.is_reverse)
        reads_length.append(read.alignment.query_length)
    pos_to_right_end = np.array(reads_length) - pos_to_left_end - 1     # 0-based
    whether_reverse_array = np.array(whether_reverse)

    """ for reads on reverse strand, use the distance to another end """
    pos_to_real_left_end = pos_to_left_end.copy()
    # print(f"pos_to_real_left_end is {pos_to_real_left_end}")
    # print(f"whether_reverse_array is {whether_reverse_array}")
    if len(pos_to_real_left_end) > 0:
        pos_to_real_left_end[whether_reverse_array] = pos_to_right_end[whether_reverse_array]
    return pos_to_real_left_end, whether_reverse_array      # two array




