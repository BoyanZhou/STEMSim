import os

# main function


def map_simulated_reads(unaligned_reads_abs_path_list, bowtie_reference, thread, bowtie2_path, samtools_path, input_type, my_logger):
    """
    :param unaligned_reads_abs_path_list: abs path of input, len of list should be 1 or 2
    :param bowtie_reference:
    :param thread:
    :param bowtie2_path:
    :param samtools_path:
    :param input_type: "fq" or "bam"
    :param my_logger:
    :return:
    """
    if unaligned_reads_abs_path_list[0].endswith("gz"):
        # "XXX.fq.gz"
        abs_path_prefix = ".".join(unaligned_reads_abs_path_list[0].split(".")[:-2])
    else:
        # "XXX.fq" or "XXX.bam"
        abs_path_prefix = ".".join(unaligned_reads_abs_path_list[0].split(".")[:-1])

    print(f"unaligned_reads_abs_path_list is {unaligned_reads_abs_path_list}")
    # """
    ####################
    # single end input #
    ####################
    if len(unaligned_reads_abs_path_list) == 1:
        bowtie2_map_command = ""
        if input_type == "fq":
            bowtie2_map_command = f"{bowtie2_path} -p {thread} -x {bowtie_reference} -U {unaligned_reads_abs_path_list[0]} -S {abs_path_prefix}_aligned.sam"

        elif input_type == "bam":
            samtools_bam_to_fq = f"{samtools_path} bam2fq {unaligned_reads_abs_path_list[0]} > {abs_path_prefix}_bam.fq"
            bowtie2_map_command = f"{bowtie2_path} -p {thread} -x {bowtie_reference} -U {abs_path_prefix}_bam.fq -S {abs_path_prefix}_aligned.sam"
            my_logger.info(samtools_bam_to_fq)
            os.system(samtools_bam_to_fq)
        my_logger.info(bowtie2_map_command)
        os.system(bowtie2_map_command)
    ####################
    # paired end input #
    ####################
    elif len(unaligned_reads_abs_path_list) == 2:
        pass

    samtools_view = f"{samtools_path} view -bS {abs_path_prefix}_aligned.sam > {abs_path_prefix}_aligned.bam"
    samtools_sort = f"{samtools_path} sort {abs_path_prefix}_aligned.bam -o {abs_path_prefix}_aligned_sorted.bam"
    samtools_index = f"{samtools_path} index {abs_path_prefix}_aligned_sorted.bam"

    my_logger.info(samtools_view)
    os.system(samtools_view)
    my_logger.info(samtools_sort)
    os.system(samtools_sort)
    my_logger.info(samtools_index)
    os.system(samtools_index)
    # """
    return f"{abs_path_prefix}_aligned_sorted.bam"


"""
samtools bam2fq Genome1.0.bam > Genome1.0.fq
bowtie2 -p 4 -x /gpfs/data/lilab/home/zhoub03/teddy/ncbi/dbGaP-21556/systematic_simulation/species_reference_genomes/
Bifidobacterium_breve/Bifidobacterium_breve_12L/Bifidobacterium_breve -U Genome1.0.fq -S Genome1.0_aligned.sam
"""