#!/usr/bin/python3
# -*- coding: utf-8 -*-

__author__ = ("Boyan Zhou (boyanzhou1992@gmail.com)")
__version__ = '1.3.0'
__date__ = '13 March 2023'


import argparse
import os


def arg_parsed():
    parser = argparse.ArgumentParser(prog="StemSim",
                                     description=
                                     f"DESCRIPTION\n"
                                     f" StemSim version {__version__} ({__date__}): \n"
                                     f" Short-Term Evolution Mutations Simulator\n\n"
                                     f"AUTHORS: {__author__}\n\n"
                                     f"COMMON COMMANDS\n\n"
                                     f"It depends on Bowtie2, samtools, and CAMISIM \n"
                                     f"\n================= StemSim Process ================= \n\n"
                                     f"Model1: Process the output of CAMISIM \n"
                                     f"Model2: Process raw simulated sequencing data \n",
                                     formatter_class=argparse.RawTextHelpFormatter
                                     )

    parser.add_argument("-m", dest="model", required=True, type=str,
                        choices=["camisim", "reads"],
                        help="choice of working model, must be one of these two options.")
    parser.add_argument('-f', dest="config_file", type=str, help="Config file of input parameters")
    parser.add_argument("-s", dest="subject_id", type=str, help="id of input subject")
    # parser.add_argument('-i', dest="input_dir", type=str, help="Absolute path of input directory")
    parser.add_argument("-o", dest="output_dir", type=str, help="Absolute path of output directory")
    parser.add_argument("-l", dest="log_file", type=str, help="Absolute path of log file")
    parser.add_argument("-p", dest="snpeff_dir", type=str, help="Absolute path of directory of SnpEff")
    parser.add_argument("-a", dest="annotation_file", type=str, help="Absolute path of annotation match table\n"
                                                                     "Each line contains a genome ID and its matching ID in snpEff database.")

    # option group "camisim"
    group1 = parser.add_argument_group("Process CAMISIM output", "Process the output of CAMISIM")
    arg = group1.add_argument
    arg('--input_from_camisim', dest="camisim_output", type=str, help="The output directory of CAMISIM")
    arg('--genome_to_id', dest="genome_to_id", type=str, help="The the absolute path of genome_to_id.tsv file of CAMISIM\n"
                                                              "Each line is a genome_id and a absolute path of fasta separated by Tab")
    parser.add_argument("--nonmutated_list", dest="nonmutated_list", type=str,
                        help="File list of regions where no mutation is generated for each genome\n"
                             "Each line is a genome_id and a absolute path of annotation file separated by Tab\n"
                             "The annotation file should be file in the format of .bed/.gtf/.gff")

    # option group "reads"
    group2 = parser.add_argument_group("Process raw reads", "Process raw reads simulated by other simulators")
    arg = group2.add_argument
    arg('--input_reads', dest="input_reads", type=str, help="Input raw read files from longitudinal samples separated by ':';\n"
                                                            "for paired reads from the same sample, they are separated by ',';\n"
                                                            "example: sample0_R1.fq,sample0_R2.fq:sample1.fq:sample2.fq")
    arg('--reference_genome', dest="genome_ref", type=str, help="The absolute path of reference genome (.fas) from which raw reads are simulated")
    arg('--genome_id', dest="genome_id", type=str, help="The ID of genome from which raw reads are simulated")
    parser.add_argument("--nonmutated_file", dest="nonmutated_file", type=str,
                        help="Absolute path of annotation file containing regions where no mutation should be generated\n"
                             "The annotation file should be file in the format of .bed/.gtf/.gff")

    args = parser.parse_args()
    return args


def file_exists(file):
    # type for checking file exists
    if not os.path.exists(file):
        raise argparse.ArgumentTypeError(f"{file} does not exist.")
    return file


