#!/usr/bin/python3
# -*- coding: utf-8 -*-

__author__ = ("Boyan Zhou (boyanzhou1992@gmail.com)")
__version__ = '1.0'
__date__ = '01 Nov 2022'


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
                                     f"\n========== StemSim Process ================= \n\n"
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
    # option group "camisim"
    group1 = parser.add_argument_group("Process CAMISIM output", "Process the output of CAMISIM")
    arg = group1.add_argument
    arg('--input_from_camisim', dest="camisim_output", type=str, help="The output directory of CAMISIM")

    # option group "reads"
    group2 = parser.add_argument_group("Process raw reads", "Combine Relative Abundance Of All Sample")
    arg = group2.add_argument
    arg('--input_reads', dest="input_reads", type=str, help="Input raw read files from longitudinal samples separated by ':'; "
                                                            "for paired reads from the same sample, they are separated by ','; "
                                                            "example: sample0_R1.fq,sample0_R2.fq:sample1.fq:sample2.fq")
    arg('--reference_genome', dest="genome_ref", type=str, help="The absolute path of reference genome (.fas) from which raw reads are simulated")
    arg('--genome_id', dest="genome_id", type=str, help="The ID of genome from which raw reads are simulated")

    args = parser.parse_args()
    return args


def file_exists(file):
    # type for checking file exists
    if not os.path.exists(file):
        raise argparse.ArgumentTypeError(f"{file} does not exist.")
    return file


