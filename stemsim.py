#!/usr/bin/python3
# -*- coding: utf-8 -*-


""" Copyright (c) 2019  Boyan Zhou

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"),
to deal in the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, sublicense,
and/or sell copies of the Software, and to permit persons to whom
the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE
OR OTHER DEALINGS IN THE SOFTWARE.


:Authors: Boyan Zhou
:Contact: boyanzhou1992@gmail.com
:Date: Aug 2022
:Type: tool
:Input: simulated NGS of longitudinal metagenomic data
:Output: vcf and estimations of proportions of identified strains

------------------------------------------------------------------------------------------------------------------------
************************
* Version of software: *
************************
bowtie2         Version: 2.3.5.1
samtools        Version: 1.9

------------------------------------------------------------------------------------------------------------------------
**************
* Update log *
**************
Date:   2022/08/23
By:     Boyan Zhou

"""

import numpy as np
import stemsim
import logging
import time
import sys
import os
import re
import random


def parse_config(config_path):
    all_parameter_dict = {}

    # get all lines of the config file
    with open(config_path, "r") as config_f:
        config_content = config_f.readlines()

    substitution_model_start = False
    # substitution_matrix_q = False
    # base_type_proportion_dict = {}
    base_mutation_freq_dict = {}
    allele_trajectory_dict = {}
    use_manually_fixed_trajectory = False
    fixed_trajectory_start = False
    trajectory_beta_start = False
    short_insertion_beta_start = False
    short_deletion_beta_start = False
    long_insertion_normal_start = False
    long_deletion_normal_start = False

    for line_index, line in enumerate(config_content):
        if line.startswith("#"):
            # proportion_of_original_bases = False
            # substitution_matrix_q = False
            fixed_trajectory_start = False
            trajectory_beta_start = False
            substitution_model_start = False
            short_insertion_beta_start = False
            short_deletion_beta_start = False
            long_insertion_normal_start = False
            long_deletion_normal_start = False
            continue
        if len(line.strip()) == 0:
            continue
        """
        if line.startswith("Subject_ID"):
            all_parameter_dict.update({"Subject_ID": "=".join(line.split("=")[1:]).strip()})
            continue
        if line.startswith("CAMISIM_output_directory"):
            all_parameter_dict.update({"CAMISIM_output_directory": "=".join(line.split("=")[1:]).strip()})
            continue
        """
        if line.startswith("Samtools_path"):
            all_parameter_dict.update({"Samtools_path": "=".join(line.split("=")[1:]).strip()})
            continue
        if line.startswith("Bowtie2_path"):
            all_parameter_dict.update({"Bowtie2_path": "=".join(line.split("=")[1:]).strip()})
            continue
        if line.startswith("N_longitudinal_samples"):
            all_parameter_dict.update({"N_longitudinal_samples": int(line.split("=")[1].strip())})
            continue
        if line.startswith("threads"):
            all_parameter_dict.update({"threads": int(line.split("=")[1].strip())})
            continue
        if line.startswith("Random_seed"):
            all_parameter_dict.update({"Random_seed": int(line.split("=")[1].strip())})
            continue
        if line.startswith("Total_mutations"):
            all_parameter_dict.update({"Total_mutations": int(line.split("=")[1].strip())})
            continue
        if line.startswith("Use_mutation_rate"):
            all_parameter_dict.update({"Use_mutation_rate": str(line.split("=")[1].strip())})   # "True" or "False"
            continue
        if line.startswith("Mutation_rate"):
            all_parameter_dict.update({"Mutation_rate": float(line.split("=")[1].strip())})     # like 5e-05
            continue

        if line.startswith("Substitution_model"):
            substitution_model = line.split("=")[1].strip()
            all_parameter_dict.update({"Substitution_model": substitution_model})
            substitution_model_start = True
            continue

        if substitution_model_start:
            if line.startswith(all_parameter_dict["Substitution_model"]):
                if all_parameter_dict["Substitution_model"] == "JC69":
                    base_mutation_freq_dict.update({"A": {"alt_base_list": ["C", "G", "T"],
                                                          "alt_freq_list": [1/3, 1/3, 1/3]},
                                                    "C": {"alt_base_list": ["A", "G", "T"],
                                                          "alt_freq_list": [1/3, 1/3, 1/3]},
                                                    "G": {"alt_base_list": ["A", "C", "T"],
                                                          "alt_freq_list": [1/3, 1/3, 1/3]},
                                                    "T": {"alt_base_list": ["A", "C", "G"],
                                                          "alt_freq_list": [1/3, 1/3, 1/3]}})

                elif all_parameter_dict["Substitution_model"] == "K80" or \
                        all_parameter_dict["Substitution_model"] == "HKY85":
                    alpha = float(line.split("=")[-1].strip().split("\"")[1])
                    base_mutation_freq_dict.update({"A": {"alt_base_list": ["C", "G", "T"],
                                                          "alt_freq_list": [alpha, 1.0, alpha]},
                                                    "C": {"alt_base_list": ["A", "G", "T"],
                                                          "alt_freq_list": [alpha, alpha, 1.0]},
                                                    "G": {"alt_base_list": ["A", "C", "T"],
                                                          "alt_freq_list": [1.0, alpha, alpha]},
                                                    "T": {"alt_base_list": ["A", "C", "G"],
                                                          "alt_freq_list": [alpha, 1.0, alpha]}})
                elif all_parameter_dict["Substitution_model"] == "TN93":
                    alpha_transversion = float(line.split(",")[1].split("=")[-1].strip().split("\"")[1])
                    alpha4 = float(line.split(",")[-1].split("=")[-1].strip().split("\"")[1])
                    base_mutation_freq_dict.update({"A": {"alt_base_list": ["C", "G", "T"],
                                                          "alt_freq_list": [alpha_transversion, 1.0, alpha_transversion]},
                                                    "C": {"alt_base_list": ["A", "G", "T"],
                                                          "alt_freq_list": [alpha_transversion, alpha_transversion, alpha4]},
                                                    "G": {"alt_base_list": ["A", "C", "T"],
                                                          "alt_freq_list": [1.0, alpha_transversion, alpha_transversion]},
                                                    "T": {"alt_base_list": ["A", "C", "G"],
                                                          "alt_freq_list": [alpha_transversion, alpha4, alpha_transversion]}})

                elif all_parameter_dict["Substitution_model"] == "REV":
                    alpha1, alpha2, alpha3, alpha4, alpha5 = [float(i.split("=")[-1].strip().split("\"")[1]) for i in
                                                              line.split(",")]
                    base_mutation_freq_dict.update({"A": {"alt_base_list": ["C", "G", "T"],
                                                          "alt_freq_list": [alpha1, 1.0, alpha2]},
                                                    "C": {"alt_base_list": ["A", "G", "T"],
                                                          "alt_freq_list": [alpha1, alpha3, alpha4]},
                                                    "G": {"alt_base_list": ["A", "C", "T"],
                                                          "alt_freq_list": [1.0, alpha3, alpha5]},
                                                    "T": {"alt_base_list": ["A", "C", "G"],
                                                          "alt_freq_list": [alpha2, alpha4, alpha5]}})
                # normalize the sum of proportion to 1
                for base in list(base_mutation_freq_dict.keys()):
                    alt_freq_list_temp = base_mutation_freq_dict[base]["alt_freq_list"]
                    base_mutation_freq_dict[base]["alt_freq_list"] = [i/sum(alt_freq_list_temp) for i in alt_freq_list_temp]

        if line.startswith("Use_manually_fixed_trajectory"):
            if line.split("=")[1].strip() == "True":
                use_manually_fixed_trajectory = True

        #################################
        # use manually fixed trajectory #
        #################################
        if line.startswith("Manually_fixed_trajectory"):
            if use_manually_fixed_trajectory:
                fixed_trajectory_start = True
                continue
        if fixed_trajectory_start:
            cols = line.strip().split("\t")
            if len(cols) == 3:
                allele_trajectory_dict.update(
                    {cols[0]: {"prob": float(cols[1]), "longitudinal_prop": [float(i) for i in cols[2].split(",")]}})
            continue

        #########################################
        # use trajectory from beta distribution #
        #########################################
        if line.startswith("Trajectory_from_beta_distribution"):
            if not use_manually_fixed_trajectory:
                trajectory_beta_start = True
                continue
        if trajectory_beta_start:
            cols = line.strip().split("\t")
            if len(cols) == 5:
                allele_frequency_list = np.random.beta(float(cols[2]), float(cols[3]),
                                                       all_parameter_dict["N_longitudinal_samples"])
                if cols[4] == "increase":
                    allele_frequency_list = sorted(allele_frequency_list)
                elif cols[4] == "decrease":
                    allele_frequency_list = sorted(allele_frequency_list, reverse=True)
                allele_trajectory_dict.update({cols[0]: {"prob": float(cols[1]),
                                                         "longitudinal_prop": allele_frequency_list}})
        ########################
        # get indel parameters #
        ########################
        # get short insertion parameter
        if line.startswith("Generate_short_insertion"):
            if line.split("=")[1].strip() == "True":
                all_parameter_dict.update({"Generate_short_insertion": True})
            else:
                all_parameter_dict.update({"Generate_short_insertion": False})
            continue
        if line.startswith("Short_insertion_number"):
            all_parameter_dict.update({"Short_insertion_number": int(line.split("=")[1].strip())})
            continue
        if line.startswith("Short_insertion_use_mutation_rate"):
            if line.split("=")[1].strip() == "True":
                all_parameter_dict.update({"Short_insertion_use_mutation_rate": True})
            else:
                all_parameter_dict.update({"Short_insertion_use_mutation_rate": False})
            continue
        if line.startswith("Short_insertion_mutation_rate"):
            all_parameter_dict.update({"Short_insertion_mutation_rate": float(line.split("=")[1].strip())})
            continue
        if line.startswith("Short_insertion_from_beta_distribution"):
            short_insertion_beta_start = True
            continue
        if short_insertion_beta_start:
            cols = line.strip().split("\t")
            if len(cols) == 2:
                all_parameter_dict.update({"Short_insertion_beta_parameter": [float(cols[0]), float(cols[1])]})
            continue

        # get short deletion parameter
        if line.startswith("Generate_short_deletion"):
            if line.split("=")[1].strip() == "True":
                all_parameter_dict.update({"Generate_short_deletion": True})
            else:
                all_parameter_dict.update({"Generate_short_deletion": False})
            continue
        if line.startswith("Short_deletion_number"):
            all_parameter_dict.update({"Short_deletion_number": int(line.split("=")[1].strip())})
            continue
        if line.startswith("Short_deletion_use_mutation_rate"):
            if line.split("=")[1].strip() == "True":
                all_parameter_dict.update({"Short_deletion_use_mutation_rate": True})
            else:
                all_parameter_dict.update({"Short_deletion_use_mutation_rate": False})
            continue
        if line.startswith("Short_deletion_mutation_rate"):
            all_parameter_dict.update({"Short_deletion_mutation_rate": float(line.split("=")[1].strip())})
            continue
        if line.startswith("Short_deletion_from_beta_distribution"):
            short_deletion_beta_start = True
            continue
        if short_deletion_beta_start:
            cols = line.strip().split("\t")
            if len(cols) == 2:
                all_parameter_dict.update({"Short_deletion_beta_parameter": [float(cols[0]), float(cols[1])]})
            continue

        # get long insertion parameter
        if line.startswith("Generate_long_insertion"):
            if line.split("=")[1].strip() == "True":
                all_parameter_dict.update({"Generate_long_insertion": True})
            else:
                all_parameter_dict.update({"Generate_long_insertion": False})
            continue
        if line.startswith("Long_insertion_number"):
            all_parameter_dict.update({"Long_insertion_number": int(line.split("=")[1].strip())})
            continue
        if line.startswith("Long_insertion_from_normal_distribution"):
            long_insertion_normal_start = True
            continue
        if long_insertion_normal_start:
            cols = line.strip().split("\t")
            if len(cols) == 2:
                all_parameter_dict.update({"Long_insertion_normal_parameter": [float(cols[0]), float(cols[1])]})
            continue
        # get long deletion parameter
        if line.startswith("Generate_long_deletion"):
            if line.split("=")[1].strip() == "True":
                all_parameter_dict.update({"Generate_long_deletion": True})
            else:
                all_parameter_dict.update({"Generate_long_deletion": False})
            continue
        if line.startswith("Long_deletion_number"):
            all_parameter_dict.update({"Long_deletion_number": int(line.split("=")[1].strip())})
            continue
        if line.startswith("Long_deletion_from_normal_distribution"):
            long_deletion_normal_start = True
            continue
        if long_deletion_normal_start:
            cols = line.strip().split("\t")
            if len(cols) == 2:
                all_parameter_dict.update({"Long_deletion_normal_parameter": [float(cols[0]), float(cols[1])]})
            continue

    # {"A": {"alt_base_list": ["G", "C", "T"], "alt_freq_list": [0.25, 0.5, 0.25]}, ...}
    all_parameter_dict.update({"Matrix_Q": base_mutation_freq_dict})
    # {"combination1": {"prob": 0.3, "longitudinal_prop": [0.1, 0.7, 0.8]}}
    all_parameter_dict.update({"Trajectory": allele_trajectory_dict})
    return all_parameter_dict


def get_genome_id_fas_dict(genome_to_id_path):
    genome_id_fas_dict = {}
    with open(genome_to_id_path, "r") as genome_to_id_f:
        for line in genome_to_id_f:
            cols = line.strip().split("\t")
            genome_id_fas_dict.update({cols[0]: cols[1]})
    return genome_id_fas_dict


def main():
    start_time = time.time()
    # get parsed arguments
    options = stemsim.argument_parse.arg_parsed()

    ################
    # set log file #
    ################
    logging.basicConfig(filename=options.log_file, format='%(asctime)s\t%(levelname)s\t%(name)s: %(message)s',
                        level=logging.DEBUG)
    handler = logging.StreamHandler()
    handler.setLevel(logging.INFO)
    logging.getLogger().addHandler(handler)

    my_logger = logging.getLogger("main")
    my_logger.info("Started with the command: " + " ".join(sys.argv) + "\n")

    parameters_dict = parse_config(options.config_file)
    random.seed(parameters_dict["Random_seed"])

    """ Annotation by SnpEff """
    geno_id_snpeff_id_dict = {}     # {"genome_id": "standard_id_in_snpeff_database"}
    if options.annotation_file:
        with open(options.annotation_file, "r") as annotation_match_f:
            for line in annotation_match_f:
                cols = line.strip().split("\t")
                geno_id_snpeff_id_dict.update({cols[0]: cols[1]})
    snpeff_dir = ""
    if options.snpeff_dir:
        # check the existence of snpeff
        snpeff_dir = options.snpeff_dir
        snpeff_path = os.path.join(snpeff_dir, "snpEff.jar")
        if not os.path.exists(snpeff_path):
            print(f"Error! snpEff.jar is not found under the directory {snpeff_dir}. "
                  f"Please install snpEff under the given path if you want to annotate. "
                  f"Otherwise, please do not use -a and -p option.")
            sys.exit()

    # *********************************
    # model1: CAMISIM output as input *
    # *********************************
    if options.model == "camisim":
        # camisim_all_files = os.listdir(parameters_dict["CAMISIM_output_directory"])
        camisim_all_files = os.listdir(options.camisim_output)

        # ## get longitudinal sample directories ###
        # camisim_all_subdir = ["2022.11.01_16.52.10_sample_0", "2022.11.01_16.52.10_sample_1", "distributions"]
        camisim_all_subdir = [i for i in camisim_all_files if os.path.isdir(os.path.join(options.camisim_output, i))]
        camisim_sample_dir_pattern = re.compile(".+_sample_")
        # like 2022.11.01_16.55.37_sample_0
        camisim_sample_dirs = sorted([i for i in camisim_all_subdir if
                                      re.search(camisim_sample_dir_pattern, i, flags=0)])
        if len(camisim_sample_dirs) == 0:
            print(f"Error! No sample directory is found under the output directory of camisim {options.camisim_output}")
            sys.exit()
        if len(camisim_sample_dirs) == len(parameters_dict["N_longitudinal_samples"]):
            print(f"Error! The number of longitudinal samples is {parameters_dict['N_longitudinal_samples']},"
                  f"does not match with the number of samples by CAMISIM. Please change N_longitudinal_samples in "
                  f"config.txt")
            sys.exit()

        # ## get genome_id_fas_dict ###
        # genome_to_id_path = os.path.join(options.output_dir, "genome_to_id.tsv")
        genome_id_fas_dict = get_genome_id_fas_dict(options.genome_to_id)
        # print(f"The genome id and reference fas are {genome_id_fas_dict}")
        genome_id_list = sorted(list(genome_id_fas_dict.keys()))

        # ## get genome_id_read_dict ###
        # for CAMISIM we use unaligned bam, {"geno0": [bam0_path, bam1_path, bam2_path]}
        genome_id_read_dict = {i: [] for i in genome_id_list}
        for camisim_sample_dir in camisim_sample_dirs:
            unaligned_bam_dir = os.path.join(options.camisim_output, camisim_sample_dir, "bam")
            for genome_id in genome_id_list:
                genome_id_read_dict[genome_id].append([os.path.join(unaligned_bam_dir, f"{genome_id}.bam")])
        print(genome_id_read_dict)

        # ## get genome_id_nonmutated_region
        # each line is: genome_id   path_of_annotation_file (gff, gtf, bed)
        genome_id_nonmutated_file_dict = {}     # "geno1": "path\nonmutated_region.gff"
        if options.nonmutated_list:
            with open(options.nonmutated_list, "r") as nonmutated_list:
                for line in nonmutated_list:
                    if line.startswith("#"):
                        continue
                    cols = line.strip().split("\t")
                    # check whether geno_id in the genome_id_list above
                    if cols[0] in genome_id_list:
                        if os.path.exists(cols[1]):
                            genome_id_nonmutated_file_dict.update({cols[0]: cols[1]})
                        else:
                            my_logger.info(f"Warning! The annotation file {cols[1]} of genome {cols[0]} does not exist."
                                           f" Skip!")
                    else:
                        # this geno_id is not in previous file
                        my_logger.info(f"Warning! The genome {cols[0]} does not exist in genome_to_id file. Skip!")

        # __init__(self, model, subject_id, genome_id_list, genome_id_fas_dict, genome_id_read_dict, bowtie2_path,
        # samtools_path, input_type, thread, my_logger)
        input_type = "bam"      # because we use CAMISIM, otherwise usually fq
        subject_microbial_community = stemsim.longitudinal_samples.Subject(options.model, options.subject_id,
                                                                           genome_id_list, genome_id_fas_dict,
                                                                           genome_id_read_dict,
                                                                           genome_id_nonmutated_file_dict,
                                                                           parameters_dict["Bowtie2_path"],
                                                                           parameters_dict["Samtools_path"],
                                                                           input_type,
                                                                           parameters_dict["Total_mutations"],
                                                                           parameters_dict["Mutation_rate"],
                                                                           parameters_dict["Substitution_model"],
                                                                           parameters_dict["threads"],
                                                                           geno_id_snpeff_id_dict, snpeff_dir,
                                                                           my_logger)

    # ********************************************
    # model2: other simulated raw reads as input *
    # ********************************************
    elif options.model == "reads":
        read_sample_list = [i.split(",") for i in options.input_reads.split(":")]
        genome_id_read_dict = {options.genome_id: read_sample_list}
        genome_id_fas_dict = {options.genome_id: options.genome_ref}
        genome_id_nonmutated_file_dict = {}

        if len(read_sample_list) != parameters_dict["N_longitudinal_samples"]:
            print(f"Error! The number of longitudinal samples is {parameters_dict['N_longitudinal_samples']},"
                  f"does not match with the number of samples ({len(read_sample_list)}) "
                  f"in provided raw reads. Please change N_longitudinal_samples in config.txt")
            sys.exit()

        if options.nonmutated_file:
            genome_id_nonmutated_file_dict = {options.genome_id: options.nonmutated_file}
        if read_sample_list[0][0].endswith(".bam"):
            input_type = "bam"
        else:
            input_type = "fq"
        subject_microbial_community = stemsim.longitudinal_samples.Subject(options.model, options.subject_id,
                                                                           [options.genome_id], genome_id_fas_dict,
                                                                           genome_id_read_dict,
                                                                           genome_id_nonmutated_file_dict,
                                                                           parameters_dict["Bowtie2_path"],
                                                                           parameters_dict["Samtools_path"],
                                                                           input_type,
                                                                           parameters_dict["Total_mutations"],
                                                                           parameters_dict["Mutation_rate"],
                                                                           parameters_dict["Substitution_model"],
                                                                           parameters_dict["threads"],
                                                                           geno_id_snpeff_id_dict, snpeff_dir,
                                                                           my_logger)
    else:
        print("Error! The working model must be camisim or reads.")
        sys.exit()

    # add the parameter of short/long indel
    if parameters_dict["Generate_short_insertion"]:
        alpha, beta = parameters_dict["Short_insertion_beta_parameter"]
        subject_microbial_community.add_insertion_parameter(alpha, beta, parameters_dict["Short_insertion_number"],
                                                            parameters_dict["Short_insertion_mutation_rate"],
                                                            parameters_dict["Short_insertion_use_mutation_rate"])
    if parameters_dict["Generate_short_deletion"]:
        alpha, beta = parameters_dict["Short_deletion_beta_parameter"]
        subject_microbial_community.add_deletion_parameter(alpha, beta, parameters_dict["Short_deletion_number"],
                                                           parameters_dict["Short_deletion_mutation_rate"],
                                                           parameters_dict["Short_deletion_use_mutation_rate"])
    if parameters_dict["Generate_long_insertion"]:
        mu, sigma = parameters_dict["Long_insertion_normal_parameter"]
        subject_microbial_community.add_long_insertion_parameter(mu, sigma, parameters_dict["Long_insertion_number"])
    if parameters_dict["Generate_long_deletion"]:
        mu, sigma = parameters_dict["Long_deletion_normal_parameter"]
        subject_microbial_community.add_long_deletion_parameter(mu, sigma, parameters_dict["Long_deletion_number"])

    subject_microbial_community.check_bowtie_ref()
    subject_microbial_community.align_simulated_reads()
    subject_microbial_community.calculate_ref_base_prop()
    if not os.path.exists(options.output_dir):
        os.system(f"mkdir {options.output_dir}")
    if not os.path.exists(options.output_dir):
        print(f"Error! The output folder {options.output_dir} does not exist. Please create the folder first.")
        sys.exit()
    subject_microbial_community.generate_mutations(parameters_dict["Trajectory"],
                                                   parameters_dict["Use_mutation_rate"],
                                                   parameters_dict["Substitution_model"],
                                                   parameters_dict["Matrix_Q"],
                                                   options.output_dir,
                                                   options.model)

    end_time = time.time()
    my_logger.info(f"It totally took {end_time-start_time}s. End of this job.")


if __name__ == "__main__":
    sys.exit(main())
