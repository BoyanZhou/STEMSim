import numpy as np
import stemsim.mapping as reads_mapping
import stemsim.reads_info_to_modify as ritm
import stemsim.modify_fastq as mf
import os
import pysam
from collections import Counter


class Subject:
    """
    Subject has multiple microbiome samples (from CAMISIM) at several time points (or repeated measures)
    last update 11/14/2022
    """

    def __init__(self, model, subject_id, genome_id_list, genome_id_fas_dict, genome_id_read_dict, bowtie2_path,
                 samtools_path, input_type, substitution_model, thread, my_logger):
        self.model = model                                              # "camisim" or "reads"
        self.ID = subject_id
        self.genome_id = genome_id_list                                 # [geno0, geno1, geno2]
        self.genome_id_fas_dict = genome_id_fas_dict                    # {"geno0": "/data/genomes/ref1.fas"}
        self.genome_id_bowtie_ref_dict = {genome_id: ".".join(fas.split(".")[:-1]) for genome_id, fas in
                                          genome_id_fas_dict.items()}   # {"geno0": "/data/genomes/ref1"}
        self.genome_id_read_dict = genome_id_read_dict                  # {"geno0": [[fq0_path], [fq1_path], [fq2_path]]} or unaligned bam
        self.genome_id_aligned_bam_dict = {}                            # {"geno0": [bam0_path, bam1_path, bam2_path]}
        self.genome_id_base_prop = {}                                   # {"geno0": {"A": 1/4, "G": 1/4, "C": 1/4, "T": 1/4}}
        self.thread = thread
        self.bowtie2_path = bowtie2_path
        self.samtools_path = samtools_path
        self.input_type = input_type                                    # "fq" or "bam"
        self.substitution_model = substitution_model    # must be JC69, K80, HKY85, TN93, or REV
        self.logger = my_logger

    def check_bowtie_ref(self):
        for genome_id, fas_path in self.genome_id_fas_dict.items():
            # fas_dir, fas_name = os.path.split(fas_path)
            fas_prefix = ".".join(fas_path.split(".")[:-1])     # /path_to/Bifidobacterium_breve_12L
            if not os.path.exists(f"{fas_prefix}.1.bt2"):
                # need to generate the bowtie reference
                build_bowtie_ref = f"{self.bowtie2_path}-build {fas_path} {fas_prefix}"
                os.system(build_bowtie_ref)
                self.logger.info(build_bowtie_ref)

    def calculate_ref_base_prop(self):
        """
        :return: the dict of base proportion, like geno_id: {"A": 1/4, "G": 1/4, "C": 1/4, "T": 1/4}
        """
        base_frequency_dict = {}
        for genome_id in self.genome_id:
            self.logger.info(f"Calculate substitution frequency for reference of {genome_id} ... ...")
            base_count_dict = {"A": 0, "C": 0, "G": 0, "T": 0}
            pysam_ref = pysam.FastaFile(self.genome_id_fas_dict[genome_id])
            seq1 = pysam_ref.fetch(reference=pysam_ref.references[0])
            base_count_original = Counter(seq1)
            for base in ["A", "C", "G", "T"]:
                if base in base_count_original:
                    base_count_dict[base] += base_count_original[base]
                if base.lower() in base_count_original:
                    base_count_dict[base] += base_count_original[base.lower()]

            self.logger.info(f"The number of four bases in genome {genome_id} is {base_count_dict}")
            four_bases_total_count = sum(list(base_count_dict.values()))    # total number of A G C T, excluding N
            if four_bases_total_count == 0:
                # no A, G, C, T is found
                base_frequency_dict.update({genome_id: None})
            else:
                # calculate frequency of four bases excluding N
                base_frequency_dict.update({genome_id: {i: base_count_dict[i]/four_bases_total_count for i in ["A", "C", "G", "T"]}})
        self.genome_id_base_prop = base_frequency_dict

    def align_simulated_reads(self):
        """
        Step1: aligns the simulated reads (fq or unaligned bam) to corresponding original fas, get aligned bam list
        :return:
        """
        for genome_id in self.genome_id:
            self.logger.info(f"Processing the simulated reads from genome {genome_id} ... ...")
            bowtie_ref = self.genome_id_bowtie_ref_dict[genome_id]
            aligned_bam_path_list = []
            ########################
            # alignment processing #
            ########################
            for unaligned_reads_abs_path_list in self.genome_id_read_dict[genome_id]:
                # unaligned_reads_abs_path_list = [fq0_path] or [fq0_0_path, fq0_1_path]
                aligned_bam_path = reads_mapping.map_simulated_reads(unaligned_reads_abs_path_list, bowtie_ref,
                                                                     self.thread, self.bowtie2_path, self.samtools_path,
                                                                     self.input_type, self.logger)
                aligned_bam_path_list.append(aligned_bam_path)

            # store the path of aligned bam
            self.genome_id_aligned_bam_dict.update({genome_id: aligned_bam_path_list})

    def generate_mutations(self, mutation_traj_combination_dict, mutation_number, substitution_model,
                           original_q_matrix, output_dir, working_model, pool_reads=True):
        """
        Step2: Generate reads info that needs to be modified
        @ like [1/6, 2/6, 1/6, 2/6], in the order of "A", "C", "G", "T"
        :param mutation_traj_combination_dict: {"combination1": {"prob": 0.3, "longitudinal_prop": [0.1, 0.7, 0.8]}}
        :param mutation_number: total number of mutation to generate
        :param substitution_model: JC69, K80, HKY85, TN93, or REV
        :param original_q_matrix: {"A": {"alt_base_list": ["C", "G", "T"], "alt_freq_list": [0.25, 0.5, 0.25]}, ...}
        :param working_model: camisim or reads
        :param pool_reads: whether combine the generated fq.gz
        :return:
        """
        genome_id_modified_fq_dict = {}     # {"geno_id1": [, , ,]}
        for genome_id in self.genome_id:
            self.logger.info(f"Generate mutations for reads from {genome_id} ... ...")
            input_bam_list = self.genome_id_aligned_bam_dict[genome_id]
            ref_fas_path = self.genome_id_fas_dict[genome_id]
            # record all created mutations, {"read_name1": {distance_to_end1: mutation_base}}, distance is 0-based

            ###########################
            # get four bases mutation #
            ###########################
            # according to give parameters and proportion of four types of bases
            if self.genome_id_base_prop[genome_id]:
                # proportion of four bases exists
                ref_base_prop_dict = self.genome_id_base_prop[genome_id]    # {"A": 1/4}
                if substitution_model == "JC69" or substitution_model == "K80":
                    # in the order of "A", "C", "G", "T"
                    ref_base_prop_vec = [0.25, 0.25, 0.25, 0.25]
                else:
                    ref_base_prop_vec = [ref_base_prop_dict[i] for i in ["A", "C", "G", "T"]]

            else:
                self.logger.info(f"Warning! A, G, C, T are not contained in the genome of {genome_id}, skip this.")
                continue
            #####################################
            """ mutation information in reads """
            #####################################
            # reads_mutations_record is list, length is the number of longitudinal samples
            # all_mutation_info is dict, {pos: {"base":[ref, alt], "ref_count": [], "alt_count": []}}
            reads_mutations_record, all_mutation_info = ritm.summarize_reads_to_modify(input_bam_list, ref_fas_path,
                                                                                       ref_base_prop_vec,
                                                                                       mutation_traj_combination_dict,
                                                                                       mutation_number,
                                                                                       original_q_matrix,
                                                                                       depth_threshold=5,
                                                                                       bin_len=10000)

            """
            Output true mutations
            """
            self.output_truth_of_mutation(all_mutation_info, output_dir, genome_id)
            # ----------------------------------------------------------------------------------------------------------
            ################################
            """ Step3: modified fq files """
            ################################
            genome_id_modified_fq_dict.update({genome_id: []})      # record the path of modified fqs
            input_fq_path_list = self.genome_id_read_dict[genome_id]    # [[], [], []]
            for fq_index in range(len(input_fq_path_list)):
                # for each longitudinal sample
                processed_fq_path = []      # [XXX_1.fq, XXX_2.fq] or [XXX.fq]
                """ Preprocess fq according to different types """
                input_fq_path = input_fq_path_list[fq_index]    # can be list of fq or fq.gz or unaligned bam
                if input_fq_path[0].endswith(".gz"):
                    for input_fq_gz in input_fq_path:
                        os.system(f"gunzip -c {input_fq_gz} > {input_fq_gz[:-3]}")
                        processed_fq_path.append(input_fq_gz[:-3])

                elif input_fq_path[0].endswith(".bam"):
                    processed_fq_path = [i[:-4] + "_bam.fq" for i in input_fq_path]
                else:
                    processed_fq_path = input_fq_path_list[fq_index]

                """ Output modified fq """
                if not os.path.exists(output_dir):
                    os.system(f"mkdir {output_dir}")
                if len(processed_fq_path) == 1:
                    output_fq_path = os.path.join(output_dir, f"{self.ID}_{genome_id}_t{fq_index}_mutated.fq")
                    mf.generate_modified_fq(processed_fq_path[0], output_fq_path, reads_mutations_record[fq_index],
                                            whether_gzip=True)
                    genome_id_modified_fq_dict[genome_id].append([f"{output_fq_path}.gz"])

                elif len(processed_fq_path) == 2:
                    fq_gz_list_temp = []
                    for fq_i, input_fq_i_path in enumerate(sorted(processed_fq_path)):
                        # fq_i is 0 or 1, represent R1 or R2 in paired fqs
                        output_fq_path = os.path.join(output_dir,
                                                      f"{self.ID}_{self.genome_id}_t{fq_index}_mutated_R{fq_i + 1}.fq")
                        mf.generate_modified_fq(input_fq_i_path, output_fq_path, reads_mutations_record[fq_index],
                                                whether_gzip=True)
                        fq_gz_list_temp.append(f"{output_fq_path}.gz")
                    genome_id_modified_fq_dict[genome_id].append(fq_gz_list_temp)
        # --------------------------------------------------------------------------------------------------------------
        # whether pool the modified reads from each genome
        if pool_reads and working_model == "camisim":
            print(f"\tgenome_id_modified_fq_dict is {genome_id_modified_fq_dict}\t")
            fq_gz_list_temp2 = list(genome_id_modified_fq_dict.values())[0]     # [[t0_fq], [t1_fq], [t2_fq]]
            if len(fq_gz_list_temp2[0]) == 1:
                ########################################################
                """ check the first sample, whether single end reads """
                ########################################################
                longitudinal_sample_lists = [[] for i in range(len(fq_gz_list_temp2))]  # genome ids at all time point
                for fq_gz_list_temp3 in list(genome_id_modified_fq_dict.values()):
                    for sample_index, fq_gz in enumerate(fq_gz_list_temp3):
                        longitudinal_sample_lists[sample_index].append(fq_gz[0])

                # pool reads from all genome ids, each element is all genome_id at a time point
                for sample_index, fq_from_genome_ids in enumerate(longitudinal_sample_lists):
                    pooled_fq = os.path.join(output_dir, f"{self.ID}_all_pooled_t{sample_index}_mutated.fq.gz")
                    os.system(f"cat {' '.join(fq_from_genome_ids)} > {pooled_fq}")

            elif len(fq_gz_list_temp2[0]) == 2:
                ##############################
                """ whether pair end reads """
                ##############################
                longitudinal_sample_r1_lists = [[] for i in range(len(fq_gz_list_temp2))]
                longitudinal_sample_r2_lists = [[] for i in range(len(fq_gz_list_temp2))]
                for fq_gz_list_temp3 in list(genome_id_modified_fq_dict.values()):
                    for sample_index, fq_gz in enumerate(fq_gz_list_temp3):
                        longitudinal_sample_r1_lists[sample_index].append(fq_gz[0])
                        longitudinal_sample_r2_lists[sample_index].append(fq_gz[1])

                # pool reads from all genome ids, each element is all genome_id at a time point
                for sample_index, fq_from_genome_ids in enumerate(longitudinal_sample_r1_lists):
                    pooled_fq_r1 = os.path.join(output_dir, f"{self.ID}_all_pooled_t{sample_index}_mutated_R1.fq.gz")
                    os.system(f"cat {' '.join(fq_from_genome_ids)} > {pooled_fq_r1}")
                for sample_index, fq_from_genome_ids in enumerate(longitudinal_sample_r2_lists):
                    pooled_fq_r2 = os.path.join(output_dir, f"{self.ID}_all_pooled_t{sample_index}_mutated_R2.fq.gz")
                    os.system(f"cat {' '.join(fq_from_genome_ids)} > {pooled_fq_r2}")

            else:
                print(f"Error! The modified reads of {self.ID} are neither single end nor pair end!")

    def output_truth_of_mutation(self, mutation_info_dict, output_dir, geno_id):
        """
        :param mutation_info_dict: {pos: {"base":[ref, alt], "ref_count": [], "alt_count": []}}
        :param output_dir:
        :param geno_id: the genome id of species/strain, e.g. Geno1.0
        :return:
        """
        with open(os.path.join(output_dir, f"{self.ID}_{geno_id}_true_mutation.txt"), "w") as true_mutation_f:
            for pos in sorted(mutation_info_dict.keys()):
                mutation_info_at_pos = mutation_info_dict[pos]
                # pos + 1 because pysam is 0-based
                ref_count = [str(i) for i in mutation_info_at_pos["ref_count"]]
                alt_count = [str(i) for i in mutation_info_at_pos["alt_count"]]
                true_mutation_f.write(str(pos + 1) + "\t" + "\t".join(mutation_info_at_pos["base"]) + "\t" + ":".join(
                    ref_count) + "\t" + ":".join(alt_count) + "\n")
