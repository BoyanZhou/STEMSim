import numpy as np
import stemsim.mapping as reads_mapping
import stemsim.reads_info_to_modify as ritm
import stemsim.modify_fastq as mf
import stemsim.ani_calculation
import os
import pysam
import stemsim.indel as indel
from collections import Counter


class Subject:
    """
    Subject has multiple microbiome samples (from CAMISIM) at several time points (or repeated measures)
    last update 11/14/2022
    """

    def __init__(self, model, subject_id, genome_id_list, genome_id_fas_dict, genome_id_read_dict,
                 genome_id_nonmutated_file_dict, bowtie2_path, samtools_path, input_type, fixed_mutation_number,
                 mutation_rate, substitution_model, thread, geno_id_snpeff_id_dict, snpeff_dir, my_logger):
        self.model = model                                              # "camisim" or "reads"
        self.ID = subject_id
        self.genome_id = genome_id_list                                 # [geno0, geno1, geno2]
        self.genome_id_fas_dict = genome_id_fas_dict                    # {"geno0": "/data/genomes/ref1.fas"}
        self.genome_id_bowtie_ref_dict = {genome_id: ".".join(fas.split(".")[:-1]) for genome_id, fas in
                                          genome_id_fas_dict.items()}   # {"geno0": "/data/genomes/ref1"}
        self.genome_id_read_dict = genome_id_read_dict                  # {"geno0": [[fq0_path], [fq1_path], [fq2_path]]} or unaligned bam
        self.genome_id_nonmutated_file_dict = genome_id_nonmutated_file_dict        # {"geno0": "PATH\nonmutated_file.gtf"}
        self.genome_id_chr_dict = {}                                    # {"geno0": ["chr1", "chr2"]}

        self.genome_id_aligned_bam_dict = {}                            # {"geno0": [bam0_path, bam1_path, bam2_path]}
        self.genome_id_sample_name_dict = {}                            # {"geno0": [sample0_name, ...]}
        self.genome_id_base_prop = {}                                   # {"geno0": {"A": 1/4, "G": 1/4, "C": 1/4, "T": 1/4}}
        self.thread = thread
        self.bowtie2_path = bowtie2_path
        self.samtools_path = samtools_path
        self.input_type = input_type                                    # "fq" or "bam"

        self.fixed_mutation_number = fixed_mutation_number
        self.mutation_rate = mutation_rate
        self.substitution_model = substitution_model    # must be JC69, K80, HKY85, TN93, or REV

        self.insertion_alpha = 0
        self.insertion_beta = 0
        self.whether_fixed_insertion_number = True
        self.fixed_insertion_number = 0
        self.insertion_mutation_rate = 0

        self.deletion_alpha = 0
        self.deletion_beta = 0
        self.whether_fixed_deletion_number = True
        self.fixed_deletion_number = 0
        self.deletion_mutation_rate = 0

        self.long_insertion_mu = 0          # mean of normal distribution
        self.long_insertion_sigma = 0       # sd of normal distribution
        self.long_insertion_number = 0

        self.long_deletion_mu = 0          # mean of normal distribution
        self.long_deletion_sigma = 0       # sd of normal distribution
        self.long_deletion_number = 0

        self.genome_id_snpeff_id_dict = geno_id_snpeff_id_dict
        self.snpeff_dir = snpeff_dir
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

    def add_insertion_parameter(self, alpha, beta, fix_num, mutation_rate, whether_use_mutation_rate):
        self.insertion_alpha = alpha
        self.insertion_beta = beta
        self.whether_fixed_insertion_number = not whether_use_mutation_rate
        self.fixed_insertion_number = fix_num
        self.insertion_mutation_rate = mutation_rate

    def add_deletion_parameter(self, alpha, beta, fix_num, mutation_rate, whether_use_mutation_rate):
        self.deletion_alpha = alpha
        self.deletion_beta = beta
        self.whether_fixed_deletion_number = not whether_use_mutation_rate
        self.fixed_deletion_number = fix_num
        self.deletion_mutation_rate = mutation_rate

    def add_long_insertion_parameter(self, mu, sigma, number):
        self.long_insertion_mu = mu
        self.long_insertion_sigma = sigma
        self.long_insertion_number = number

    def add_long_deletion_parameter(self, mu, sigma, number):
        self.long_deletion_mu = mu
        self.long_deletion_sigma = sigma
        self.long_deletion_number = number

    def generate_insertion(self, num, genome_id):
        # array of bases of insertion
        insertion_length_array = indel.generate_indel_lens(self.insertion_alpha, self.insertion_beta, num)
        insertion_array = indel.generate_insertion_array(insertion_length_array, self.genome_id_base_prop[genome_id])
        return insertion_array

    def generate_deletion(self, num, genome_id):
        # array of lengths of insertion
        insertion_length_array = indel.generate_indel_lens(self.insertion_alpha, self.insertion_beta, num)
        insertion_array = indel.generate_insertion_array(insertion_length_array, self.genome_id_base_prop[genome_id])
        return insertion_array

    @staticmethod
    def merge_and_sort_bin(bin_list):
        # bin_list = [[0, 100], [500, 1000], [3000, 3500], [150, 750]]
        # [[start0, end0], [start1, end1], [start2, end2]]
        bin_list_sorted = sorted(bin_list, key=lambda x: x[0])
        merge_sorted_bin_list = []
        current_bin = []
        for bin_i in bin_list_sorted:
            if len(current_bin) == 0:
                current_bin = bin_i
            else:
                # current_bin = [0, 100]
                if bin_i[0] <= current_bin[1]:
                    current_bin[1] = max(current_bin[1], bin_i[1])
                else:
                    merge_sorted_bin_list.append(current_bin)
                    current_bin = bin_i
        if len(current_bin) > 0:
            merge_sorted_bin_list.append(current_bin)
        return merge_sorted_bin_list

    def get_nonmutated_region(self, geno_id):
        """
        Get the dict of nonmutated_region of each chromosome from stored gtf/gff/bed path
        :param geno_id:
        :return: {"chr1": [[0, 100], [200, 300]], "chr2": ...}
        """
        # self.genome_id_nonmutated_file_dict = genome_id_nonmutated_file_dict  # {"geno0": ""PATH\nonmutated_file.gtf}
        if geno_id not in self.genome_id_nonmutated_file_dict:
            return {}
        abs_path_region_file = self.genome_id_nonmutated_file_dict[geno_id]
        nonmutated_region_dict = {}     # {chr_name: [[bin_start, bin_end]]}, 0-based
        if abs_path_region_file.endswith("gff") or abs_path_region_file.endswith("gtf"):
            # col0: chr name, col3: start (1-based), col4: end (1-based)
            with open(abs_path_region_file, "r") as region_f:
                for line in region_f:
                    if line.startswith("#"):
                        continue
                    cols = line.strip().split("\t")
                    if cols[0] not in nonmutated_region_dict:
                        nonmutated_region_dict.update({cols[0]: [[int(cols[3]) - 1, int(cols[4]) - 1]]})
                    else:
                        nonmutated_region_dict[cols[0]].append([int(cols[3]) - 1, int(cols[4]) - 1])

        elif abs_path_region_file.endswith("bed"):
            # col0: chr name, col1: start (0-based), col2: end (0-based)
            with open(abs_path_region_file, "r") as region_f:
                for line in region_f:
                    if line.startswith("#"):
                        continue
                    cols = line.strip().split("\t")
                    if cols[0] not in nonmutated_region_dict:
                        nonmutated_region_dict.update({cols[0]: [[int(cols[1]), int(cols[2])]]})
                    else:
                        nonmutated_region_dict[cols[0]].append([int(cols[1]), int(cols[2])])
        else:
            print(f"Warning! For genome {geno_id}, the format of nonmutated file is not accepted! "
                  f"Must have the suffix of bed or gtf or gff!")
            self.logger.info(f"Warning! For genome {geno_id}, the format of nonmutated file is not accepted! "
                             f"Must have the suffix of bed or gtf or gff!")
        for key, value in nonmutated_region_dict.items():
            # sort and merge
            nonmutated_region_dict[key] = self.merge_and_sort_bin(value)
        return nonmutated_region_dict

    def align_simulated_reads(self):
        """
        Step1: aligns the simulated reads (fq or unaligned bam) to corresponding original fas, get aligned bam list
        :return:
        """
        for genome_id in self.genome_id:
            self.logger.info(f"Processing the simulated reads from genome {genome_id} ... ...")
            bowtie_ref = self.genome_id_bowtie_ref_dict[genome_id]
            aligned_bam_path_list = []
            # sample_name_list = []
            ########################
            # alignment processing #
            ########################
            for unaligned_reads_abs_path_list in self.genome_id_read_dict[genome_id]:
                # unaligned_reads_abs_path_list = [fq0_path] or [fq0_0_path, fq0_1_path]
                aligned_bam_path = reads_mapping.map_simulated_reads(unaligned_reads_abs_path_list, bowtie_ref,
                                                                     self.thread, self.bowtie2_path, self.samtools_path,
                                                                     self.input_type, self.logger)
                aligned_bam_path_list.append(aligned_bam_path)
                # sample_name_list.append(sample_name)

            # store the path of aligned bam, and sample names
            self.genome_id_aligned_bam_dict.update({genome_id: aligned_bam_path_list})
            # self.genome_id_sample_name_dict.update({genome_id: sample_name_list})

    def generate_mutations(self, mutation_traj_combination_dict, whether_use_mutation_rate_number, substitution_model,
                           original_q_matrix, output_dir, working_model, pool_reads=True):
        """
        Step2: Generate reads info that needs to be modified
        @ like [1/6, 2/6, 1/6, 2/6], in the order of "A", "C", "G", "T"
        :param mutation_traj_combination_dict: {"combination1": {"prob": 0.3, "longitudinal_prop": [0.1, 0.7, 0.8]}}
        :param whether_use_mutation_rate_number: total number of mutation to generate
        :param substitution_model: JC69, K80, HKY85, TN93, or REV
        :param original_q_matrix: {"A": {"alt_base_list": ["C", "G", "T"], "alt_freq_list": [0.25, 0.5, 0.25]}, ...}
        :param output_dir: output_dir
        :param working_model: camisim or reads
        :param pool_reads: whether combine the generated fq.gz
        :return:
        """
        genome_id_modified_fq_dict = {}     # {"geno_id1": [, , ,]}
        for genome_id in self.genome_id:
            self.logger.info(f"Generate mutations for reads from {genome_id} ... ...")
            input_bam_list = self.genome_id_aligned_bam_dict[genome_id]
            ref_fas_path = self.genome_id_fas_dict[genome_id]
            sample_name_list = [f"{genome_id}_t{i}" for i in range(len(input_bam_list))]
            self.genome_id_sample_name_dict.update({genome_id: sample_name_list})
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

            ##################################
            # assign mutations to chromosome #
            ##################################
            # according to the length of chromosomes/contigs
            pysam_ref = pysam.FastaFile(ref_fas_path)
            self.genome_id_chr_dict.update({genome_id: pysam_ref.references})
            total_length_of_chrs = sum(pysam_ref.lengths)
            chr_prop_list = [i/total_length_of_chrs for i in pysam_ref.lengths]   # [1/6, 2/6, 1/6, 2/6], proportions of four chromosomes according to lengths
            # # 1. dict of mutation number assigned to each chromosome by fixed number or relative rate
            if whether_use_mutation_rate_number:
                mutation_number_dict = {i: int(j * self.mutation_rate) for i, j in
                                        zip(pysam_ref.references, pysam_ref.lengths)}
            else:
                mutation_number_dict = {i: j for i, j in zip(pysam_ref.references,
                                                             np.random.multinomial(n=self.fixed_mutation_number,
                                                                                   pvals=chr_prop_list))}

            # # 2. dict of insertion number assigned to each chromosome by fixed number or relative rate
            if self.whether_fixed_insertion_number:
                insertion_dict = {i: self.generate_insertion(j, genome_id) for i, j in zip(pysam_ref.references, np.random.multinomial(n=self.fixed_insertion_number, pvals=chr_prop_list))}
            else:
                insertion_dict = {i: self.generate_insertion(int(j*self.insertion_mutation_rate), genome_id) for i, j in zip(pysam_ref.references, pysam_ref.lengths)}

            # # 3. dict of deletion number assigned to each chromosome by fixed number or relative rate
            if self.whether_fixed_deletion_number:
                deletion_length_dict = {i: indel.generate_indel_lens(self.deletion_alpha, self.deletion_beta, j) for i, j in zip(pysam_ref.references, np.random.multinomial(n=self.fixed_deletion_number, pvals=chr_prop_list))}
            else:
                deletion_length_dict = {i: indel.generate_indel_lens(self.deletion_alpha, self.deletion_beta, int(j*self.deletion_mutation_rate)) for i, j in zip(pysam_ref.references, pysam_ref.lengths)}

            # # 4. dict of nonmutated_region {"chr1":[[200, 500], [700, 1000]]}, no overlap, from small to large
            nonmutated_region_dict = self.get_nonmutated_region(genome_id)      # can be empty {}

            # # 5. dict of long insertion assigned to each chromosome
            long_insertion_dict = {}  # {"chr1": [[start0, end0], []]}
            long_deletion_dict = {}
            if self.long_insertion_number > 0 or self.long_deletion_number > 0:
                long_insertion_number_each_chr = np.random.multinomial(n=self.long_insertion_number, pvals=chr_prop_list)   # list
                long_deletion_number_each_chr = np.random.multinomial(n=self.long_deletion_number, pvals=chr_prop_list)     # list
                for chr_name, chr_len, chr_long_insertion_num, \
                    chr_long_deletion_num in zip(pysam_ref.references, pysam_ref.lengths,
                                                 long_insertion_number_each_chr,long_deletion_number_each_chr):
                    long_insertion_length_array = np.random.normal(loc=self.long_insertion_mu,
                                                                   scale=self.long_insertion_sigma,
                                                                   size=chr_long_insertion_num).astype(int)
                    long_deletion_length_array = np.random.normal(loc=self.long_deletion_mu,
                                                                  scale=self.long_deletion_sigma,
                                                                  size=chr_long_deletion_num).astype(int)
                    long_indel_length_list = list(long_insertion_length_array) + list(long_deletion_length_array)

                    # ## get long indel range for this chromosome ###
                    long_insertion_dict.update({chr_name: []})
                    long_deletion_dict.update({chr_name: []})
                    # if long_indel exist in this chr
                    if len(long_indel_length_list) > 0:
                        avoided_region_list = []
                        if chr_name in nonmutated_region_dict:
                            avoided_region_list = nonmutated_region_dict[chr_name]
                        # combined list of start pos of insertion and deletion
                        indel_start_pos_list = indel.random_long_indel_from_blocks(chr_len, avoided_region_list, long_indel_length_list)
                        # 5.1 record insertion
                        insertion_range_list = []
                        for tmp_index, indel_length in enumerate(long_insertion_length_array):
                            start_pos = indel_start_pos_list[tmp_index]
                            if start_pos >= 0:
                                insertion_range_list.append([start_pos, start_pos+indel_length])
                        long_insertion_dict[chr_name] = sorted(insertion_range_list, key=lambda x: x[0])

                        # 5.2 record deletion
                        deletion_range_list = []
                        for tmp_index, indel_length in enumerate(long_deletion_length_array):
                            tmp_index += len(long_insertion_length_array)
                            start_pos = indel_start_pos_list[tmp_index]
                            if start_pos >= 0:
                                deletion_range_list.append([start_pos, start_pos + indel_length])
                        long_deletion_dict[chr_name] = sorted(deletion_range_list, key=lambda x: x[0])

            #####################################
            """ mutation information in reads """
            #####################################
            # 1. reads_mutations_record is list of dict, length is the number of longitudinal samples, e.g.
            # {"read_name1": {"range_list":[[distance_to_end_start, distance_to_end_end]], "alt_list": [mutation_base]}}
            # 2. list of reads set
            # 3. all_mutation_info is dict, {pos: {"base":[ref, alt], "ref_count": [], "alt_count": []}}
            reads_mutations_record, reads_to_remove_record, all_mutation_info, long_indel_info = ritm.summarize_reads_to_modify(
                input_bam_list, ref_fas_path, ref_base_prop_vec, mutation_traj_combination_dict, nonmutated_region_dict,
                long_insertion_dict, long_deletion_dict, insertion_dict, deletion_length_dict, mutation_number_dict,
                original_q_matrix, depth_threshold=5, bin_len=10000, space_between_mutations=150)

            """
            Output true mutations in VCF, and ANI, long insertion/deletion
            """
            ani_output_path = os.path.join(output_dir, f"{genome_id}_pairwise_ANI_between_longitudinal_samples.txt")
            stemsim.ani_calculation.calculate_ani(all_mutation_info, pysam_ref.references, pysam_ref.lengths,
                                                  self.genome_id_sample_name_dict[genome_id], ani_output_path)
            # Output true mutations in VCF
            self.output_truth_of_mutation(all_mutation_info, output_dir, genome_id, pysam_ref.references,
                                          self.genome_id_snpeff_id_dict, self.snpeff_dir)
            # if long indel need to be reported
            if (len(long_insertion_dict) + len(long_deletion_dict)) > 0:
                self.output_truth_of_long_indel(long_indel_info, output_dir, genome_id, pysam_ref.references)

            # ----------------------------------------------------------------------------------------------------------
            ################################
            """ Step3: modified fq files """
            ################################
            genome_id_modified_fq_dict.update({genome_id: []})      # record the path of modified fqs
            input_fq_path_list = self.genome_id_read_dict[genome_id]    # [[], [], []]
            for fq_index in range(len(input_fq_path_list)):
                # for each longitudinal sample
                processed_fq_path = []      # [XXX_1.fq, XXX_2.fq] or [XXX.fq]
                """ Pre-process fq according to different types """
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
                # 1. for single end reads
                if len(processed_fq_path) == 1:
                    output_fq_path = os.path.join(output_dir, f"{self.ID}_{genome_id}_t{fq_index}_mutated.fq")
                    mf.generate_modified_fq(processed_fq_path[0], output_fq_path, reads_mutations_record[fq_index],
                                            reads_to_remove_record[fq_index], whether_gzip=True)
                    genome_id_modified_fq_dict[genome_id].append([f"{output_fq_path}.gz"])
                # 2. for paired end reads
                elif len(processed_fq_path) == 2:
                    output_fq1_path = os.path.join(output_dir, f"{self.ID}_{genome_id}_t{fq_index}_mutated_R1.fq")
                    output_fq2_path = os.path.join(output_dir, f"{self.ID}_{genome_id}_t{fq_index}_mutated_R2.fq")
                    mf.generate_modified_paired_fqs(processed_fq_path[0], processed_fq_path[1], output_fq1_path, output_fq2_path,
                                                    reads_mutations_record[fq_index],
                                                    reads_to_remove_record[fq_index], whether_gzip=True)
                    genome_id_modified_fq_dict[genome_id].append([f"{output_fq1_path}.gz", f"{output_fq2_path}.gz"])
        # --------------------------------------------------------------------------------------------------------------
        # whether pool the modified reads from each genome/species/strain
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

    def output_truth_of_mutation(self, mutation_info_dict, output_dir, geno_id, chr_name_list, geno_id_snpeff_id_dict, snpeff_dir=""):
        """
        :param mutation_info_dict: {chr_name: {{position: {"type": "deletion", "base": [ref_base, mutated_base],
                                    "ref_count": [0] * len(input_bam_list), "alt_count": [0] * len(input_bam_list)}}}}
        :param output_dir:
        :param geno_id: the genome id of species/strain, e.g. Geno1.0
        :param chr_name_list: chr name of this genome
        :param geno_id_snpeff_id_dict: the id in SnpEff that matches with the genome id
        :param snpeff_dir: the directory of SnpEff
        :return:

        """
        ##############
        # output vcf #
        ##############
        output_vcf_path = os.path.join(output_dir, f"{self.ID}_{geno_id}_true_mutation.vcf")
        with open(output_vcf_path, "w") as true_mutation_f:
            # write header
            sample_ids = "\t".join(self.genome_id_sample_name_dict[geno_id])
            true_mutation_f.write(f"##fileformat=VCFv4.1\n"
                                  f"##source=STEMSim\n"
                                  f"##reference={self.genome_id_fas_dict[geno_id]}\n"
                                  f'##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">\n'
                                  f'##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
                                  f'##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">\n'
                                  f'##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">\n'
                                  f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_ids}\n")
            for chr_name in chr_name_list:
                if chr_name in mutation_info_dict:
                    mutation_info_in_chr_dict = mutation_info_dict[chr_name]
                    for pos in sorted(mutation_info_in_chr_dict.keys()):
                        # {position: {"type": "deletion", "base": [ref_base, mutated_base],
                        # "ref_count": [0] * len(input_bam_list), "alt_count": [0] * len(input_bam_list)}}
                        mutation_info_in_chr_at_pos_dict = mutation_info_in_chr_dict[pos]
                        ref_base, mutated_base = mutation_info_in_chr_at_pos_dict["base"]
                        ref_count_list = mutation_info_in_chr_at_pos_dict["ref_count"]
                        alt_count_list = mutation_info_in_chr_at_pos_dict["alt_count"]
                        total_dp = sum(ref_count_list) + sum(alt_count_list)
                        INFO = f"DP={total_dp}"             # total depth of all samples
                        FORMAT = "GT:DP:AD"     # DP for this sample
                        sample_info_list = []
                        # get GT
                        for ref_count, alt_count in zip(ref_count_list, alt_count_list):
                            if ref_count == 0:
                                if alt_count == 0:
                                    genotype = "./."
                                else:
                                    genotype = "1/1"
                            else:
                                if alt_count == 0:
                                    genotype = "0/0"
                                else:
                                    genotype = "0/1"
                            sample_info = f"{genotype}:{ref_count + alt_count}:{ref_count},{alt_count}"
                            sample_info_list.append(sample_info)
                        sample_info_list = "\t".join(sample_info_list)
                        # pos + 1 because pysam is 0-based
                        vcf_line = f"{chr_name}\t{pos + 1}\t.\t{ref_base}\t{mutated_base}\t100\t.\t{INFO}\t{FORMAT}\t" \
                                   f"{sample_info_list}\n"
                        true_mutation_f.write(vcf_line)

        ########################
        # annotation by SnpEff #
        ########################
        if len(snpeff_dir) > 0 and geno_id in geno_id_snpeff_id_dict:
            # if annotation match profile
            snpEff_path = os.path.join(snpeff_dir, "snpEff.jar")
            snpEff_config_path = os.path.join(snpeff_dir, "snpEff.config")
            snpeff_id = geno_id_snpeff_id_dict[geno_id]
            output_annotated_vcf_path = os.path.join(output_dir, f"{self.ID}_{geno_id}_annotated_true_mutation.vcf")
            snpeff_command = f"java -Xmx4g -jar {snpEff_path} -c {snpEff_config_path} -v {snpeff_id} {output_vcf_path} > {output_annotated_vcf_path}"
            self.logger.info(f"Annotation by snpEff:\n{snpeff_command}\n")
            os.system(snpeff_command)

    def output_truth_of_long_indel(self, long_indel_info_dict, output_dir, genome_id, chr_name_list):
        """

        :param long_indel_info_dict: {"chr1": {pos: {"range": [,], "type": "long_insertion", "prop": [0.1, 0.2, 0.5, 0.9]}}}
        :param output_dir:
        :param genome_id:
        :param chr_name_list:
        :return:
        """
        with open(os.path.join(output_dir, f"{self.ID}_{genome_id}_long_indel_record.txt"), "w") as long_indel_f:
            # write header
            long_indel_f.write(f"#CHROM\tStart\tEnd\tType\tLongitudinal_proportion\n")
            for chr_name in chr_name_list:
                if chr_name in long_indel_info_dict:
                    long_indel_info_in_chr_dict = long_indel_info_dict[chr_name]
                    for pos in sorted(long_indel_info_in_chr_dict.keys()):
                        indel_range = long_indel_info_in_chr_dict[pos]["range"]
                        longitudinal_prop = [str(i) for i in long_indel_info_in_chr_dict[pos]["prop"]]
                        longitudinal_prop_str = ",".join(longitudinal_prop)
                        long_indel_f.write(f"{chr_name}\t{indel_range[0]}\t{indel_range[1]}\t"
                                           f"{long_indel_info_in_chr_dict[pos]['type']}\t"
                                           f"{longitudinal_prop_str}\n")
