import numpy as np
import stemsim.mapping as reads_mapping
import stemsim.reads_info_to_modify as ritm
import stemsim.modify_fastq as mf
import os


class Subject:
    """
    Subject has multiple microbiome samples (from CAMISIM) at several time points (or repeated measures)
    last update 11/14/2022
    """

    def __init__(self, model, subject_id, genome_id_list, genome_id_fas_dict, genome_id_read_dict, bowtie2_path, samtools_path, input_type, thread, my_logger):
        self.model = model                                              # "camisim" or "reads"
        self.ID = subject_id
        self.genome_id = genome_id_list                                 # [geno0, geno1, geno2]
        self.genome_id_fas_dict = genome_id_fas_dict                    # {"geno0": "/data/genomes/ref1.fas"}
        self.genome_id_bowtie_ref_dict = {genome_id: ".".join(fas.split(".")[:-1]) for genome_id, fas in
                                          genome_id_fas_dict.items()}   # {"geno0": "/data/genomes/ref1"}
        self.genome_id_read_dict = genome_id_read_dict                  # {"geno0": [[fq0_path], [fq1_path], [fq2_path]]} or unaligned bam
        self.genome_id_aligned_bam_dict = {}                            # {"geno0": [bam0_path, bam1_path, bam2_path]}
        self.thread = thread
        self.bowtie2_path = bowtie2_path
        self.samtools_path = samtools_path
        self.input_type = input_type                                    # "fq" or "bam"
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

    def generate_mutations(self, ref_base_prop_vec, mutation_traj_combination_dict, mutation_number,
                           base_mutation_freq_dict, output_dir):
        """
        Step2: Generate reads info that needs to be modified
        :param ref_base_prop_vec:
        @ like [1/6, 2/6, 1/6, 2/6], in the order of "A", "G", "C", "T"
        :param mutation_traj_combination_dict: {"combination1": {"prob": 0.3, "longitudinal_prop": [0.1, 0.7, 0.8]}}
        :param mutation_number: total number of mutation to generate
        :param base_mutation_freq_dict: {"A": {"alt_base_list": ["G", "C", "T"], "alt_freq_list": [0.25, 0.5, 0.25]}, ...}
        :return:
        """
        for genome_id in self.genome_id:
            self.logger.info(f"Generate mutations for reads from {genome_id} ... ...")
            input_bam_list = self.genome_id_aligned_bam_dict[genome_id]
            ref_fas_path = self.genome_id_fas_dict[genome_id]
            # record all created mutations, {"read_name1": {distance_to_end1: mutation_base}}, distance is 0-based
            #####################################
            """ mutation information in reads """
            #####################################
            # reads_mutations_record is list, length is the number of longitudinal samples
            # all_mutation_info is dict, {pos: {"base":[ref, alt], "ref_count": [], "alt_count": []}}
            reads_mutations_record, all_mutation_info = ritm.summarize_reads_to_modify(input_bam_list, ref_fas_path,
                                                                                       ref_base_prop_vec,
                                                                                       mutation_traj_combination_dict,
                                                                                       mutation_number,
                                                                                       base_mutation_freq_dict,
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

                elif len(processed_fq_path) == 2:
                    for fq_i, input_fq_i_path in enumerate(sorted(processed_fq_path)):
                        # fq_i is 0 or 1, represent R1 or R2 in paired fqs
                        output_fq_path = os.path.join(output_dir,
                                                      f"{self.ID}_{self.genome_id}_t{fq_index}_mutated_R{fq_i + 1}.fq")
                        mf.generate_modified_fq(input_fq_i_path, output_fq_path, reads_mutations_record[fq_index],
                                                whether_gzip=True)

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
