#####################################
## Step1: Simulation of raw sequencing reads ##
#####################################
### users can simulate microbial communities by CAMISIM or raw reads of a single strain by ART
### users can choose one of the following option (a or b)

### Step1a: simulate metagenomic data of complex microbial community by CAMISIM (need to be installed by the user)
## detailed instruction of CAMISIM can be found: https://github.com/CAMI-challenge/CAMISIM/wiki/User-manual#metagenome-simulation
python metagenomesimulation.py {PATH}/CAMISIM_example_config.ini
#-------------------------------------------------------------------------------------
### Step1b: simulate raw sequencing reads by ART (need to be installed by the user)
##  detailed instruction can be found: https://github.com/CAMI-challenge/CAMISIM/wiki/User-manual#metagenome-simulation
## simulate three longitudinal samples (t0, t1, t2) using the reference {genome1.0_fas}
art_illumina -ss HS20 -i {genome1.0_fas} -p -l 100 -f 20 -m 200 -s 10 -o genome1.0_ART_t0
art_illumina -ss HS20 -i {genome1.0_fas} -p -l 100 -f 20 -m 200 -s 10 -o genome1.0_ART_t1
art_illumina -ss HS20 -i {genome1.0_fas} -p -l 100 -f 20 -m 200 -s 10 -o genome1.0_ART_t2

#############################################
## Step2: Simulation of  short term evolution mutations ##
#############################################
### Step2a: follow Step1a; No line breaks are required between parameters
## Under this mode, stemsim takes the output directory of CAMISIM as input directory, and generate longitudinal mutations on the simulated raw sequencing data.   
## # the absolute path of genome_to_id.tsv file of CAMISIM, each line are a genome_id and a absolute path of fasta separated by Tab (see example below);
python {PATH_stemsim}/stemsim.py -m camisim -f {PATH_to}/config.txt -s test_rep0 --input_from_camisim {OUTPUT_directory_of_CAMISIM} --genome_to_id {PATH_to}/genome_to_id.tsv -o {OUTPUT_directory_stemsim} -l stemsim_test1.log

## Step2a optional parameters
# the absolute path of nonmutated_region_file_list.txt, in which each line stores the genome ID and its corresponding genomic-region file (the genome ID and genomic-region file are separated by Tab). Three formats (.gtf, .gff, and .bed) of genomic-region file are acceptable. For each genome listed, no mutation will be generated in the genomic regions listed in the genomic-region file;
--nonmutated_list {PATH_to}/example_nonmutated_region_file_list.txt 
# the absolute path of the directory of SnpEff, where SnpEff.jar is stored
-p {SnpEff_diretory}
#  absolute path of ID-matching table, in which each line contains a genome ID and its matching ID in snpEff database.
-a {PATH_matching_table}/example_genomeID_match_snpEff_ID.txt

#-------------------------------------------------------------------------------------
### Step2b: follow Step1b
python {PATH_stemsim}/stemsim.py -m reads -f STEMSIM_example_config.txt -s test_rep0 --input_reads genome1.0_ART_t0.fq:genome1.0_ART_t1.fq:genome1.0_ART_t2.fq -o {OUTPUT_directory_stemsim} --reference_genome {genome1.0_fas} --genome_id Genome.1.0 -l stemsim_test1.log

## Step2b optional parameters
# the absolute of genomic-region file (.gtf or .gff or .bed), in which no mutation will be generated for this genome
--nonmutated_file example_nonmutated_region_in_Bifidobacterium_breve_12L.gtf 	
# the absolute path of the directory of SnpEff, where SnpEff.jar is stored
-p {SnpEff_diretory}
#  absolute path of ID-matching table, in which each line contains a genome ID and its matching ID in snpEff database.
-a {PATH_matching_table}/genomeID_match_snpEff_ID.txt

