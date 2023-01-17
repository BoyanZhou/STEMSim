**STEMSim for the simulation of short-term evolutionary mutations**
===========================================

This code was used for the paper "STEMSIM: a simulator of short-term evolutionary mutations for longitudinal metagenomic data".


## Requirements

1.  Bowtie2 (tested with 0.16.0.1)  
http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
2.  samtools (tested with 1.9)  
http://www.htslib.org/doc/
3.  CAMISIM (tested with 1.3, optional, only required for "camisim" mode)  
https://github.com/CAMI-challenge/CAMISIM
4.  This code was written and tested on python 3.6.5, and requires the following packages:
    - numpy (tested with 1.19.5)
    - pysam (tested with 0.16.0.1)
    
    **Note:** CAMISIM is only necessary for "camisim" mode. But we recommend using the output of CAMISIM as input.
	      Because CAMISIM is a well-established pipeline for the simulation of microbial community.

## Install

1. Install the required software and packages

2. Download all files and deposit them under a directory.

## Usage

There are following steps for analysis of longitudinal metagenomic data using LongStrain.

### Mode1: camisim
Under this mode, stemsim takes the output directory of CAMISIM as input directory, and generate longitudinal mutations on the simulated raw sequencing data.   
```python /path_to_stemsim_directory/stemsim.py -m camisim -f PATH/config.txt -s test_subject1 --input_from_camisim PATH_to_output_directory_of_CAMISIM/ -o PATH_to_output_directory_of_STEMSim/ -l Path/stemsim_test1_log```

* ```-m camisim``` is the working mode for processing CAMISIM output;
* ```-f PATH/config.txt``` is the absolute path to config file (the template is provided in config.txt);
* ```-s test_subject1``` is the ID of current subject;
* ```--input_from_camisim PATH_to_output_directory_of_CAMISIM/``` is the path to the output of CAMISIM;
* ```--genome_to_id PATH_to_genome_to_id/genome_to_id.tsv``` is the absolute path of genome_to_id.tsv file of CAMISIM, each line are a genome_id and a absolute path of fasta separated by Tab;
* ```-o PATH_to_output_directory_of_STEMSim/``` is the output directory of STEMSim;
* ```-l Path/stemsim_test1_log``` is the absolute path of log file.

### Mode2: simulated raw reads
Under this mode, stemsim takes the simulated raw reads from other software (such as ART) as input, and generate longitudinal mutations on the simulated raw sequencing data.  
```python /path_to_stemsim_directory/stemsim.py -m reads -f PATH/config.txt -s test_subject1 --input_reads PATH/sample0_R1.fq,PATH/sample0_R2.fq:PATH/sample1.fq:PATH/sample2.fq -o PATH_to_output_directory_of_STEMSim/ --reference_genome PATH/species0_strain0.fas --genome_id species0_strain0 -l Path/stemsim_test1_log```

* ```-m reads``` is the working mode for processing raw sequencing data generated by other reads simulator;
* ```-f PATH/config.txt``` is the absolute path to config file (the template is provided in config.txt);
* ```-s test_subject1``` is the ID of current subject;
* ```--input_reads PATH/sample0_R1.fq,PATH/sample0_R2.fq:PATH/sample1.fq:PATH/sample2.fq``` is input raw read files from longitudinal samples separated by ':'; for paired reads from the same sample, they are separated by ',';
* ```-o PATH_to_output_directory_of_STEMSim/``` is the output directory of STEMSim;
* ```--reference_genome PATH/species0_strain0.fas``` is the absolute path of reference genome (.fas) from which raw reads are simulated;
* ```--genome_id species0_strain0``` is the genome id of reference genome;
* ```-l Path/stemsim_test1_log``` is the absolute of log file.

### Config file: the template is provided in config.txt

**Warnning:** STEMSim is designed to work on longitudinal or concurrent simulated raw sequencing reads of microbial communities. Although it can also works on real metagenomic data, it only works well on species with high depth and well-characterized genomes.


**Please contact boyanzhou1992@gmail.com for any questions or bug reporting.**