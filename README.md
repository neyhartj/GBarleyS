# GBarleyS
###### This software package is a compilation of scripts to run a pipeline designed to call variants (SNPs) from raw Genotyping-by-Sequencing (GBS) data. The scripts rely heavily on batch job submission, paralleled processing, and high-performance computation. Much of the workflow/content/ideas were inspired from `sequence_handling` from [MorrellLAB/sequence_handling](https://github.com/MorrellLAB/sequence_handling). The package is optimized for the genomic resources available for barley, but could be adjusted for any GBS data analysis involving reference-based variant discovery.
___

## Pipeline Steps
The GBarleyS pipeline is divided into three parts:

#### Raw Data Handling and Quality Control
This is the initial stage and involves demultiplexing, quality assessment, and quality trimming of sequence reads.
#### Alignment and Variant Calling
This is the second stage and involves aligning reads to the reference genome, and using the Genome Analysis Toolkit (GATK) to call variants.
#### Post-Variant Processing
This is the last stage and involves analysis and conversion of variants to other file formats for subsequent analysis.

The following document will guide you through using the pipeline and describes the order in implementing scripts.
___

## Dependencies
The `GBarleyS` pipeline relies on a number of other software tools to run correctly. These tools are all publicly available by clicking on the links.

Software | Version
-------- | --------
[Bcftools](http://www.htslib.org/download/) | 1.2
[Bedtools](http://bedtools.readthedocs.org/en/latest/index.html) | 2.17.0
[Bowtie2](http://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.4/) | 2.2.4
[Cutadapt](https://cutadapt.readthedocs.org/en/stable/) | 1.8.1
[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) | 0.11.2
[FastX_Toolkit](http://hannonlab.cshl.edu/fastx_toolkit/download.html) | 0.0.14
[Genome Analysis Toolkit (GATK)](https://www.broadinstitute.org/gatk/download/) | 3.4-46 
[GNU Parallel](http://www.gnu.org/software/parallel/) | NA
[SAMtools](http://www.htslib.org/download/) | 0.1.18
[VCFtools](https://vcftools.github.io/index.html) | 0.1.14

> **NOTE:** The GATK requires you to sign up for a license, but this is free for academic use.

___

## GBarleyS Installation
Before using this software package, you must install it first! To clone the `GBarleyS` repository, use the following commands:
```
git clone https://github.umn.edu/neyha001/GBarleyS.git
```
This will copy the entire repository in the current working directory. It is recommended that the repository not be copied into a directory in which you plan to run the pipeline. The first step is to execute the `run_installer.sh` script as such:
```
cd GBarleyS/Pipeline_Scripts
sh run_installer.sh /path/to/directory
```
This will create a number of sub-directories in the the specified directory. The script also copies all of the `GBarleyS` pipeline scripts into the specified directory. Last, the script checks to see if the required software is present. See **Dependencies** for a list of required software.
___

## A General Primer on How the GBarleyS Scripts Operate
The bulk of the GBarleyS pipeline is written in shell in a way that calls upon the programs listed above. Each script relies on the Portable Batch System (PBS) for job scheduling. Initially designed to work on the Minnesota Supercomputing Institute (MSI), further development should make the pipeline more widely useable.

#### Editing Scripts
To run the pipeline, a user opens up the intended script for editing with a text editor, exemplified here with VIM. The user then provides information for the requested variables. For example, let's say we want to run the `run_fastx_barcode_splitter.sh` script. First we would open the script for editing:
```
vim run_fastx_barcode_splitter.sh
```
The first few lines of the script give a description of its function, followed by version information for the software that is called by the script (in this case, the FastX_Toolkit):
```
#!/bin/bash

# Description
# This script demultiplexes a number of multiplexed fastq files based on the barcode sequences that 
## appears at the beginning of the read. This script calls on the FastX_Barcode_Splitter, 
## which only allows barcodes of a constant length. Therefore, this script will run n jobs where n 
## is the number of different barcode lengths.

# Version info
## FastX_Toolkit: 0.0.14
```

#### User Input
Next, there is a section where the user provides information for different variables. This should be the only section in which the user makes edits. Of course, one is free to make additional edits, but at one's own risk. Here are the first four variables that require user input for this script:
```
# The directory where the multiplexed samples are located or where their symbolic links are located
INDIR=

# The directory to output the demultiplexed fastq files
OUTDIR=

# The key file with columns: Flowcell Lane Barcode Sample PlateName Row Column
KEYFILE=

# The name of the current project (e.g. 2row_TP)
PROJECT=''
```
Most scripts in the GBarleyS pipeline require the `INDIR`, `OUTDIR`, and `PROJECT` variables. In this case, the `KEYFILE` variable is specific to this script. **The important thing to note is that all of the GBarleyS scripts will have a section at the beginning for the user to provide inputs to variables.**

#### Software Settings
In many scripts, one of the variables to be defined is a string of settings for the software called by the script. For instance, in the `run_fastx_barcode_splitter.sh` script, the settings for the FastX_Toolkit are already defined:
```
FASTX_SETTINGS='--bol --mismatches 1'
```
These sort of settings are defined to optimize the pipeline, however they are flexible enough to allow for user modifications. In most cases, using the default settings is acceptable.

#### Running the Scripts
All GBarleyS pipeline scripts are executed using shell. Once a user is done editing a script and has verified the correctness of their input, the user may run the script with the following command (again we are using the `run_fastx_barcode_splitter.sh` script as an example):
```
sh run_fastx_barcode_splitter.sh
```
The script will perform its operations, call upon software, and submit a computing job via PBS. The user will receive a notification of success such as the following:
```
Jobs away!
```
Additionally, a simple text log file is created that includes the date, username, name of the script, and inputs from the user. Once the job is submitted, the user can move on to other activities (such as getting coffee or playing raquetball) while the job runs!

> **NOTE:** When defining variables in shell, there should be no space between the name of the variable, the '=' sign, and the input for said variable. For instance, `KEYFILE=keyfile.txt` is acceptable, but `KEYFILE = keyfile.txt` is not. This may be new for users familiar with other languages (e.g. Python, Perl, R).
___

## Raw Data Handling and Quality Control
Before any data processing can occur, it is imperative to assess the quality of sequence data that is present. This will allow you to determine if technical errors occurred in sequencing (i.e. bubbles in the flowcell) or if there is bias or batch effect. 

##### 1. `run_fastx_barcode_splitter.sh`
This script uses the `fastx_barcode_splitter.pl` from the FastX Toolkit to separate a multiplexed `.fastq` or `.fastq.gz` file (typical of GBS raw data) to a single `.fastq` file for each sample. The number of resulting files depends on the plex level, or, the number of sample assessed on a single lane on a flowcell.

##### 2. `run_post_demultiplex_cleanup.sh`
This script runs some cleanup operations after the demutiplexing step. The `fastx_barcode_splitter.pl` tool results in uncompressed `.fastq` files, which can take up substatial disk space. The script re-compresses the `.fastq` files using `gzip` to `.fastq.gz` files.

##### 3. `run_quality_control.sh`
This script uses both the `FastX Toolkit` and `Cutadapt` to make quality adjustments to the raw read data. This includes trimming reads based on Phred quality score (scale of base call likelihood), trimming reads based on length, removal of barcode sequences, and removal of adapter sequences. The script also outputs logs for each program which detail the number and proportion of reads/bases that were trimmed.

##### `run_assess_quality.sh`
Assessing the quality of your data is critically important to a successful workflow. Therefore, it is encouraged that the user runs the `run_assess_quality.sh` script **both before and after** performing quality control. This script uses the `FastQC` program to assess the quality of all of the demultiplexed `.fastq.gz` files. For each file, an individual FastQC report is generated. The `.html` version of these reports can be visualized using a browser. Additionally, raw data is saved as a `.zip` file. If you want to get an idea of the read quality across all of your samples (say from multiple lanes or multiple flowcells), you can specify the `STATS` variable as 'true.' This will call upon other scripts ( `fastqc_stats_parser.py` and `pipeline_graphing_functions.R`) to generate a `.pdf` file with graphs of some important read qualities parameters.
___

## Alignment and Variant Calling
After raw reads have been assessed and trimmed for quality, the processing steps can begin. The next stage involves aligning reads to the reference genome and calling variants using GATK. This is by far the most computationally-intensive stage, and, depending on the number of samples, will likely required the use of high-performance computing resources.

### Read Mapping

##### 1. `run_read_mapping.sh`
This script calls on the `Bowtie2` or `BWA` programs to align the reads to the reference genome. The standard output from `Bowtie2` and `BWA` is a `.sam` file, however this script uses `SAMtools` to conver the `.sam` file into a binary `.bam` file.

##### 2. `run_samtools_processing.sh`
This script does some editing and processing of the `.bam` files generated by the `run_read_mapping.sh` script. First, stats are generated on the raw `.bam` files. Next, the files are filtered to remove reads that did not align and reads that aligned more than once. Additionally, the files are sorted on tile position and stats are generated for the processed files. Finally, the `.bam` files are indexed for use in the GATK pipeline.

### GATK Pipeline
The GATK pipeline includes several steps to process the `.bam` files, call variants, and process the variants. The GATK pipeline uses a known set of variants to recalibrate base and variant quality scores. If you do not have a known set of SNP and indel variants, you have the option of implementing the pipeline without this prior set of data, imposing strict filters, and then using these high-confidence SNPs as the known set of variants.

___
#### Follow these steps if you do not have a set of known variants
##### 1. `run_indel_realignment.sh`
This script runs the [IndelRealigner](https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_indels_IndelRealigner.php) program in the GATK pipeline. Briefly, it performs local realignment of reads around indels, since alignment mismatches in these regions can result in erroneous SNP calls. The script outputs `.bam` files with realignment completed.

##### 2. `run_haplotypecaller.sh`
This script runs the [HaplotypeCaller](https://github.umn.edu/neyha001/GBarleyS/blob/master/Pipeline_Scripts/run_haplotypecaller.sh) program in the GATK pipeline. Briefly, it identifies active regions (regions that may contain variants) and performs de-novo reassembly of the active region, detecting possible haplotypes. The program then performs pairwise alignment of reads against the possible haplotypes in the active region to determine the likelihood of all possible alleles at a given site. The program uses the likelihood of alleles given the read data to calculate the likelihood of each genotype per sample given the read data per sample. The most likely genotype is assigned to the sample. This script uses the gVCF option of the HaploytypeCaller, which is useful for cohort genotyping.

##### 3. `run_combine_gvcfs.sh`
This script runs the [CombineGVCFs](https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_variantutils_GenotypeGVCFs.php) program to merge gVCF files for each sample into hierarchical cohorts of `.g.vcf` files. The script impements a two-step process: First, `.g.vcf` files are combined by each flowcell-lane combination (e.g. C613JACXX_1) into 4 sets (for 96-plex, this means that 4 sets of 24 `.g.vcf` files are combined). Second, each set of `.g.vcf` files is combined into a single `.g.vcf` file per flowcell-lane combination. For instance, say you have one flowcell with 8 lanes (therefore 8 flowcell-lane combinations) and each lane contains 96 samples. From the `run_haplotypecaller.sh` script, you should have 8 * 96 = 768 `.g.vcf` files. The first part of `run_combine_gvcf.sh` will combine the 768 files into 4 sets of 24 files for each of 8 flowcell-lane combinations, therefore creating 4 * 8 = 32 `.g.vcf` files. The second part combines those sets into one `cohort.g.vcf` file per flowcell-lane combination. Therefore in this instance, 8 `cohort.g.vcf` files will be created. Since the intermediate set `g.vcf` files will not be needed, one might consider deleting them to save disk space. The `cohort.g.vcf` files, however, are very important and will be needed for additional genotyping as the population increases.

##### 4. `run_genotype_gvcf.sh`
This script runs the [GenotypeGVCFs](https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_variantutils_GenotypeGVCFs.php) program in the GATK pipeline. It aggregates all of the genotype likelihoods in each `.g.vcf` file and assigns genotypes at each site. The output is a single `.vcf` file of raw variants ready for filtering.

##### 5. `filter_variants.sh`
This script uses VCFtools to implement hard-filters on the VCF file. Filtering parameters may be adjusted as needed. This step is intended to create a `.vcf` file of high-quality SNP variants which may be used as a set of known variants when running the GATK pipeline using recalibration.
___
#### Continue with these steps if you have a set of known variants
##### 1. `run_indel_realignment.sh`
This script runs the [IndelRealigner](https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_indels_IndelRealigner.php) program in the GATK pipeline. Briefly, it performs local realignment of reads around indels, since alignment mismatches in these regions can result in erroneous SNP calls. The script outputs `.bam` files with realignment completed. When running the GATK pipeline using recalibration, this program requires `.vcf` files of known indels. 

##### 2. `run_baserecalibrator.sh`
This script runs the [BaseRecalibrator](https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_bqsr_BaseRecalibrator.php) program in the GATK pipeline on the `.bam` files from the `run_indel_realignment.sh` step. Briefly, the program uses a set of known SNPs to readjust the base quality scores in the input `.bam` files. Since the base quality scores can sometimes be over-estimated, this program calculates a probability that any reference mismatches are errors instead of true SNPs. This script requires a `.vcf` file of known indels and a `.vcf` file of known SNPs.

##### 3. `run_haplotypecaller.sh`
This script runs the [HaplotypeCaller](https://github.umn.edu/neyha001/GBarleyS/blob/master/Pipeline_Scripts/run_haplotypecaller.sh) program in the GATK pipeline using the `.bam` files from the `run_baserecalibrator.sh` step. It identifies active regions (regions that may contain variants) and performs de-novo reassembly of the active region, detecting possible haplotypes. The program then performs pairwise alignment of reads against the possible haplotypes in the active region to determine the likelihood of all possible alleles at a given site. The likelihood of alleles given the read data are used to calculate the likelihood of each genotype per sample given the read data per sample. The most likely genotype is assigned to the sample. This script uses the gVCF option of the HaploytypeCaller, which is useful for cohort genotyping.

##### 4. `run_combine_gvcfs.sh`
This script runs the [CombineGVCFs](https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_variantutils_GenotypeGVCFs.php) program to merge gVCF files for each sample into hierarchical cohorts of `.g.vcf` files. The script impements a two-step process: First, `.g.vcf` files are combined by each flowcell-lane combination (e.g. C613JACXX_1) into 4 sets (for 96-plex, this means that 4 sets of 24 `.g.vcf` files are combined). Second, each set of `.g.vcf` files is combined into a single `.g.vcf` file per flowcell-lane combination. For instance, say you have one flowcell with 8 lanes (therefore 8 flowcell-lane combinations) and each lane contains 96 samples. From the `run_haplotypecaller.sh` script, you should have 8 * 96 = 768 `.g.vcf` files. The first part of `run_combine_gvcf.sh` will combine the 768 files into 4 sets of 24 files for each of 8 flowcell-lane combinations, therefore creating 4 * 8 = 32 `.g.vcf` files. The second part combines those sets into one `cohort.g.vcf` file per flowcell-lane combination. Therefore in this instance, 8 `cohort.g.vcf` files will be created. Since the intermediate set `g.vcf` files will not be needed, one might consider deleting them to save disk space. The `cohort.g.vcf` files, however, are very important and will be needed for additional genotyping as the population increases.

##### 5. `run_genotype_gvcf.sh`
This script runs the [GenotypeGVCFs](https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_variantutils_GenotypeGVCFs.php) program in the GATK pipeline. It aggregates all of the genotype likelihoods in each `.g.vcf` file and assigns genotypes at each site. The output is a single `.vcf` file of raw variants ready for filtering.



#### The README file is not yet complete. Please check back later for updates.
