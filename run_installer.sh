#!/bin/bash

# This script begins the setup process for the GBarleyS pipeline. It creates directories and
## subdirectories needed for other parts of the pipeline.

# Usage run_installer.sh DIR

# Where DIR is the directory in which you want to place the GBarleyS head directory




#######################################
##### DO NOT EDIT BELOW THIS LINE #####
#######################################

set -e
set -u
set -o pipefail

# Test if any arguments are passed
if [ $# -eq 0 ]
	then echo "The target directory was not specified. Please specify a target directory, but not the current directory (./)." && exit 1
fi

# Set the directory variable as the first argument (excluding the script name)
DIR="$1"

# Test if the DIR variable is empty
if [ $DIR == "." ]
	then echo "Do not specify the working directory (./) as the target directory." && exit 1
fi

# Check if the current working directory contains the pipeline scripts
if ! basename $(pwd) | grep -q "GBarleyS";
	then echo "Please move to the cloned "GBarleyS" directory before running the run_installer.sh script." && exit 1
fi

# Make directories
mkdir -vp $DIR/GBarleyS && echo "Head directory created!"

# Copy the scripts to the GBarleyS directory
cp -r $(pwd)/Pipeline_Scripts $DIR/GBarleyS && echo "Scripts copied"
# Copy the resources folder
cp -r $(pwd)/.resources $DIR/GBarleyS && echo "Resources copied"

# Change directory
cd $DIR/GBarleyS

# Make other directories
mkdir -vp Pipeline \
Pipeline/Demultiplex \
Pipeline/Demultiplex/Barcode_Splitter \
Pipeline/Demultiplex/Post-Splitting \
Pipeline/Assess_Quality \
Pipeline/Assess_Quality/Pre_QC \
Pipeline/Assess_Quality/Post_QC \
Pipeline/Quality_Control \
Pipeline/Read_Mapping \
Pipeline/BAM_Processing \
Pipeline/GATK_Pipeline \
Pipeline/GATK_Pipeline/IndelRealignment \
Pipeline/GATK_Pipeline/HaplotypeCaller \
Pipeline/GATK_Pipeline/CombineGVCF \
Pipeline/GATK_Pipeline/GenotypeGVCF \
Pipeline/GATK_Pipeline/BaseRecalibration \
Pipeline/GATK_Pipeline/VariantRecalibration \
Pipeline/VCF_Processing \
&& echo "Sub-directories created!"


# Check to see if software is installed.
# Load software modules
#if ! module load bowtie2/2.2.4 
#	then 
#	echo "The bowtie2/2.2.4 module is not available. Please install from http://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.4/" && exit 1
#elif ! module load bedtools/2.17.0
#	then
#	echo "The module bedtools/2.17.0 is not available. Please install from
module load bcftools/1.2 && \
module load bedtools/2.17.0 && \
module load bowtie2/2.2.4 && \
module load cutadapt/1.8.1 && \
module load fastqc/0.11.2 && \
module load fastx_toolkit/0.0.14 && \
module load parallel && \
module load samtools/0.1.18 && \
echo "All modules are present!"

# Check if the commands are in the $PATH
((command -v bcftools) > /dev/null) || (echo "Error: bcftools not installed. Please install it from http://www.htslib.org/download/" && exit 1)

((command -v bedtools) > /dev/null) || (echo "Error: bedtools not installed. Please install it from http://bedtools.readthedocs.org/en/latest/index.html" && exit 1)

((command -v bowtie2) > /dev/null) || (echo "Error: bowtie2 not installed. Please install it from http://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.4/" && exit 1)

((command -v cutadapt) > /dev/null) || (echo "Error: cutadapt not installed. Please install it from https://cutadapt.readthedocs.org/en/stable/" && exit 1)

((command -v fastqc) > /dev/null) || (echo "Error: fastqc not installed. Please install it from http://www.bioinformatics.babraham.ac.uk/projects/fastqc/" && exit 1)

((command -v fastq_quality_trimmer) > /dev/null) || (echo "Error: fastx_toolkit not installed. Please install it from http://hannonlab.cshl.edu/fastx_toolkit/download.html" && exit 1)

((command -v fastx_barcode_splitter.pl) > /dev/null) || (echo "Error: fastx_toolkit not installed. Please install it from http://hannonlab.cshl.edu/fastx_toolkit/download.html" && exit 1)

((command -v parallel) > /dev/null) || (echo "Error: GNU parallel not installed. Please install it from http://www.gnu.org/software/parallel/" && exit 1)

((command -v samtools) > /dev/null) || (echo "Error: samtools not installed. Please install it from http://samtools.sourceforge.net/" && exit 1)

((command -v vcftools) > /dev/null) || (echo "Error: vcftools not installed. Please install it from https://vcftools.github.io/index.html" && exit 1)

echo "All software is in the PATH!"

echo -e "\nLooks like you're all set!"
