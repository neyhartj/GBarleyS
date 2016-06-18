#!/bin/bash

# Description
# This submission runs a script that uses the GATK GenotypeGVCF function to 
## use the combined .g.vcf files and generate variant calls based on the 
## genotype liklihoods present in those files.

# Version Info
# GATK: 3.4-46

# Inputs
## g.vcf files located in the "Combined_GVCFs" directory, with the extension
## "cohort.g.vcf"
# Outputs
# A single .vcf file with genotypes for all samples included in the combined 
## g.vcf files

##### Make changes below this line #####

# Directory in which the combined g.vcf files are located
## This directory will be "Combined_GVCFs"
INDIR=

# Directory in which the vcf file should be placed
OUTDIR=

# The name of the current project
## (e.g. "2row_GBS")
PROJECT=''

# The path to the GenomeAnalysisTK.jar file for calling GATK programs
# (e.g. /shared/software/GenomeAnalysisTK-3.4-46/GenomeAnalysisTK.jar)
GATK=

# The reference genome for GATK processing. The file should end in a FASTA
# extension (i.e. .fasta, .fa, .fas, etc)
REFERENCE=

# Set the queue settings
# It is recommended that this step be performed on a HPC node (e.g. Mesabi)
QUEUE_SETTINGS='-l walltime=06:00:00,mem=62g,nodes=1:ppn=16'

# Set the node for job submission
QUEUE=''

# Email address to receive job notifications
EMAIL=''

# Set settings for GATK IndelRealigner
# -drf DuplicateRead: disable the DuplicateRead filter (important for GBS)
GATK_SETTINGS='-drf DuplicateRead'

# File path up to and including the GBarleyS directory
## (i.e. path/to/GBarleyS)
VCPWD=

# Variable of which stage the pipeline is in
## Acceptable inputs are 'Raw' for a first-time use of the GATK pipeline without
## recalibrated data, or 'Recal' for a use of the GATK pipeline with recalibrated data
STAGE=''


#######################################
##### DO NOT EDIT BELOW THIS LINE #####
#######################################

set -e
set -u
set -o pipefail

# Error reporting
# Check if variables are empty
if [[ -z $INDIR ]] || [[ -z $OUTDIR ]] || [[ -z $PROJECT ]] || [[ -z $VCPWD ]] || [[ -z $GATK ]] || [[ -z ${GATK_SETTINGS} ]] || [[ -z $REFERENCE ]]; then
	echo -e "\nERROR: One or more variables was not specified. Please check the script and re-run." && exit 1
fi
# Check the STAGE input
if [[ $STAGE != "Raw" ]] && [[ $STAGE != "Recal" ]]; then
        echo -e "\nERROR: An improper input was given to STAGE. Acceptable inputs are 'Recal' or 'Raw'." && exit 1
fi

INPUTS=( INDIR:$INDIR OUTDIR:$OUTDIR VCPWD:$VCPWD QUEUE:${QUEUE_SETTINGS} GATK:${GATK_SETTINGS} PROJECT:$PROJECT NODE:$NODE STAGE:$STAGE )


# Change the working directory
cd $VCPWD/Pipeline/GATK_Pipeline/GenotypeGVCF 

# Make directories in the OUTDIR directory
mkdir -p $OUTDIR/$PROJECT $OUTDIR/$PROJECT/${STAGE}_Variants
OUTDIR=$OUTDIR/$PROJECT/${STAGE}_Variants

# Find the combined g.vcf files and store as an array
GVCF_LIST=( $(find $INDIR -name "*cohort.g.vcf") )

# Check to make sure that the .g.vcf files were found
# If there aren't any files in the list, exit
if [[ ${#GVCF_LIST[@]} == 0 ]]; then
        echo -e "\nNo .g.vcf files were found in the INDIR. Please make sure you have specified the correct directory." && exit 1
else
        echo -e "\nThere are ${#GVCF_LIST[@]} .g.vcf files in the INDIR."
fi

# For GATK to work, the sample list needs to be concatenated with "-V" before each file path
GVCFS=$(echo ${GVCF_LIST[@]} | sed 's/ / -V /g')

# Set a date variable for the logs
YMD=$(date +"%H%M-%m%d%y")
# Find the number of cores to use based on the QUEUE_SETTINGS
NCORES=$(echo ${QUEUE_SETTINGS} | grep -o "ppn=[0-9]*" | grep -o "[0-9]*")

# Set info for standard error and output
INFO=GenotypeGVCFs_${PROJECT}

# Send the job
# Notify the user
echo -e "Launching the job.\n"
# Launched GATK with specified inputs and settings, then pipes it to the queue
echo "cd $OUTDIR && module load java && \
java -Xmx60g -jar $GATK \
-T GenotypeGVCFs \
-R $REFERENCE \
-V $GVCFS \
${GATK_SETTINGS} \
-nt $NCORES \
-o ${PROJECT}_${STAGE}_variants.vcf \
" | qsub ${QUEUE_SETTINGS} -M $EMAIL -m abe -N $INFO -r n -q $NODE && \

echo -e "\nJob away!"

# Print a user log for records
echo -e "
Basic Info
Script: ${0}
Date: $(date)
User: $(whoami)
Project: $PROJECT

Inputs
$(echo ${INPUTS[@]} | tr ' ' '\n')

" > $(basename ${0} ".sh")_log_${YMD}.out
