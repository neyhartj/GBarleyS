#!/bin/bash

# Description
# This script launches jobs to run the HaplotypeCaller function from GATK.
## Briefly, HaplotypeCaller performs de novo genome assembly, then calls SNPs
## in the active region. This function is run in GVCF mode, where genotype
## liklihoods are assigned. This is amenable for large, scalable genotyping
## (i.e. reapeated addition of samples in cohorts). Following steps output
## genotype calls for the whole collection of GVCF files.
# This script runs 2n jobs, where n is the number of flowcell-lane combinations

# Inputs
## Indel realigned or base recalibrated .bam files
# Outputs
## g.vcf files for each sample
# Intermediate/Temporary files
##

# Version info
# GATK: 3.4-46

##### Make changes below this line #####

# The directory in which the indel-realigned or base-recalibrated .bam files
## are located
INDIR=

# The directory in which the .g.vcf files should be placed
## Note: a scratch directory is recommended
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

# Email address for queue notifications
EMAIL=''

# Set the queue settings
## It is recommended that this step be performed on an HPC node (e.g. Mesabi)
QUEUE_SETTINGS='-l walltime=24:00:00,mem=62gb,nodes=1:ppn=12'

# The computing node to submit the jobs to
# Note: for this script, using HPC resources is highly recommended
QUEUE=''

# Set the settings for GATK HaplotypeCaller
# --genotyping_mode DISCOVERY: the most likely alternate allele will be chosen
# --emitRefConfidence GVCF: HaplotypeCaller will produce a .g.vcf file with 
## confidence levels for each site. This is done per-sample, then the g.vcf 
## files will be merged downstream.
# -drf DuplicateRead: disable the DuplicateRead filter (this is very important
## for GBS data)
GATK_SETTINGS='--genotyping_mode DISCOVERY --emitRefConfidence GVCF -drf DuplicateRead'

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
if [[ -z $INDIR ]] || [[ -z $OUTDIR ]] || [[ -z $PROJECT ]] || [[ -z $VCPWD ]] || [[ -z $GATK ]] || [[ -z $REFERENCE ]] || [[ -z $EMAIL ]] || [[ -z ${QUEUE_SETTINGS} ]] || [[ -z ${GATK_SETTINGS} ]] || [[ -z $NODE ]] || [[ -z $STAGE ]]; then
        echo -e "\nERROR: One or more variables was not specified. Please check the script and re-run." && exit 1
fi
# Check the STAGE input
if [[ $STAGE != "Raw" ]] && [[ $STAGE != "Recal" ]]; then
	echo -e "\nERROR: An improper input was given to STAGE. Please correct." && exit 1
fi

INPUTS=( INDIR:$INDIR OUTDIR:$OUTDIR VCPWD:$VCPWD REFERENCE:$REFERENCE QUEUE:${QUEUE_SETTINGS} \
GATK_SETTINGS:${GATK_SETTINGS} GATK:$GATK PROJECT:$PROJECT NODE:$NODE STAGE:$STAGE )

# Change directory
cd $VCPWD/Pipeline/GATK_Pipeline/HaplotypeCaller

# Find the .bam files in the INDIR and create an array
BAM_LIST=( $(find $INDIR -name "*.bam") )

# Check to make sure that the .bam files were found
# If there aren't any files in the list, exit
if [[ ${#BAM_LIST[@]} == 0 ]]; then
        echo -e "\nERROR: No .bam files were found in the INDIR. Please make sure you have specified the correct directory." && exit 1
else
        echo -e "\nThere are ${#BAM_LIST[@]} .bam files in the BAMDIR."
fi

# Make output directories for the project depending on the stage
mkdir -p $OUTDIR/$PROJECT $OUTDIR/$PROJECT/${STAGE}_Variants
OUTDIR=$OUTDIR/$PROJECT/${STAGE}_Variants

# Set the YMD variable
YMD=$(date +%m%d%y-%H%M%S)
# Find the number of cores to use based on the QUEUE_SETTINGS
NCORES=$(echo ${QUEUE_SETTINGS} | grep -o "ppn=[0-9]*" | grep -o "[0-9]*")

# Subdivide the .bam list by flowcell-lane and run jobs based on this subdivision
for fl in $(echo "${BAM_LIST[@]}" | grep -o '[A-Z0-9]\{9\}_[0-9]' | sort -u)
do

        # Create an array of all bam files in the larger list that are from
        ## the specific flowcell-lane
        LANEBAMS=( $(echo "${BAM_LIST[@]}" | tr ' ' '\n' | grep $fl) )

        # Set the job info for the queue
        INFO=${fl}_haplotypecaller_${YMD}

        # Launch the jobs
        echo "cd $OUTDIR && \
module load java && \
module load parallel; \
function haplocall () { SAMPLE="'"$1"'"; \
OUTFILE="'$(basename $SAMPLE ".bam")raw.g.vcf'"; \
java -Xmx4g -jar $GATK -T HaplotypeCaller -R $REFERENCE -I "'$SAMPLE'" ${GATK_SETTINGS} -o "'$OUTFILE'"; }; \
export -f haplocall; \
parallel --no-notice -j $NCORES \"haplocall {}\" ::: ${LANEBAMS[@]}" \
| qsub ${QUEUE_SETTINGS} -M $EMAIL -m abe -N $INFO -r n -q $NODE

done && echo "Jobs away!"

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
