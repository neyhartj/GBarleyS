#!/bin/bash

# Description
# This script launches jobs that run the first step of the GATK variant calling
# pipeline: indel realignment. Briefly, reads that align near indels often result
# in artifacts that resemble evidence for SNPs. This step realigns around those 
# regions to clean up artifacts. The number of jobs launched is 2n, where n is 
# the number of flowcell-lane combinations. The processing is broken up into
# part a and part b for each flowcell-lane. This helps to save time.

# Inputs
# Processed .bam files
# GATK-formatted reference (see https://www.broadinstitute.org/gatk/guide/article?id=2798)
# A VCF file of known indels (if unavailable, or if running the pipeline for the
# first time to obtain such information, leave the input blank)

# Version info
# GATK: 3.4-46

##### Make changes below this line #####

# The directory in which the processed bam files are located
# Note: a scratch directory is recommended.
# Also Note: for HPC nodes, you need to specify the scratch space by using /scratch.global/username/path/to/indir
INDIR=

# The directory in which the resulting realigned bam files should be placed
# Note: for HPC nodes, you need to specify the scratch space by using /scratch.global/username/path/to/outdir
OUTDIR=

# The name of the current project
# e.g. "2row_GBS"
PROJECT=''

# The reference genome for GATK processing (the file should end in a FASTA 
# extension (i.e. .fasta, .fa, .fas, etc)
REFERENCE=

# The path to the GenomeAnalysisTK.jar file for calling GATK programs
# (e.g. /shared/software/GenomeAnalysisTK-3.4-46/GenomeAnalysisTK.jar)
GATK=

# A vcf file of known indels
# Leave blank if you are running the pipeline for the first time
KNOWN_INDELS=

# Email address for queue notification
EMAIL=''

# Set the queue settings
# It is recommended that this step be performed on a HPC node (e.g. Mesabi)
QUEUE_SETTINGS='-l walltime=24:00:00,mem=32gb,nodes=1:ppn=1'

# The computing node to submit the jobs to
# Note: for this script, using HPC resources is highly recommended
NODE=''

# Set settings for GATK IndelRealigner
# -drf DuplicateRead: disable the DuplicateRead filter (important for GBS)
GATK_SETTINGS='-drf DuplicateRead'

# File path up to and including the GBarleyS directory
# (i.e. path/to/GBarleyS)
VCPWD=


#############################
##### DO NOT EDIT BELOW #####
#############################

set -e
set -u
set -o pipefail

# Error reporting
# Check if variables are empty
if [[ -z $INDIR ]] || [[ -z $OUTDIR ]] || [[ -z $PROJECT ]] || [[ -z $VCPWD ]] || [[ -z $GATK ]] || [[ -z $REFERENCE ]] \
|| [[ -z $EMAIL ]] || [[ -z ${QUEUE_SETTINGS} ]] || [[ -z ${GATK_SETTINGS} ]] || [[ -z $NODE ]]; then
        echo -e "\nERROR: One or more variables was not specified. Please check the script and re-run." && exit 1
fi

# Combine the inputs into one array
INPUTS=( INDIR:$INDIR OUTDIR:$OUTDIR VCPWD:$VCPWD REFERENCE:$REFERENCE QUEUE:${QUEUE_SETTINGS} GATK_SETTINGS:${GATK_SETTINGS} GATK:$GATK \
PROJECT:$PROJECT )

# Change workding directory
cd $VCPWD/Pipeline/GATK_Pipeline/IndelRealignment

# Find the .bam files in the INDIR and store as an array
BAM_LIST=( $(find $INDIR -name "*.bam") )

# Check to make sure that the .bam files were found
# If there aren't any files in the list, exit
if [[ ${#BAM_LIST[@]} == 0 ]]; then
        echo -e "\nERROR: No .bam files were found in the INDIR. Please make sure you have specified the correct directory." && exit 1
else
        echo -e "\nThere are ${#BAM_LIST[@]} .bam files in the BAMDIR."
fi

# Make out directories for the project
mkdir -p $OUTDIR/$PROJECT/Target_Lists $OUTDIR/$PROJECT/IndelRealigned_BAMs

# Set the YMD variable for the date
YMD=$(date +%m%d%y-%H%M%S)

# Notify the user whether indel realignment will be performed with or without
# a set of known indels
if [[ -z ${KNOWN_INDELS} ]]; then
	echo -e "\nYou have not provided a set of known indels. IndelRealignment will proceed without such input."
else
	echo -e "\nYou have provided a set of known indels from the file ${KNOWN_INDELS}. IndelRealignment will proceed with this input."
fi


# Subdivide the .bam list by flowcell-lane and run jobs based on this subdivision
for fl in $(echo "${BAM_LIST[@]}" | grep -o '[A-Z0-9]\{9\}_[0-9]' | sort -u)
do

        # Create an array of all bam files in the larger list that are from
	## the specific flowcell-lane
        LANEBAMS=( $(echo "${BAM_LIST[@]}" | tr ' ' '\n' | grep $fl) )
	
	# Set the job info for the queue
        INFO=${fl}_indelrealignment_${YMD}

	# If the KNOWN_INDELS variable is blank, launch jobs without that option
	if [[ -z ${KNOWN_INDELS} ]]; then

		echo "cd $OUTDIR/$PROJECT && \
module load java && \
for SAMPLE in ${LANEBAMS[@]}; \
do OUTLIST=Target_Lists/"'$(basename $SAMPLE ".bam")'"_realign_targets.list; \
OUTBAM=IndelRealigned_BAMs/"'$(basename $SAMPLE ".bam")'"_realigned.bam; \
java -Xmx31g -jar $GATK -T RealignerTargetCreator -R $REFERENCE -I "'$SAMPLE'" ${GATK_SETTINGS} -o "'$OUTLIST'"; \
java -Xmx14g -jar $GATK -T IndelRealigner -R $REFERENCE -I "'$SAMPLE'" -targetIntervals "'$OUTLIST'" -o "'$OUTBAM'"; \
done" | qsub ${QUEUE_SETTINGS} -M $EMAIL -m abe -N $INFO -r n -q $NODE

	# If the KNOWN_INDELS variable is not blank, launch jobs with that option
	else

		echo "cd $OUTDIR/$PROJECT && \
module load java && \
for SAMPLE in ${LANEBAMS[@]}; \
do OUTLIST=Target_Lists/"'$(basename $SAMPLE ".bam")'"_realign_targets.list; \
OUTBAM=IndelRealigned_BAMs/"'$(basename $SAMPLE ".bam")'"_realigned.bam; \
java -Xmx31g -jar $GATK -T RealignerTargetCreator -R $REFERENCE -I "'$SAMPLE'" ${GATK_SETTINGS} -known ${KNOWN_INDELS} -o "'$OUTLIST'"; \
java -Xmx14g -jar $GATK -T IndelRealigner -R $REFERENCE -I "'$SAMPLE'" -targetIntervals "'$OUTLIST'" -known ${KNOWN_INDELS} -o "'$OUTBAM'"; \
done" | qsub ${QUEUE_SETTINGS} -M $EMAIL -m abe -N $INFO -r n -q $NODE

	# End the if statement
	fi
done && echo -e "\nJobs away!"


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
