#!/bin/bash

# Description
# This script uses samtools to process the raw .bam files from read mapping.
# Briefly, stats are gathered for the raw .bam files, which are then filtered
# for uniquely mapped reads (i.e. reads that align only once to the reference).
# The alignments are then sorted, stats are gathered for the processed reads,
# and the processed .bam files are indexed.

# Version info
# SamTools: 0.1.18

##### Make changes below this line #####

# Directory in which the raw .bam files are located
INDIR=

# Directory in which the processed .bam files should be placed
OUTDIR=

# The name of the current pipeline project
PROJECT=''

# The path up to and including the GBarleyS directory (e.g. path/to/GBarleyS)
VCPWD=

# Settings for the job queue
QUEUE_SETTINGS='-l walltime=06:00:00,mem=4g,nodes=1:ppn=16'

# Computing node for job submission
NODE=''

# Email address to receive notifications of job start / completion
EMAIL=''



#######################################
##### DO NOT EDIT BELOW THIS LINE #####
#######################################

set -e
set -u
set -o pipefail

# Error reporting
# Check if variables are empty
if [[ -z $INDIR ]] || [[ -z $PROJECT ]] || [[ -z $VCPWD ]] || [[ -z $OUTDIR ]] || [[ -z ${QUEUE_SETTINGS} ]] || [[ -z $NODE ]] || [[ -z $EMAIL ]]; then
	echo -e "\nERROR: One or more variables was not specified. Please check the script and re-run." && exit 1
fi

# Combine the inputs into one array
INPUTS=( INDIR:$INDIR VCPWD:$VCPWD QUEUE:${QUEUE_SETTINGS} PROJECT:$PROJECT )

# Change directory
cd $VCPWD/Pipeline/BAM_Processing

# Create an array of the raw .bam files
BAMLIST=( $(find $INDIR -name "*_raw.bam") )

# Check to see if the array contains any samples
if [[ ${#BAMLIST[@]} == 0 ]]; then
        echo -e "\nERROR: No .bam files were found in the INDIR. Please make sure the directory is correct." && exit 1
else
        echo -e "\nThere are ${#BAMLIST[@]} .bam files in the INDIR."
fi

# Make directories in the output directory and change the OUTDIR variable
mkdir -p $OUTDIR/$PROJECT/Stats $OUTDIR/$PROJECT/Processed
OUTDIR=$OUTDIR/$PROJECT

# Set a date variable for the logs
YMD=$(date +"%H%M-%m%d%y")
# Find the number of cores to use based on the QUEUE_SETTINGS
NCORES=$(echo ${QUEUE_SETTINGS} | grep -o "ppn=[0-9]*" | grep -o "[0-9]*")

# Set the info for the job
INFO=${PROJECT}_SAMTools_processing_$YMD


# Launch the job
echo "cd $OUTDIR && \
module load parallel && \
module load samtools/0.1.18 && \
function process_bam () { bam="'"$1"'"; \
bamname="'$(basename $bam "_raw.bam")'"; \
samtools flagstat "'$bam'" > $OUTDIR/Stats/"'${bamname}'"_raw_stats.out; \
samtools view -hF 4 "'$bam'" | grep -v "'"XS:i"'" | samtools view -Su - | samtools sort - $OUTDIR"'/Processed/${bamname}_processed'"; \
samtools flagstat $OUTDIR"'/Processed/${bamname}_processed.bam'" > $OUTDIR"'/Stats/${bamname}_processed_stats.out'"; \
samtools index $OUTDIR"'/Processed/${bamname}_processed.bam'"; }; \
export -f process_bam; \
parallel "'"process_bam {}"'" ::: ${BAMLIST[@]}" \
| qsub ${QUEUE_SETTINGS} -M $EMAIL -m abe -N $INFO -r n -q $NODE && \

echo -e "\nJobs away!"


# # Print a user log for records
echo -e "
Basic Info
Script: ${0}
Date: $(date)
User: $(whoami)
Project: $PROJECT

Inputs
$(echo ${INPUTS[@]} | tr ' ' '\n')

" > $(basename ${0} ".sh")_log_${YMD}.out
