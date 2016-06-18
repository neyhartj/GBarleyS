#!/bin/bash

# Description
## This script runs a job that calls upon the FastQC program to assess 
# the quality of samples reads following the quality control steps. 
# The output includes html files that can be viewed on a web browser, 
# along with .zip files of the data. These .zip files are used by a 
# later script to visualize the quality

##### Make changes below this line #####

# Directory of compressed, demultiplexed fastq.gz files
INDIR=

# Name of the current pipeline project (e.g. 2row)
PROJECT=''

# Set the queue settings
QUEUE_SETTINGS='-l walltime=03:00:00,mem=4g,nodes=1:ppn=16'

# Set the node for job submission
QUEUE=''

# Email address to receive job notifications
EMAIL=''

# Path up to and including the GBarleyS directory (e.g. /path/to/GBarleyS/)
VCPWD=

# A condition whether or not to output a PDF and .txt file of combined stats on the FastQC report
## Use 'true' to output a .txt file and .pdf of FastQC stats
## Use 'false' to not output stats
STATS=''

# Variable of which stage the pipeline is in
## Acceptable inputs are 'Pre' (for running before quality control) and 'Post' (for running
## after quality control).
STAGE=''

#########################################
##### DO NOT CHANGE BELOW THIS LINE #####
#########################################

set -e
set -u
set -o pipefail

# Error reporting
# Check if variables are empty
if [[ -z $INDIR ]] || [[ -z $PROJECT ]] || [[ -z $VCPWD ]] || [[ -z $STATS ]] || [[ -z ${QUEUE_SETTINGS} ]] || [[ -z $NODE ]] || [[ -z $STAGE ]]; then
	echo -e "\nERROR: One or more variables was not specified. Please check the script and re-run." && exit 1
fi
# Check the STAGE input
if [[ $STAGE != "Pre" ]] && [[ $STAGE != "Post" ]]; then
        echo -e "\nERROR: An improper input was given to STAGE. Acceptable inputs are 'Pre' and 'Post'." && exit 1
fi

# Combine the inputs into one array
INPUTS=( INDIR:$INDIR VCPWD:$VCPWD QUEUE:${QUEUE_SETTINGS} STATS:$STATS STAGE:$STAGE PROJECT:$PROJECT )

#Change directory to the output directory
cd $VCPWD/Pipeline/Assess_Quality/${STAGE}_QC/

# Create a directory based on the project and set the outdir to that directory
mkdir -p ${PROJECT} 
OUTDIR=$(pwd)/${PROJECT}

# Make a list of the demultiplexed samples
SAMPLES=( $(find $INDIR -name "*.fastq.gz") )

# Check to see if the array contains any samples
if [[ ${#SAMPLES[@]} == 0 ]]; then
        echo -e "\nERROR: No .fastq.gz files were found in the INDIR. Please make sure the directory is correct." && exit 1
else
        echo -e "\nThere are ${#SAMPLES[@]} .fastq.gz files in the INDIR."
fi

# Set a date variable for the logs
YMD=$(date +"%H%M-%m%d%y")
# Find the number of cores to use based on the QUEUE_SETTINGS
NCORES=$(echo ${QUEUE_SETTINGS} | grep -o "ppn=[0-9]*" | grep -o "[0-9]*")

# Set info for standard error and output
INFO=Assess_quality_${STAGE}qc_${PROJECT}

# Launch the job
echo "cd $OUTDIR && \
module load fastqc/0.11.2 && \
module load parallel && \
parallel "'"fastqc {} -o $(pwd)"'" ::: ${SAMPLES[@]}; \
wait; \
find $(pwd) -name "'"*.zip"'" | sort > ${PROJECT}_FastQC_${STAGE}QC_zipfiles.txt; \
wait; \
if [[ $STATS ]]; then python $VCPWD/Pipeline_Scripts/fastqc_stats_parser.py -i ${PROJECT}_FastQC_${STAGE}QC_zipfiles.txt -o ${PROJECT}_${STAGE}qc_stats.txt -d $VCPWD -p; \
fi" \
| qsub ${QUEUE_SETTINGS} -M $EMAIL -m abe -N $INFO -r n -q $NODE && \

echo -e "\nJob away!"

# Notify user of location for zipfile .txt file
echo -e "\nThe list of ZIP files will be available at"
echo "${OUTDIR}/${PROJECT}_FastQC_${STAGE}QC_zipfiles.txt"

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
