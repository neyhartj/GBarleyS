#!/bin/bash

# Description
# The "run_barcode_splitter.sh" script results in demultiplexed, but uncompressed fastq files. This script compresses those fastq files

##### Make changes below this line #####

# The directory containing the uncompressed fastq files (this should be in scratch space)
INDIR=

# The current pipeline project (i.e. 2row_TP_092315)
PROJECT=''

# Email for queue notification
EMAIL=''

# The path to the "GBarleyS" folder (this should be in your home directory)
VCPWD=


##### Program Settings #####
# Queue settings for MSI
QUEUE_SETTINGS='-l walltime=01:00:00,mem=1gb,nodes=1:ppn=8'

# Name of the computing node to use
QUEUE=''


########################################
###### DO NOT EDIT BELOW THIS LINE #####
########################################

set -e
set -u
set -o pipefail

# Error reporting
# Check if variables are empty
if [[ -z $INDIR ]] || [[ -z $PROJECT ]] || [[ -z $VCPWD ]] || [[ -z $EMAIL ]] || [[ -z ${QUEUE_SETTINGS} ]]; then

	echo "One or more variables was not specified. Please check the script and re-run." && exit 1
fi

# Save inputs as an array
INPUTS=( INDIR:$INDIR VCPWD:$VCPWD QUEUE:${QUEUE_SETTINGS} PROJECT:$PROJECT )

# Change directory
cd $VCPWD/Pipeline/Demultiplex/Post-Splitting/

# First create an array of the uncompressed fastq files
echo -e "\nGathering .fastq files"
FASTQLIST=( $(find $INDIR -name "*BC[0-9]*.fastq") )
# Check to make sure there are files in the INDIR directory
if [[ ${#FASTQLIST[@]} == 0 ]]; then
        echo -e "\nERROR: No .fastq files were found in the INDIR. Please make sure you have specified the correct directory." && exit 1
else
        echo -e "\nThere are ${#FASTQLIST[@]} .fastq files in the INDIR."
fi

# Remove unmatched.fastq files
UNMATCHED=( $(find $INDIR -name "*unmatched.fastq") )
if [[ ${#UNMATCHED[@]} != 0 ]]; then
        echo -e "Removing ${#UNMATCHED[@]} .fastq files."
	# Loop to remove files
	for file in ${UNMATCHED[@]}; do
		rm $file
	done
fi


# Set the YMD variable
YMD=$(date +%m%d%y-%H%M%S)
# Find the number of cores to use based on the QUEUE_SETTINGS
NCORES=$(echo ${QUEUE_SETTINGS} | grep -o "ppn=[0-9]*" | grep -o "[0-9]*")

# Subdivide the .fastq list by flowcell-lane and run jobs based on this subdivision
for fl in $(echo "${FASTQLIST[@]}" | grep -o '[A-Z0-9]\{9\}_[0-9]' | sort -u); 
do

        # Create an array of all fastq files in the larger list that are from
        ## the specific flowcell-lane
        LANEFQS=( $(echo "${FASTQLIST[@]}" | tr ' ' '\n' | grep $fl) )

        # Set the job info for the queue
        INFO=${fl}_post-demultiplex_${YMD}

        echo "module load parallel; \
parallel -j $NCORES \"gzip {}\" ::: ${LANEFQS[@]}" \
| qsub "${QUEUE_SETTINGS}" -m abe -M $EMAIL -N $INFO -q $NODE

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




