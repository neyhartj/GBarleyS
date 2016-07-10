#!/bin/bash

# Description
# This script launches jobs that combine the gVCF files produced by the GATK
## HaplotypeCaller. This is a two-step process: first, in each group of samples
## belonging to a flowcell-lane combination, .g.vcf files are combined into 4 
## sets of 24. Second, the 4 sets are combined into a single .g.vcf file. The 
## script launches n jobs, where n is the number of flowcell-lane combinations.

# Version Info
# GATK: 3.4-46

# Inputs
## .g.vcf files from the output of HaplotypeCaller
# Outputs
## n .g.vcf files, where n is the number of flowcell-lane combinations. These
## files will be located in a directory named "Combined_GVCFs" and will have
## the extension "cohort.g.vcf"

##### Make changes below this line #####

# The directory in which the .g.vcf files from HaplotypeCaller are located
INDIR=

# The directory in which the resulting combined .g.vcf files should be placed
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
QUEUE_SETTINGS='-l walltime=24:00:00,mem=62gb,nodes=1:ppn=4'

# The computing node to submit the jobs to
# Note: for this script, using HPC resources is highly recommended
QUEUE=''

# Set the settings for GATK CombineGVCF
# -drf DuplicateRead: disable the DuplicateRead filter (this is very important
## for GBS data)
GATK_SETTINGS='-drf DuplicateRead'

# File path up to and including the Barley_VCP directory
## (i.e. path/to/Barley_VCP)
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
if [[ -z $INDIR ]] || [[ -z $OUTDIR ]] || [[ -z $PROJECT ]] || [[ -z $VCPWD ]] || [[ -z $GATK ]] \
	|| [[ -z $REFERENCE ]] || [[ -z $EMAIL ]] || [[ -z ${QUEUE_SETTINGS} ]] || [[ -z ${GATK_SETTINGS} ]] \
	|| [[ -z $QUEUE ]] || [[ -z $STAGE ]]; then
        echo -e "\nERROR: One or more variables was not specified. Please check the script and re-run." && exit 1
fi
# Check the STAGE input
if [[ $STAGE != "Raw" ]] && [[ $STAGE != "Recal" ]]; then
        echo -e "\nERROR: An improper input was given to STAGE. Please correct." && exit 1
fi

INPUTS=( INDIR:$INDIR OUTDIR:$OUTDIR VCPWD:$VCPWD REFERENCE:$REFERENCE QUEUE_SET:${QUEUE_SETTINGS} \
GATK_SETTINGS:${GATK_SETTINGS} GATK:$GATK PROJECT:$PROJECT QUEUE:$QUEUE STAGE:$STAGE )


# Change directory
cd $VCPWD/Pipeline/GATK_Pipeline/CombineGVCF

# Find the .g.vcf files in the INDIR and store as an array
GVCF_LIST=( $(find $INDIR -name "*.g.vcf") )

# Check to make sure that the .g.vcf files are present
## If not, print an error and exit
if [[ ${#GVCF_LIST[@]} == 0 ]]; then
        echo -e "\nERROR: No .g.vcf files were found in the INDIR. Please make sure you have specified the correct directory." && exit 1
else
        echo -e "\nThere are ${#GVCF_LIST[@]} .g.vcf files in the INDIR."
fi

# Make output directories for the project based on the stage
mkdir -p $OUTDIR/$PROJECT $OUTDIR/$PROJECT/$STAGE $OUTDIR/$PROJECT/$STAGE/Combined_GVCFs $OUTDIR/$PROJECT/$STAGE/Intermediate_GVCFs

# Set directories where the final combined .g.vcf files will go
OUTDIRC=$OUTDIR/$PROJECT/$STAGE/Combined_GVCFs
# And where the intermediate combined g.vcf files will go
OUTDIR=$OUTDIR/$PROJECT/$STAGE/Intermediate_GVCFs

# Set the YMD variable
YMD=$(date +%m%d%y-%H%M%S)
# Find the number of cores to use based on the QUEUE_SETTINGS
NCORES=$(echo ${QUEUE_SETTINGS} | grep -o "ppn=[0-9]*" | grep -o "[0-9]*")

# Subdivide the processed bam list by lane and run jobs based on this 
## separation
# For each unique flowcell-lane name
for fl in $(echo "${GVCF_LIST[@]}" | grep -o '[A-Z0-9]\{9\}_[0-9]' | sort -u); do

        # Create an array of all g.vcf files in the larger list that are from
	## the specific flowcell-lane
	LANEGVCFS=( $(echo "${GVCF_LIST[@]}" | tr ' ' '\n' | grep $fl) )

	# Set info for the name of the job
	INFO=${fl}_combinedgvcf_${YMD}
	
	# Divide the array of .gvcf files into 4 sets
	## Remember, this needs to be in base 0
	## First find the length of the array
	NUMGVCFS=${#LANEGVCFS[@]}
	# Declare empty array
	declare -a JOBARRAY

	## If the number of .gvcf files is less than 4, simply continue because 
	# it is not worth combining so few
	if [ $NUMGVCFS -lt 4 ]; then
		echo -e "\nThere are less than 4 .gvcf files from this Flowcell-Lane combination \
and it is therefore not worth combining so few.\n"

		# Move the files to the "Combined_GVCF" directory
		for file in ${NUMGVCFS[@]}; do
			mv $file $OUTDIRC
		# Continue
		continue

	else 

		## Next find the appropriate length of the interval
		INT=$(expr $NUMGVCFS / 4)
		## Reset the number of GVCFS to base 0 (minus 1)
		NUMGVCFS=$(expr $NUMGVCFS - 1)
		## Create an array of the cutoffs
		GVCFSEQ=( $(seq 0 $INT $NUMGVCFS) )	
		## Replace the 5th array component with the number of GVCFS (minus 1)
		GVCFSEQ[4]=$NUMGVCFS
		## For loop to divide the array
		for i in $(seq 0 3); do
			FIRST=${GVCFSEQ[$i]}
			# Select the GVCFs starting at FIRST and going INT in length
			## and add the GATK -V flag
			# Add a number to the beginning to designate the set	
			JOBARRAY[$i]=${i}__$(echo ${LANEGVCFS[@]:$FIRST:$INT} | sed 's\ \--variant\g')
		done

		# Launch the jobs
		echo "cd $OUTDIR && \
module load parallel && \
module load java && \
function gvcf_combine () { SETNUM="'$(grep -o ^[0-9] <<< $1)'"; \
GVCFS="'$(sed '"'s/[0-9]__//'"' <<< $1 | sed '"'s/--variant/ --variant /g'"')'"; \
OUTFILE=${fl}_set"'${SETNUM}'".g.vcf; \
java -Xmx15g -jar $GATK -T CombineGVCFs -R $REFERENCE --variant "'$GVCFS'" ${GATK_SETTINGS} -o "'$OUTFILE'"; }; \
export -f gvcf_combine; \
parallel --no-notice -j $NCORES \"gvcf_combine {}\" ::: ${JOBARRAY[@]}; \
wait; \
cd $OUTDIRC; \
GVCFSC="'$(echo $(find '"$OUTDIR -name \"${fl}*g.vcf\""') | sed '"'s/ / --variant /g'"')'"; \
OUTFILE=${fl}_cohort.g.vcf; \
java -Xmx60g -jar $GATK -T CombineGVCFs -R $REFERENCE --variant "'$GVCFSC'" ${GATK_SETTINGS} -o "'$OUTFILE'"" \
| qsub ${QUEUE_SETTINGS} -M $EMAIL -m abe -N $INFO -r n -q $QUEUE

	fi

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
