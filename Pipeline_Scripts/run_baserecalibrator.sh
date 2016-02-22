#!/bin/bash

# Description
# This script launches a series of jobs that run the GATK BaseRecalibrator 
## on the indel-realigned bam files. Briefly, the base quality scores from a
## sequencing machine often overestimate the quality of the bases. This function
## uses known SNPs to readjust those scores to achieve better variant calling
## results. Note that this step should not be performed if a file of known
## SNPs and Indels is not available. If this is the case, the GATK pipeline
## should be completed without the BaseRecalibrator (i.e. IndelRealignment > 
## HaplotypeCaller > CombinedGVCF > GenotypeGVCF), the variants should be 
## hard-filtered based on some criteria, and the GATK pipeline should be
## completed again, but using the filtered variants as "known variants."
## This script launches n jobs, where n is the number of flowcell-lane
## combinations.

# Version Info
# GATK: 3.4-46

# Inputs
## realigned .bam file from the IndelRealigner (using known Indels)
## a file with known SNPs
## a file with known Indels
# Outputs
## recalibrated .bam files for use in HaplotypeCaller

##### Make changes below this line #####

# The directory in which the indel-realigned .bam filed are located
INDIR=


# The directory in which the recalibrated .bam files should be placed
## Note: it is recommended that this be a scratch directory
OUTDIR=

# The name of the current project
# (e.g. "2row_GBS")
PROJECT=''

# The reference genome for GATK processing. The file should end in a FASTA
## extension (i.e. .fasta, .fa, .fas, etc)
REFERENCE=

# The path to the GenomeAnalysisTK.jar file for calling GATK programs
# (e.g. /shared/software/GenomeAnalysisTK-3.4-46/GenomeAnalysisTK.jar)
GATK=

# Input file(s) of known indels OR known SNPs
## More than one file may be included by separating the filepaths
## with a comma (,). 
## See https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_bqsr_BaseRecalibrator.php#--knownSites for accepted file types
KNOWN_VARIANTS=

# Email address for queue notification
EMAIL=''

# Set the queue settings
## It is recommended that this program be run on an HPC node (e.g. Mesabi)
QUEUE_SETTINGS='-l walltime=12:00:00,mem=62gb,nodes=1:ppn=8'

# The computing node to submit the jobs to
# Note: for this script, using HPC resources is highly recommended
NODE=''

# Set the settings for GATK BaseRecalibrator
# -drf DuplicateRead: disable the DuplicateRead filter (this is very important
## for GBS data)
GATK_SETTINGS='-drf DuplicateRead'

# Path up to and including the Barley_VCP directory
## (i.e. /path/to/Barley_VCP)
VCPWD=


#######################################
##### DO NOT EDIT BELOW THIS LINE #####
#######################################

set -e
set -u
set -o pipefail

# Error reporting
# Check if variables are empty
if [[ -z $INDIR ]] || [[ -z $OUTDIR ]] || [[ -z $PROJECT ]] || [[ -z $VCPWD ]] || [[ -z $GATK ]] || [[ -z $REFERENCE ]] || [[ -z ${KNOWN_VARIANTS} ]] || [[ -z $EMAIL ]] || [[ -z ${QUEUE_SETTINGS} ]] || [[ -z ${GATK_SETTINGS} ]] || [[ -z $NODE ]]; then
        echo "\nERROR: One or more variables was not specified. Please check the script and re-run." && exit 1
fi

# Combine the inputs into one array
INPUTS=( INDIR:$INDIR OUTDIR:$OUTDIR VCPWD:$VCPWD REFERENCE:$REFERENCE QUEUE:${QUEUE_SETTINGS} \
KNOWN_VARIANTS:${KNOWN_VARIANTS} GATK_SETTINGS:${GATK_SETTINGS} GATK:$GATK \
PROJECT:$PROJECT )

# Change directory
cd $VCPWD/Pipeline/GATK_Pipeline/BaseRecalibration

# If multiple files of known variants are present, concatenate them with
## the correct flag for GATK
KNOWN_VARIANTS=$(echo ${KNOWN_VARIANTS} | sed 's\,\ -knownSites \g')

# Create a list of the realigned bam files in the INDIR directory
BAM_LIST=( $(find $INDIR -name "*.bam") )

# Check to make sure that the .bam files are present
## If not, print an error and exit
if [[ ${#BAM_LIST[@]} == 0 ]]; then
        echo -e "\nERROR: No .bam files were found in the INDIR. Please make sure you have specified the correct directory." && exit 1
else
        echo -e "\nThere are ${#BAM_LIST[@]} .bam files in the INDIR."
fi

# Make output directories for the project
mkdir -p $OUTDIR/$PROJECT $OUTDIR/$PROJECT/DataTables $OUTDIR/$PROJECT/Recalibrated_BAMs
# Set the OUTDIR
OUTDIR=$OUTDIR/$PROJECT

# Set the YMD variable
YMD=$(date +%m%d%y-%H%M%S)
# Find the number of cores to use based on the QUEUE_SETTINGS
NCORES=$(echo ${QUEUE_SETTINGS} | grep -o "ppn=[0-9]*" | grep -o "[0-9]*")


# Subdivide the processed bam list by lane and run jobs based on this 
## separation
# For each unique flowcell-lane name
for fl in $(echo "${BAM_LIST[@]}" | grep -o '[A-Z0-9]\{9\}_[0-9]' | sort -u)
do
        # Create an array of all bam files in the larger list that are from
	## the specific flowcell-lane
        LANEBAMS=( $(echo "${BAM_LIST[@]}" | tr ' ' '\n' | grep $fl) )
	
	# Set the job info for the queue
        INFO=${fl}_baserecalibrator_${YMD}

	# Launch the jobs
	echo "cd $OUTDIR && \
module load java && \
module load parallel && \
function baserecal () { SAMPLE="'"$1"'"; \
OUTLIST=DataTables/"'$(basename $SAMPLE "processed_realigned.bam")recal_data.table'"; \
OUTFILE=Recalibrated_BAMs/"'$(basename $SAMPLE "processed_realigned.bam")recal.bam'"; \
java -Xmx15g -jar $GATK -T BaseRecalibrator -R $REFERENCE -I "'$SAMPLE'" -knownSites ${KNOWN_VARIANTS} ${GATK_SETTINGS} -o "'$OUTLIST'"; \
java -Xmx15g -jar $GATK -T PrintReads -R $REFERENCE -I "'$SAMPLE'" ${GATK_SETTINGS} -BQSR "'$OUTLIST'" -o "'$OUTFILE'"; \
}; \
export -f baserecal; \
parallel -j $NCORES --no-notice "'"baserecal {}"'" ::: ${LANEBAMS[@]}" \
| qsub ${QUEUE_SETTINGS} -M $EMAIL -m abe -N $INFO -r n -q $NODE

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
