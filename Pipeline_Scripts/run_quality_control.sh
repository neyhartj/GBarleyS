#!/bin/bash

# Description
## This script launches the quality control procedures to trim reads.
# The first step is running cutadapt, which removes the provided adapter
# from the reads and removes the barcode sequence from the beginning of the
# read. The next step is the fastx_quality_trimmer, which will trim reads
# based on a provided quality threshold and minimum read length.

# Version info
# FastX_Toolkit: version 0.0.14
# Cutadapt: version 1.8.1

##### Make changes below this line #####

# The directory in which the raw fastq.gz files are located.
# These fastq.gz files should be demultiplexed and still retain their barcodes and adapters
INDIR=

# The directory in which the post-qc fastq files should be placed
OUTDIR=

# The name of the current pipeline project
PROJECT=

# Email address for queue notification
EMAIL=''

# Set the settings for the Qsub submission
QUEUE_SETTINGS='-l walltime=04:00:00,mem=22gb,nodes=1:ppn=16'

# Set the computing node to use
QUEUE=''

# Name of the rare restriction enzyme (e.g. PstI or MspI)
RES='PstI'

# Set the settings for cutadapt
# The adapter sequence. The listed adapted is the solexa-reverse adapter
# Be sure to check with your genotyping facility for the correct adapter
# Multiple adapters may be added by separating them with a comma (',')
# For example, ADAPTER='TATA','CAAT'
ADAPTER='AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAG'
# Other cutadapt settings
# -f fastq: the file format is fastq
# -m 30: minimum read length of 30
# --quality-base=33: the base quality score is phred+33 scale
CUTADAPT_SETTINGS='-f fastq -m 30 --quality-base=33'

# Set the fastq_quality_trimmer settings
# -Q 33: the base quality score is phred+33 scale
# -t 30: trim bases with a quality score less than 30
# -l 30: minimum read length of 30
# -v: verbose
# -z: compress output with GZIP
FASTX_SETTINGS='-Q 33 -t 30 -l 30 -v -z'

# The path up to and including the GBarleyS folder (e.g. path/to/GBarleyS)
VCPWD=


#########################################
##### DO NOT CHANGE BELOW THIS LINE #####
#########################################

set -e
set -u
set -o pipefail

# Error reporting
# Check if variables are empty
if [[ -z $INDIR ]] || [[ -z $OUTDIR ]] || [[ -z $PROJECT ]] || [[ -z $RES ]] || \
	[[ -z $VCPWD ]] || [[ -z $QUEUE ]] || [[ -z $EMAIL ]] || [[ -z $ADAPTER ]] || \
	[[ -z ${QUEUE_SETTINGS} ]] || [[ -z ${FASTX_SETTINGS} ]] || [[ -z ${CUTADAPT_SETTINGS} ]]; then
	echo "One or more variables was not specified. Please check the script and re-run." && exit 1
fi

# Combine the inputs into one array
INPUTS=( INDIR:$INDIR OUTDIR:$OUTDIR RES:$RES VCPWD:$VCPWD QUEUE:${QUEUE_SETTINGS} FASTX:${FASTX_SETTINGS} CUTADAPT:${CUTADAPT_SETTINGS} PROJECT:$PROJECT )


# Change working directory
cd $VCPWD/Pipeline/Quality_Control


# Read in the enzyme cut info file
RESINFO=$VCPWD/.resources/.enzyme_site_info
# Make sure the user-provided restriction site is in the file
if ! cut -f 1 $RESINFO | grep -q -m 1 $RES; then 
	echo "The provided RES is not acceptable. Please refer to the GBarleyS documentation for \
	acceptable restriction enzyme names." && exit 1
fi

# Find the length of the restriction site
RESLEN=$(awk -v res=$RES '{ if ($1 == res) print length($2) }' $RESINFO)
# Calculate the number of additional bases to remove from the beginning of the read
# The script will use the length of the barcode + the length of the restriction site - 1
RESLEN=$((${RESLEN} - 1))

# Deal with multiple adapter input
ADAPTER=$(echo $ADAPTER | sed 's/,/ -a /g')


# Create an array of gzipped and compressed samples
FASTQLIST=( $(find $INDIR -name "*BC[0-9]*.fastq.gz") )
# Check to make sure there are files in the INDIR directory
if [[ ${#FASTQLIST[@]} == 0 ]]; then
        echo -e "\nERROR: No _fastq.gz files were found in the INDIR. Please make sure you have specified the correct directory." && exit 1
else
        echo -e "\nThere are ${#FASTQLIST[@]} _fastq.gz files in the INDIR."
fi


# Read in the samples+barcodes+id file from the run_fastx_barcode_splitter.sh script
BCODEFILE=$VCPWD/Pipeline/Demultiplex/Barcode_Splitter/Samples+barcode+ID_${PROJECT}.txt
# Create an array of the barcode lengths
BCODE_LENGTH=( $(awk '{print length($3)}' $BCODEFILE | sort -u) )


# Create directories for the log files
mkdir -p $(pwd)/$PROJECT $(pwd)/$PROJECT/Cutadapt_Logs $(pwd)/$PROJECT/FastX_QT_Logs
# Create a output directory if it doesn't already exist
mkdir -p $OUTDIR/Quality_Control_Fastq
OUTDIR=$OUTDIR/Quality_Control_Fastq

# Set a date variable for the logs
YMD=$(date +"%H%M-%m%d%y")
# Find the number of cores to use based on the QUEUE_SETTINGS
NCORES=$(echo ${QUEUE_SETTINGS} | grep -o "ppn=[0-9]*" | grep -o "[0-9]*")


# Iterate through each barcode length...
for len in ${BCODE_LENGTH[@]}; do
	# Extract all barcode IDs with length "len" and store as an array
	bcodeids=( $(awk -v l=$len 'NR > 1 { if ( length($3) == l ) print $5}' $BCODEFILE | sort -u) )
	fastqs=()
	# For each barcode ID...
	for bcid in ${bcodeids[@]}; do
		# Find all fastq files with that barcode ID and add them to the array
		fastqs+=( $(echo ${FASTQLIST[@]} | tr ' ' '\n' | grep ${bcid}.fastq.gz) )
	done

	# Submit the job to perform quality control based on the fastqs list
	# Calculate the cut length for cutadapt
	CUT_LENGTH=$((${len} + ${RESLEN}))
	# Standard out and err
	COUT=$VCPWD/Pipeline/Quality_Control/$PROJECT/Cutadapt_Logs/Cutadapt_log_length_${CUT_LENGTH}_${PROJECT}_${YMD}.out
        # Fastq_qualitry_trimmer log output is standard out - capture it in a file
        FOUT=$VCPWD/Pipeline/Quality_Control/$PROJECT/FastX_QT_Logs/Fastq_quality_trimmer_log_length_${CUT_LENGTH}_${PROJECT}_${YMD}.out

	# Begin job submission
	echo "cd $OUTDIR; \
module load parallel; \
module load cutadapt/1.8.1; \
module load fastx_toolkit/0.0.14; \
function qual_control () { SAMPLE="'"$1"'"; \
SAMPLE_OUT="'$(basename $SAMPLE ".fastq.gz")_QC.fastq.gz'"; \
cutadapt -a $ADAPTER ${CUTADAPT_SETTINGS} -u ${CUT_LENGTH} "'$SAMPLE'" | fastq_quality_trimmer ${FASTX_SETTINGS} -o "'${SAMPLE_OUT}'"; \
}; \
export -f qual_control; \
parallel -j $NCORES \"qual_control {}\" ::: ${fastqs[@]}" \
| qsub ${QUEUE_SETTINGS} -M $EMAIL -m abe -N QC_barcode_length_${CUT_LENGTH} -e $COUT -o $FOUT -r n -q $QUEUE


done && echo -e "Jobs away!\n\nThe quality trimmed fastq.gz files can be found here: ${OUTDIR}/Quality_Control_Fastq.\n\nThe logs from Cutadapt can be found here: $VCPWD/Pipeline/Quality_Control/$PROJECT/Cutadapt_Logs.\n\nThe logs from the Fastq_Quality_Trimmer can be found here: $VCPWD/Pipeline/Quality_Control/$PROJECT/FastX_QT_Logs."

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
