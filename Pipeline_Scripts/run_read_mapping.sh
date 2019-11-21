#!/bin/bash

#Description
# This script lauches a number of jobs that each take a number of samples and
# map the reads to a reference genome. The jobs break up a large number of 
# samples into flowcell-lane combinations and run in parallel.

# Version information
# BWA: 0.7.12
# Bowtie2: 2.2.4
# Samtools: 0.1.18

##### Make changes below this line #####

# Directory of post-qc fastq files. These fastq files should be compressed and
# have their barcodes and adapter sequences removed
INDIR=

# The directory in which the output .bam files will be placed
OUTDIR=

# Name of the current pipeline project
PROJECT=''

# The path to the basename of the indexed reference genome
REFINDEX=

# Email for queue notifications
EMAIL=''

# Qsub job settings
QUEUE_SETTINGS='-l walltime=08:00:00,mem=22gb,nodes=1:ppn=16'

# Set the computing node to use
QUEUE=''

# Choose the software to use for read-mapping
# Acceptable options are "Bowtie" or "BWA"
MAPPER=''

# Bowtie settings
# -q: input files are fastq format
# -p 8: thread the process to 8 cores
# --phred33: the input files are encoded with phred +33 quality scores
# --seed 3: sets the seed to 3 for the pseudo-random number generator (this is important for reproducibility)
# all other settings default
BOWTIE_SETTINGS='-q --phred33 --seed 3'

# BWA settings
# -r Trigger re-seeding if the MEM is longer than r * 
BWA_SETTINGS='-r 1.0 -M'

# SAMtools settings
# -S: input is .sam format
# -b: output is .bam format
SAMTOOLS_SETTINGS='-bS'

# The path up to and including the GBarleyS directory (e.g. path/to/GBarleyS)
VCPWD=



#######################################
##### DO NOT EDIT BELOW THIS LINE #####
#######################################

set -e
set -u
set -o pipefail

# Error reporting
# Check if variables are empty
if [[ -z $INDIR ]] || [[ -z $OUTDIR ]] || [[ -z $PROJECT ]] || [[ -z $REFINDEX ]] || [[ -z $VCPWD ]] \
|| [[ -z $EMAIL ]] || [[ -z $QUEUE ]] || [[ -z $MAPPER ]] || [[ -z ${BWA_SETTINGS} ]] || [[ -z ${QUEUE_SETTINGS} ]] \
|| [[ -z ${BOWTIE_SETTINGS} ]] || [[ -z ${SAMTOOLS_SETTINGS} ]]; then
	echo "One or more variables was not specified. Please check the script and re-run." && exit 1
fi

# Check the MAPPER input
if [[ $MAPPER != "Bowtie" ]] && [[ $MAPPER != "BWA" ]]; then
        echo -e "\nERROR: An improper input was given to MAPPER. Acceptable inputs are 'Bowtie' and 'BWA'." && exit 1
fi

# Combine the inputs into one array
INPUTS=( INDIR:$INDIR OUTDIR:$OUTDIR REFERENCE:$REFINDEX VCPWD:$VCPWD MAPPER:$MAPPER QUEUE:${QUEUE_SETTINGS} \
BOWTIE:${BOWTIE_SETTINGS} BWA:${BWA_SETTINGS} SAMTOOLS:${SAMTOOLS_SETTINGS} PROJECT:$PROJECT )


# Change directory
mkdir -p $VCPWD/Pipeline/Read_Mapping/$PROJECT/Read_Map_Logs/
cd $VCPWD/Pipeline/Read_Mapping/$PROJECT

# Create an array of post-quality control fastq.gz files
FASTQLIST=( $(find $INDIR -name "*BC[0-9]*.fastq.gz") )
# Check to make sure there are files in the INDIR directory
if [[ ${#FASTQLIST[@]} == 0 ]]; then
        echo -e "\nERROR: No _fastq.gz files were found in the INDIR. Please make sure you have specified the correct directory." && exit 1
else
        echo -e "\nThere are ${#FASTQLIST[@]} _fastq.gz files in the INDIR."
fi


# Create output directories
mkdir -p ${OUTDIR}/ReadMap_Output/${PROJECT} # For the .bam files
# Set a new OUTDIR path
OUTDIR=${OUTDIR}/ReadMap_Output/${PROJECT}


# Set the date
YMD=$(date +%m%d%y-%H%M%S)
# Find the number of cores to use based on the QUEUE_SETTINGS
NCORES=$(echo ${QUEUE_SETTINGS} | grep -o "ppn=[0-9]*" | grep -o "[0-9]*")


# Launch a job for each flowcell-lane combination
for fl in $(echo "${FASTQLIST[@]}" | grep -o '[A-Z0-9]\{9\}_[0-9]' | sort -u); do

	# Create an array of the .fastq files belonging to that flowcell-lane
	LANEFASTQS=( $(echo "${FASTQLIST[@]}" | tr ' ' '\n' | grep $fl) )

	# Set the INFO for the job queue
	INFO=${fl}_readmapping_${YMD}

	if [[ $MAPPER == "Bowtie" ]]; then

	# The submission runs read-mapping on each sample in each list,
	# then pipes the output to samtools view to convert the .sam file 
	# to a compressed .bam file

		echo "cd $OUTDIR && \
module load bowtie2/2.2.4 && \
module load samtools/0.1.18 && \
for fastq in ${LANEFASTQS[@]}; do \
RGID="'$(basename $fastq "QC.fastq.gz")'"; \
SMID=SM:"'$RGID'"; \
OUTPUT="'${RGID}'"raw.bam; \
>&2 echo "'"Starting mapping for sample ${RGID}: $(date)"'"; \
bowtie2 -p $NCORES ${BOWTIE_SETTINGS} --rg-id "'$RGID'" --rg "'$SMID'" --rg PL:Illumina --rg LB:${PROJECT} -x $REFINDEX -U "'$fastq'" | samtools view ${SAMTOOLS_SETTINGS} - -o "'$OUTPUT'"; \
>&2 echo "'"Mapping complete: $(date)"'"; \
done" \
| qsub ${QUEUE_SETTINGS} -e $VCPWD/Pipeline/Read_Mapping/$PROJECT/Read_Map_Logs/${INFO}_log.out -M $EMAIL -m abe -N $INFO -r n -q $QUEUE


	# Else run BWA
	else
		echo "cd $OUTDIR && \
module load bwa/0.7.12 && \
module load samtools/0.1.18 && \
for fastq in ${LANEFASTQS[@]}; do \
RGID="'$(basename $fastq "QC.fastq.gz")'"; \
OUTPUT="'${RGID}'"raw.bam; \
>&2 echo "'"Starting mapping for sample ${RGID}: $(date)"'"; \
bwa mem -t $NCORES ${BWA_SETTINGS} -R "'@RG\tID:$RGID\tSM:$RGID\tPL:'"$PLATFORM"'\tLB:'"${PROJECT} $REFINDEX "'$fastq'" | samtools view ${SAMTOOLS_SETTINGS} - -o "'$OUTPUT'"; \
>&2 echo "'"Mapping complete: $(date)"'"; \
done" \
| qsub ${QUEUE_SETTINGS} -e $VCPWD/Pipeline/Read_Mapping/$PROJECT/Read_Map_Logs/${INFO}_log.out -M $EMAIL -m abe -N $INFO -r n -q $QUEUE

	fi

done

echo "Jobs away!"


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
