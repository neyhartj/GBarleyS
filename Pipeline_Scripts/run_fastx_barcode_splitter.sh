#!/bin/bash

# Description
# This script demultiplexes a number of multiplexed fastq files based on the barcode sequences that 
## appears at the beginning of the read. This script calls on the FastX_Barcode_Splitter, 
## which only allows barcodes of a constant length. Therefore, this script will run n jobs where n 
## is the number of different barcode lengths.

# Version info
## FastX_Toolkit: 0.0.14

##### Make changes below this line #####

# The directory where the multiplexed samples are located or where their symbolic links are located
INDIR=

# The directory to output the demultiplexed fastq files
OUTDIR=

# The key file with columns: Flowcell Lane Barcode Sample PlateName Row Column
KEYFILE=

# The name of the current project (e.g. 2row_TP)
PROJECT=''

# Email address for queue notification
EMAIL=''

# The path to the "GBarleyS" folder
VCPWD=

##### Program Settings #####
# Queue settings for MSI
QUEUE_SETTINGS='-l walltime=04:00:00,mem=16gb'

# Specify the computing node
QUEUE=''

# FastX_barcode_splitter settings
# See http://hannonlab.cshl.edu/fastx_toolkit/commandline.html#fastx_barcode_splitter_usage for more information
# bol: match barcodes at the beginning of a read
# mismatches 1: 1 mismatch allowed when matching the barcodes
FASTX_SETTINGS='--bol --mismatches 1'



#######################################
##### DO NOT EDIT BELOW THIS LINE #####
#######################################

set -e
set -u
set -o pipefail

# Error reporting
# Check if variables are empty
if [[ -z $KEYFILE ]] || [[ -z $INDIR ]] || [[ -z $OUTDIR ]] || [[ -z $PROJECT ]] || [[ -z $VCPWD ]] || \
	[[ -z $EMAIL ]] || [[ -z ${QUEUE_SETTINGS} ]] || [[ -z $QUEUE ]] || [[ -z ${FASTX_SETTINGS} ]]; then
	echo "One or more variables was not specified. Please check the script and re-run." && exit 1
fi

# Save inputs as an array
INPUTS=( INDIR:$INDIR OUTDIR:$OUTDIR VCPWD:$VCPWD KEYFILE:$KEYFILE QUEUE:${QUEUE_SETTINGS} FASTX:${FASTX_SETTINGS} PROJECT:$PROJECT )

# Create an array of the multiplexed fastq files (e.g. C631EAXX_1_fastq.gz)
FASTQLIST=( $(find $INDIR -name "*_fastq.gz") )
# Check to make sure there are files in the INDIR directory
if [[ ${#FASTQLIST[@]} == 0 ]]; then
        echo -e "\nERROR: No _fastq.gz files were found in the INDIR. Please make sure you have specified the correct directory." && exit 1
else
        echo -e "\nThere are ${#FASTQLIST[@]} _fastq.gz files in the INDIR."
fi


# Change working directory
cd $VCPWD/Pipeline/Demultiplex/Barcode_Splitter

# Create an array with the 96 barcode sequences	and sort by unique barcodes and the length of the barcode
BARCODES_EXT=( $(awk 'NR > 1 {print $3}' ${KEYFILE} | sort -u) ) \
&& echo "Keyfile $KEYFILE successfully read."

# Clear the new output file if it exists
truncate -s 0 Samples+barcode+ID_${PROJECT}.txt

#Create a file with just the IDs and the barcodes
for i in $(seq ${#BARCODES_EXT[@]}); do
	
	# Create a variable of the barcode identifier
	ID=BC$i
	#UNIX is base 0 so 1 must be subtracted from i
	j=$((${i}-1))
	# Extract the sequence for the ith barcode
	BCODE=${BARCODES_EXT[$j]}
	# Match the barcode to the correct row and output that row
	awk -v bc=$BCODE '{ if ($3 == bc) print }' $KEYFILE | awk -v id=$ID '{print $1"\t"$2"\t"$3"\t"$4"\t"id}' >> Samples+barcode+ID_${PROJECT}.txt

done

# Sort the file by flowcell
sort -k 1,1 < Samples+barcode+ID_${PROJECT}.txt > temp && mv temp Samples+barcode+ID_${PROJECT}.txt

# Add the header back to the file
HEADER=$(awk 'NR == 1 {print $1"\t"$2"\t"$3"\t"$4"\tBarcode_ID"}' $KEYFILE)
sed -i "1s/^/${HEADER}\n/" Samples+barcode+ID_${PROJECT}.txt \
&& echo "Sample, barcode, and ID file created"

# Create a variable of the sample+barcode+ID file
BCODEFILE=Samples+barcode+ID_${PROJECT}.txt

#Create an array of the barcode lengths
BCODE_LENGTH=( $(awk '{print length($3)}' $BCODEFILE | sort -u) )

# Make a directory to store the barcode files
mkdir -p Barcode_Files

# Create an array of fastq files based on length of the barcode sequence
# For each barcode length...
for len in ${BCODE_LENGTH[@]}; do

	# Use awk to find the barcodes in the sample+barcode+id file that match the length "len"
	# Additionally, output only the barcode ID and the barcode sequence, separated by a tab
	# Write this output to a new file
	awk -v l=$len 'NR > 1 { if (length($3) == l) print $5"\t"$3 }' $BCODEFILE > Barcode_Files/Barcode_file_length${len}.txt
	
	echo "Barcode file created for length $len"

done && echo "All barcode files created."

# Create an array of the barcode files
BCODEFILES=( $(find $(pwd)/Barcode_Files -name "*.txt") )


# Create the output head directory
mkdir -p ${OUTDIR}/Barcode_Splitter_Output/${PROJECT} && echo "Output directories created."
# Set the output directory as a variable
OUTDIR=${OUTDIR}/Barcode_Splitter_Output/${PROJECT}

# Determine the number of cores to use
# We will use the number of unique barcode lengths as a proxy
NCORES="${#BCODE_LENGTH[@]}"

# Adjust the queue settings
QUEUE_SETTINGS="${QUEUE_SETTINGS},nodes=1:ppn=${NCORES}"

# Set date
YMD=$(date +%m%d%y-%H%M%S)

# This code will start qsub submissions, one for each unparsed sample (e.g. C631EAXX_1_fastq.gz). Each job iterates through the list of barcode files and runs the fastx_barcode_splitter.pl script. 

# For each fastq file...
for fastq in ${FASTQLIST[@]}; do
	
	# Set a prefix variable to be placed in the beginning of each output fastq file name
	prefix=$(basename $fastq "_fastq.gz")

	echo "mkdir -p ${OUTDIR}/${prefix}; \
cd ${OUTDIR}/${prefix}; \
module load parallel; \
module load fastx_toolkit/0.0.14; \
function barcode_split () { zcat $fastq | fastx_barcode_splitter.pl --bcfile "'"$1"'" ${FASTX_SETTINGS} --prefix ${prefix}_ --suffix "'".fastq"'" >> ${prefix}_log.txt; \
}; \
export -f barcode_split; \
parallel --no-notice -j ${NCORES} \"barcode_split {}\" ::: ${BCODEFILES[@]}" \
| qsub "${QUEUE_SETTINGS}" -M $EMAIL -N ${prefix}_Barcode_Splitting_$YMD -m abe -r n -q $QUEUE

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
