#!/bin/bash

## Script to subset a directory of files in the GBarleyS pipeline
# based on the barcode ID attributed to the file (e.g. BC10). The script
# requires the sample_keyfile_barcode_id.txt file created in the demultiplexing stage,
# along with a txt file of the desired samples.

## USAGE
#directory_subset.sh [SAMPLE+BC_ID FILE] [DESIRED SAMPLES] [DIR]

# where
# DIR is the directory in which to move the files
# SAMPLE+BC_ID FILE is a txt file of sample names and barcode IDs created
## during the demultiplexing step
# DESIRED SAMPLES is a txt file where each line is a sample name to extract


bcidfile="$1"
samples="$2"
dir="$3"


# Iterate over the samples in the provided file
for sam in $(cat $samples); do

	# Find the sample in the bcid file and iterate over the resulting lines
	grep $sam $bcidfile | cut -f 1,2,5 | while read line; do

		# Concatenate the line elements into the beginning of the file
		filehead=$(echo $line | tr ' ' '_')_

		# Search for the filehead in the current directory
		for file in $(ls | grep $filehead); do

			# Move the file
			mv $file $dir

		done

	done

done
