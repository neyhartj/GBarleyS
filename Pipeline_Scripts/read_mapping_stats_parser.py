#!/usr/bin/env python

# Python script to iterate through the bowtie2 read mapping log to determine the proportion of aligned reads

import sys
import subprocess
import re
import argparse
import os

## Define functions
# Define a function to deal with the bowtie logs
def bowtie_log(bowlog):
        # Establish the list first
        toprintline = []

        for l in bowlog:
                l = l.strip()

		# Retrieve the sample name
		if l.startswith('Starting mapping for sample'):
			# This is the first line of information, so append the data
			## from the previous iteration, if present
			if toprintline:
				toprint.append(toprintline)
			toprintline = []
			tmp = l.split()
			sample = tmp[4].strip("_:")
			toprintline.append(sample)

		# Retrieve the number of reads
                if re.search('reads; of these', l):
                        tmp = l.split()
                        reads = tmp[0]
                        toprintline.append(reads)
		
		# Retrieve the number of reads that aligned zero times
		if re.search('aligned 0 times', l):
			tmp = l.split()
			reads = tmp[0]
			toprintline.append(reads)

		# Retrieve the number of reads that aligned once
                if re.search("aligned exactly 1 time", l):
                        tmp = l.split()
                        reads = tmp[0]
                        toprintline.append(reads)

		# Retrieve the number of reads that aligned once
                if re.search("aligned >1 time", l):
                        tmp = l.split()
                        reads = tmp[0]
                        toprintline.append(reads)

# Argument parser
parser = argparse.ArgumentParser()

# Add arguments
# File containing the list of zipped FastQC files
inputgroup = parser.add_mutually_exclusive_group()
inputgroup.add_argument('-i',
                '--log-file',
                metavar = 'LOG',
                help = 'The read mapping log file (This is the stderr from the run_read_mapping.sh script',
                required = False)
inputgroup.add_argument('-d',
                '--log_directory',
                metavar = 'DIR',
                help = 'Head directory containing the logs from the run_read_mapping.sh script',
                required = False)
parser.add_argument('-o',
                '--output',
                metavar = 'OUTFILE',
                help = 'Output filename (including extension)',
                required = True)
parser.add_argument('-g',
                '--gbarleys_directory',
                metavar = 'GDIR',
                help = 'Complete filepath to the GBarleyS directory',
                required = True)
parser.add_argument('-p',
                '--output_PDF',
                help = 'Optional flag to print a PDF of the graphs from the log file (requires that accompanying R script be present in the GBarleyS directory)',
                action = 'store_true',
                required = False)


# Parse the arguments
args = parser.parse_args()

# Get the output filename
filename = args.output

# Dealing with new data from the log files
toprintmap = [] # List of file data lists
# Add a header
toprintmap.append(['Filename', 'Total.reads', 'Reads.not.aligned', 'Reads.once.aligned', 'Reads.multiple.aligned'])

# Run the following if the input if a directory
if args.log_directory:

	# Collect the dirname of the directory containing the log files
	dirname = os.path.dirname(args.log_directory)

        # Create a list of directory files to work with
        filelist = []
        # Listing roots, subdirs, and files in the head directory containing the log files
        for root, subdir, files in os.walk(args.log_directory):
                filelist.append([root, subdir, files])

        newfilelist = []
        # Iterate through filelist to assemble filepaths
        for item in filelist:
                root = item[0]
                subdir = item[1]
                files = item[2]
                # Skip if files is empty
                if not files:
                        continue
                else:
                        for f in files:
                                newfilelist.append('/'.join([root, f]))

	# For each file in the list, parse it
	for f in newfilelist:
        	# Open the file
	        with open(f, 'r') as infile:
        	        # Run the function
                	toprint = []
	                bowtie_log(infile)
        	        # Append each line to the master list
			for line in toprint:
				toprintmap.append(line)



# If just a file is given, just use the file
elif args.log_file:
        newfilelist = args.log_file

	with open(newfilelist, 'r') as infile:
		# Run the function
		toprint = []
		bowtie_log(infile)
		toprintmap.append(toprint)

else:
        print "No input was given. Stopping."
        exit(1)


# Open a file to print
handle = open(filename, 'w')

# Print to the file
for line in toprintmap:
        handle.write('\t'.join(line) + '\n')

# Close the file
handle.close()

# Set the GBarleyS directory so the graphing function script can be found
graph_function = args.gbarleys_directory + '/Pipeline_Scripts/pipeline_graphing_functions.R'

if args.output_PDF:
        # Notify the user
        print "Sending the data to R and outputting a PDF."
        # Set the command for shell
        cmd = ['/soft/R/3.1.1/bin/Rscript', graph_function, filename, 'readmap']
        p = subprocess.Popen(cmd)
        p.wait()
