#!/usr/bin/env python

# Python script to iterate through the FastX Toolkit and Cutadapt output logs, extract relevant information, and output that
## information to a txt file that is read by R to develop graphs.

import sys
import subprocess
import re
import argparse
import os

## Define functions
# Define a function to deal with the Cutadapt logs
def cutadapt_log(cutlog):
        # Establish the list first
        toprintline = []
	
        for l in cutlog:
                l = l.strip()

                # Get the filename of the input file
                if l.startswith('Command line parameters'):
                        if toprintline: # If it's full, append to the upper list
                                toprint.append(toprintline)
                        # This is the first entry in a line, so clear the list
                        toprintline = []
                        tmp = l.split(' ') # Split on space
                        filepath = tmp[len(tmp)-1] # Extract the filepath
                        toprintline.append(os.path.basename(filepath)) # Add the file name to the name list

                # Get the total reads processed
                if l.startswith('Total reads processed'):
                        tmp = l.split() # Split on whitespace
                        num = tmp[3] # The number is the 4th index in the line
                        num = num.replace(',', '') # Remove the commas
                        toprintline.append(num)

                # Get the number of reads with adapters
                if l.startswith('Reads with adapters'):
                        tmp = l.split() # Split on whitespace
                        num = tmp[3] # The number is the 4th index in the line
                        num = num.replace(',', '') # Remove the commas
                        toprintline.append(num)

                # Get the number of reads that were too short
                if l.startswith('Reads that were too short'):
                        tmp = l.split() # Split on whitespace
                        num = tmp[5] # The number is the 4th index in the line
                        num = num.replace(',', '') # Remove the commas
                        toprintline.append(num)

                # Get the number of reads that passed
                if l.startswith('Reads written (passing filters)'):
                        tmp = l.split() # Split on whitespace
                        num = tmp[4] # The number is the 4th index in the line
                        num = num.replace(',', '') # Remove the commas
                        toprintline.append(num)

                # Get the number of base pairs processed
                if l.startswith('Total basepairs processed'):
                        tmp = l.split() # Split on whitespace
                        num = tmp[3] # The number is the 4th index in the line
                        num = num.replace(',', '') # Remove the commas
                        toprintline.append(num)

                # Get the number of base pairs written
                if l.startswith('Total written (filtered)'):
                        tmp = l.split() # Split on whitespace
                        num = tmp[3] # The number is the 4th index in the line
                        num = num.replace(',', '') # Remove the commas
                        toprintline.append(num)


# Define the fastx toolkit function
def fastx_log(x_log):
	# Establish a blank list
	toprintline = []
	
	# Iterate through the file
	for l in x_log:
		l = l.strip()

		if l.startswith('Input'):
			if toprintline:
				toprint.append(toprintline)
			toprintline = []
                        tmp = l.split()
                        toprintline.append(tmp[1])

		if l.startswith('Output'):
                        tmp = l.split()
                        toprintline.append(tmp[1])

		if l.startswith('discarded'):
			tmp = l.split()
			toprintline.append(tmp[1])


# Argument parser
parser = argparse.ArgumentParser()

# Add arguments
# File containing the list of zipped FastQC files
inputgroup = parser.add_mutually_exclusive_group()
inputgroup.add_argument('-i',
                '--log-file',
                metavar = 'LOG',
                help = 'Either the FastX Toolkit or Cutadapt log file',
                required = False)
inputgroup.add_argument('-d',
		'--log_directory',
		metavar = 'DIR',
		help = 'Head directory containing the logs from the FastX Toolkit AND Cutadapt',
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
		help = 'Optional flag to print a PDF of the FastQC stats (requires the accompanying R script)',
		action = 'store_true',
		required = False)

# Parse the arguments
args = parser.parse_args()

# Get the output filename
filename = args.output

# Collect the dirname of the directory containing the log files
dirname = os.path.dirname(args.log_directory) 

# Check to make sure the correct directory is specified
if dirname.endswith('Cutadapt_Logs') or dirname.endswith('FastX_QT_Logs'):
	print "The specified directory is too deep. Use the containing directory!"
	exit(1)

# Dealing with new data from the log files
toprintcut = [] # List of file data lists
# Add a header
toprintcut.append(['Filename', 'Reads.processed', 'Reads.with.adapters', 'Reads.too.short', 'Reads.written', 'Basepairs.processed', 'Basepairs.written'])

toprintfast = [] # FastX log list
# Add a header
toprintfast.append(['Input.reads', 'Output.reads', 'Discarded.reads'])


# Run the following if the input if a directory
if args.log_directory:

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

# If just a file is given, just use the file
elif args.log_file:
	newfilelist = args.log_file

else:
	print "No input was given. Stopping."
	exit(1)


# For each file in the list, use either the cutadapt function or the 
## fastx function to parse it
for f in newfilelist:

	with open(f, 'r') as infile:
        	for l in infile:
                	# If the file starts with the "This is cutadapt" line, it is the cutadapt log file
	                if l.startswith('This is cutadapt'):
				# Set a blank toprint list
				toprint = []
				# Run the function
        	                cutadapt_log(infile)
				# Append the list
				for line in toprint:
					toprintcut.append(line)
                	elif l.startswith('Minimum Quality Threshold'):
				# Set a blank toprint list
				toprint = []
				# Run the function
				fastx_log(infile)
				# Append the list
				for line in toprint:
					toprintfast.append(line)
			else:
                        	break

handle = open(filename, 'w')

# Print to the file
for cline in toprintcut:
	# Get the index of the line
	ind = toprintcut.index(cline)
	# Get the corresponding line from the other list
	fline = toprintfast[ind]
	handle.write('\t'.join(cline) + '\t' + '\t'.join(fline) + '\n')

# Close the file
handle.close()

# Set the GBarleyS directory so the graphing function script can be found
graph_function = args.gbarleys_directory + '/Pipeline_Scripts/.ancillary/pipeline_graphing_functions.R'

if args.output_PDF:
	# Notify the user
	print "Sending the data to R and outputting a PDF."
	# Set the command for shell
	cmd = ['/soft/R/3.1.1/bin/Rscript', graph_function, filename, 'cutfast']
	p = subprocess.Popen(cmd)
	p.wait()
