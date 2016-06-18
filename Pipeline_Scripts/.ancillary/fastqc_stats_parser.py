#!/usr/bin/env python

# Python script to iterate through zipped FastQC results, extract relevant data,
# and compile a curated file for use in R.

import sys
import subprocess
import re
import argparse
import os
import numpy

# Argument parser
parser = argparse.ArgumentParser()

# Add arguments
# File containing the list of zipped FastQC files
parser.add_argument('-i',
                '--fastqc_zip_list',
                metavar = 'LIST',
                help = 'Input .txt file containing a list of zipped FastQC files',
                required = True)
parser.add_argument('-o',
		'--outfile',
		metavar = 'OUTFILE',
		help = 'Filename (including extension) for the curated FastQC data file',
		required = True)
parser.add_argument('-d',
                '--vcpwd_directory',
                metavar = 'DIR',
                help = 'Complete filepath to the GBarleyS directory',
                required = True)
parser.add_argument('-p',
		'--output_PDF',
		help = 'Optional flag to print a PDF of the FastQC stats (requires the accompanying R script)',
		action = 'store_true',
		required = False)

# Parse the arguments
args = parser.parse_args()

dirname = os.path.dirname(args.fastqc_zip_list) # Collect the dirname of the list of fastqc zip files

filename = args.outfile # Move the output filename to another variable

handle = open(filename, 'w') # Open the output file for writing

toprint = [] # List of file data lists
toprintqualpos = [] # List of base positions as an index

# Assembled lists of information
finqual = []
finpos = []

findist = []
finlen = []

# Header lists
basic_header = [] # Create a list for the header information on the output file
grades_header = []
quals_header = []

# Empty counters
length_min = 10000000
length_max = 0

with open(args.fastqc_zip_list, 'r') as z:
	for n, entry in enumerate(z): # Each line is a zip file
		tmp = entry.strip() # Strip the newline
		
		# Set filepath of fastqc data file
		# the extension needs to be removed
		# the file basename needs to be used
		datafile = os.path.splitext(os.path.basename(tmp))[0] + '/fastqc_data.txt'
		# Set the command for shell
		cmd = ['unzip', '-p', tmp, datafile]
		p = subprocess.Popen(cmd, stdout = subprocess.PIPE)
		out, err = p.communicate() # Store stdout as a variable


		# Split the output on newline
		data = out.split('\n')
		
		filedata = []
		grades = []
		quals = [] # List for GC content

		# Base quality material
		posobs = [] # List to store the observed positions
		qualobs = [] # List to store the observed quality scores

		lenobs = []
		lenquantobs = [] # List for seq length dist

		ender = '>>END' # This line separates tests
		# Extract data
		for i, line in enumerate(data[:-1]):

			# Continue if it says "END MODULE")
			if line.startswith(ender):
				continue
			# Pull out read length information
			if line.startswith('Sequence length'): # Pull out the range of sequence lengths
				bl = line.split('\t')[1]
				lmax = int(re.split('-', bl)[-1]) # Pull out the max
				lmin = int(re.split('-', bl)[0])

				if lmax > length_max:
					length_max = lmax
				if lmin < length_max:
					length_min = lmin
			# Pull out the filename
			if line.startswith('Filename'):
				res = line.split('\t')
				if n == 0: 
					basic_header.append(res[0])
				filedata.append(res[1])
			# Pull out the number of sequences
                        if line.startswith('Total Sequences'):
                                res = line.split('\t')
				if n == 0:
					basic_header.append(res[0])
                                filedata.append(res[1])
			# Pull out pass/fail data on tests
			if line.startswith('>>'):
				res = line.split('\t')
				if n == 0:
					grades_header.append(re.split('>>', res[0])[1]) # Append the test name to the header
				grades.append(res[1]) # Append the test results
			# Pull out the percent GC
			if line.startswith('%GC'):
				gc = line.split('\t')
				if n == 0:
					quals_header.append(gc[0])
				quals.append(gc[1])
			
			# Start looking at per-base data
			if line.startswith('>>Per base sequence quality'):				
				i = i + 2
				for l in data[i:]:
					# As long as the ender is not present, continue reading lines
					if l.startswith(ender):
						break # Break the loop
					else:
						tmp = l.split('\t')
						# Record the base position if not already present
						op = int(tmp[0].split('-')[0])
						if op not in posobs:
							posobs.append(op)
						
						#mean = tmp[1] # Mean
						med = tmp[2] # Median
						per25 = tmp[3] # Lower quartile
						per75 = tmp[4] # Upper quartile
						per10 = tmp[5] # 10 percentile
						per90 = tmp[6] # 90 percentile

						# Append the data as a list
#						bqual.append(mean)
						qualperpos = [med,per25,per75,per10,per90]
						qualobs.append(qualperpos)
									
			if line.startswith('>>Sequence Length Distribution'):
				i = i + 2
				for l in data[i:]:
					if l.startswith(ender):
						break
					else:
						tmp = l.split('\t') # Split on tab
						leni = int(tmp[0].split('-')[0]) # Record the length
						if leni not in lenobs: # If the length is already in the list, skip it
							lenobs.append(leni) # Otherwise add it to the list

						cnt = tmp[1] # Record the number of reads at that length

						# Append data
						lenquantobs.append(cnt)

		# Append quality data to multi-file list
		finqual.append([qualobs])
		finpos.append(posobs)

		findist.append([lenquantobs])
		finlen.append(lenobs)

		# Append to the large print list
		datalist = [filedata] + [grades] + [quals]
		toprint.append(datalist) # Append the large list to the print list
		

# Outside the while open loop
# Create ranges
allpos = range(1, length_max+1, 1)
alllen = range(length_min, length_max+1,1)

# Create toprint lists for base quality and read length
toprintbasequal = []
toprintreadlen = []

# Iterate through each entry in posobs and qualobs
for g, group in enumerate(finpos):
	# Create base quality list
	toprintbq = []
	for p in allpos:
        	if p in group: # If the position is in the list of observed positons
                	pind = group.index(p) # Find the position in the smaller list
                        qualp = finqual[g][0][pind] # Extract the quality info
                        toprintbq.append(qualp) # Add to new list
                        # If p is not there, pad with NA
                else:
                        toprintbq.append(['NA'] * 5)
	# Flatten the list
	toprintbq = [item for sublist in toprintbq for item in sublist]
	# Append to toprint
	toprint[g].append(toprintbq)

for g, group in enumerate(finlen):
	# Create the read count list
	toprintrl = []
	for p in alllen:
		if p in group:
			pind = group.index(p)
			distp = findist[g][0][pind]
			toprintrl.append(distp)
		else:
			toprintrl.append('NA')

	# Append to toprint
	toprint[g].append(toprintrl)

# Append to the large headers list
# Base quality positions
allpos = list(numpy.repeat(allpos, 5)) # The base positions need to be replicated by length of bqual

allpos.sort() # Sort the positions

# Combine the headers
# The map(str, list) is used to conver the integers to list-suitable
headers = basic_header + grades_header + quals_header + map(str, allpos) + map(str, alllen)

newheader = [] # Create new variable for headers without spaces
for h in headers: # Remove the spaces
	i = h.replace(' ','.')
	newheader.append(i)

# Write to the file
handle.write('\t'.join(newheader) + '\n') # Write the header
for item in toprint:
	handle.write('\t'.join('\t'.join(map(str, l)) for l in item) + '\n') # Write the data


handle.close() # Close the file

# Set the GBarleyS directory so the graphing function script can be found
graph_function = args.vcpwd_directory + '/Pipeline_Scripts/.ancillary/pipeline_graphing_functions.R'

if args.output_PDF:
	# Notify the user
	print "Sending the data to R and outputting a PDF."
	# Set the command for shell
	cmd = ['/soft/R/3.1.1/bin/Rscript', graph_function, filename, 'fastqc']
	p = subprocess.Popen(cmd)
	p.wait()

