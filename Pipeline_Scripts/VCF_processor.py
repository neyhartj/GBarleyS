#!/usr/bin/env python

# This python script is meant to be an inclusive processor of the VCF file
# from the GBS GATK pipeline. The script can perform two functions:
# First, if desired, the VCF file is intersected to a bedfile. If the bedfile contains the
# positions of the morex contigs within the reference, the output is either a VCF file with 
# modified variant names to reflect the morex contig location, or a hapmap file with
# the same modified variant names. If the bedfile contains the positions of the BOPA markers within
# the reference, the output is a file containing the name of the variant, the alleles of the variant,
# the name of the intersecting BOPA marker, and the alleles of said BOPA marker. Second, the script can
# simply convert a VCF file to a hapmap file, with no change in variant or entry names.

# Import modules
import subprocess # To spawn subprocesses
import sys # To take arguments
import re # To use regular expressions
import argparse # To get the arguments

#####
# Define functions
#####

def basic_hapmap(VCF):
	print "Writing the non-reformatted hapmap file using rrBLUP encoding."
	# Create filename
        filename = str(args.outfile) + '_hmp.txt'
        # Open handle for writing file
        handle = open(filename, 'w')
	
	# Lists for handling chromosome name
	chrom_l = [] # Empty list of chromosomes
	chrom_index = [] # Index of the chromosome positions

	# Start reading through the vcf file
       	for line in VCF:

               	if line.startswith('##'):
                       	continue
                # Look for the #CHROM line as this is the line that contains sample information
       	        elif line.startswith('#CHROM'):
               	        # Split the line on tabs
                       	tmp = line.strip().split('\t')
                        # Fine the format field
       	                format_field = tmp.index('FORMAT')
               	        # Get the samples out of the list
                        # Add 1 to the index, so "FORMAT" isn't included
                        samples = tmp[format_field + 1:]
			# First the header line with the information
               	        handle.write('rs\tallele\tchrom\tpos\t' + '\t'.join(samples) + '\n')

                # Now we have the sample names; let's move on to the genotype data
       	        else:
               	        # Create a new list to store the data
                       	toprint = []

                        tmp = line.strip().split('\t')

       	                # Assigning variable
               	        chrom = tmp[0]
                       	position = tmp[1]
                        ref_allele = tmp[3]
       	                alt_allele = tmp[4]

			# Handling chromosome names
			if chrom in chrom_l:
				pass
			else:
				chrom_l.append(chrom)
			# Assign the index of the chromosome within the unique list of chromosomes
			## as the name of that chromosome
			chrom_name = str(chrom_l.index(chrom) + 1)
			
               	        # The genotype data
                       	genotypes = tmp[9:]

                        # Create variable for the output file
       	                # Create the alleles variable
               	        alleles = ref_allele + '/' + alt_allele
       	                # The position variable was already created
               	        # Create a name for the SNP
                       	snp_id = 'S' + chrom_name + '_' + position

                        # Append to the list each snp_id, alleles, etc
       	                toprint.append(snp_id)
               	        toprint.append(alleles)
                       	toprint.append(chrom_name)
                        toprint.append(position)

       	                for g in genotypes:
               	                # The genotype string is separated by :
                       	        # The first element of the genotype string is the genotype call
                               	call = g.split(':')[0]
                                # Genotypes are listed as allele1/allele2
       	                        # Assume the genotypes are unphased
               	                # 0 = ref, 1 = alt 1
                       	        # 0/0 = homo ref, 0/1 = het, 1/1 = homo alt
                                
				# Encode genotypes in rrBLUP format
				individual_call =''

	                        if call == '0/0': # If the call is 0/0, declare as 1
  	                        	individual_call += '1'
                	       	elif call == '0/1': # If the call is 0/1, declare as 0
                        	        individual_call += '0'
                                elif call == '1/1': # If the call is 1/1, declare as -1
                                      	individual_call += '-1'
	                        else:
        	                        individual_call += 'NA' # If it isn't any of the above, it its missng
               	                # Append the individual calls to the genotype matrix row
                       	        toprint.append(individual_call)

                        # Print the organized list
       	                handle.write('\t'.join(toprint) + '\n')
	
	print "File was written as " + filename
	# Close the handle
	handle.close()
##### End of function #####


# Define the hapmap formatter from intersect results
def adv_hapmap(VCF):
	print "Writing the reformatted hapmap file using rrblup encoding."
        # Create filename
        filename = str(args.outfile) + '_hmp.txt'
        # Open handle for writing file
        handle = open(filename, 'w')

        # Start reading through the vcf file
        for line in VCF[:-1]:
		if line.startswith('##'):
                        continue
                elif line.startswith('#CHROM'): # Look for the #CHROM line as this is the line that contains sample information
                        # Split the line on tabs
                        tmp = line.strip().split('\t')
                        # Fine the format field
                        format_field = tmp.index('FORMAT')

                        # Get the samples out of the list
                        # Add 1 to the index, so "FORMAT" isn't included
                        samples = tmp[format_field + 1:]
                        # Create a list to store the line names
                        # This part of the script will remove the flowcell-lane-BCID part of the name
                        linename = []
                        # For each name in the samples string...
                        for name in samples:
                                # Extract the line name and append it to the linename list
                                linename.append(re.split('BC[0-9]*_', name)[1])

                        # First the header line with the information
                        handle.write('rs\tallele\tchrom\tpos\t' + '\t'.join(linename) + '\n')

                # Now we have the sample names; let's move on to the genotype data
                else:
                        # Create a new list to store the data
                        toprint = []

                        tmp = line.strip().split('\t')
                        # Assigning variable
                        # The pseudoscaffold
                        pscaffold = tmp[0]
                        # The position of the variant in the pseudoscaffold
                        position = tmp[1]
                        # The reference allele
                        ref_allele = tmp[3]
                        # The alternate allele
                        alt_allele = tmp[4]

                        # The genotype data goes from the 9th entry until the 5th from last entry
                        # (since there is intersect data at the end)
                        genotypes = tmp[9:-6]

                        # Extract intersect data
                        match = tmp[-6:]         
			
			# The contig start position
                        start = match[1]
                        # The name of the contig
                        contig = match[3]
                        # The cM position of the contig
                        centimorgan = match[4]

                        # Create variable for the output file
                        # Calculate the contig position
                        contigpos = int(position) - (int(start) + 1)
                        # Extract the pseudoscaffold name
                        pschrom = pscaffold.split('_')[3]
                        # Extract the contig number
                        contignum = contig.split('_')[2]
                        # Create the SNPid
                        snp_id = 'hvgbs_' + pschrom + '_' +  str(contignum) + ':' + str(contigpos) + '_' + ref_allele + '_' + alt_allele
                        # Create the alleles variable
                        alleles = ref_allele + '/' + alt_allele
                        # Creat the chrom variable
                        chrom = pscaffold.split('PS')[1]
                        # Check if, upon splitting the scaffold name, the chromosome is a number
                        # If it is, keep it
                        if chrom.isdigit():
                                chrom = chrom
                        # If not, call it NA
                        else:
                                chrom = "NA"
                        # Append to the list each snp_id, alleles, etc in the order they should appear
                        toprint.append(snp_id)
                        toprint.append(alleles)
                        toprint.append(chrom)
                        toprint.append(centimorgan)
			
                        for g in genotypes:
                                # The genotype string is separated by :
                                # The first element of the genotype string is the genotype call
                                call = g.split(':')[0]
                                # Genotypes are listed as allele1/allele2
                                # Assume the genotypes are unphased
                                # 0 = ref, 1 = alt 1
                                # 0/0 = homo ref, 0/1 = het, 1/1 = homo alt

                                # Encode genotypes in rrBLUP format
                                individual_call =''

                                if call == '0/0': # If the call is 0/0, declare as 1
	                                individual_call += '1'
                                elif call == '0/1': # If the call is 0/1, declare as 0
                                        individual_call += '0'
                                elif call == '1/1': # If the call is 1/1, declare as -1
                                        individual_call += '-1'
                                else:
                                        individual_call += 'NA' # If it isn't any of the above, it its missing
                                # Append the individual calls to the genotype matrix row
                                toprint.append(individual_call)

                        # Print the organized list
                        handle.write('\t'.join(toprint) + '\n')

	print "File was written as " + filename
	# Close the file
	handle.close()
##### End of function #####


# Define the contig intersect parsing function
def contig_int(inter):
	print "Writing the Morex contig intersect results."
        # Create filename
        filename = str(args.outfile) + '_contig_int.txt'
        # Open handle for writing file
        handle = open(filename, 'w')

        # Run through the output line-by-line
        for line in inter[:-1]:
                if line.startswith('##'): # Skip the header and format lines
                        continue
                elif line.startswith('#CHROM'):
                        # Print header information
                        handle.write('GBS_variant\tChrom\tMorex_contig\tPos_on_contig\n')
                else:
                        # Create a new list to store the data
                        toprint = []
                        tmp = line.strip().split('\t')

                        # Assigning variables
                        # The pseudoscaffold
                        pscaffold = tmp[0]
                        # Extract the PS number
                        chrom = pscaffold.split('PS')[1]
                        # The position of the variant in the pseudoscaffold
                        position = tmp[1]

                        # Extract intersect data
                        match = tmp[-6:]
                        # The contig start position
                        start = match[1]
                        # The name of the contig
                        contig = match[3]

                        # Create the variant name
                        name = 'S' + str(chrom) + '_' + str(position)
                        # Calculate position in the contig
                        contigpos = int(position) - (int(start) + 1)

                        # Append to the toprint list
                        toprint.append(str(name))
                        toprint.append(str(pscaffold))
                        toprint.append(str(contig))
                        toprint.append(str(contigpos))

                        # Print the toprint list
                        handle.write('\t'.join(toprint) + '\n')

	print "File was written as " + filename
        handle.close()
##### End of function #####


# Define the bopa marker intersect parsing function
def bopa_int(inter):
	print "Writing the BOPA marker intersect results."
        # Create filename
        filename = str(args.outfile) + '_bopa_int.txt'
        # Open handle for writing file
        handle = open(filename, 'w')

        # Run through the output line-by-line
        for line in inter[:-1]:
                if line.startswith('##'): # Skip the header and format lines
                        continue
                elif line.startswith('#CHROM'):
                        # Print header information
                        handle.write('GBS_variant\tChrom\tGBS_alleles\tBOPA_marker\tBOPA_alleles\n')
                else:
                        # Create a new list to store the data
                        toprint = []
                        tmp = line.strip().split('\t')

                        # Assigning variables
                        # The pseudoscaffold
                        pscaffold = tmp[0]
                        # Extract the PS number
                        chrom = pscaffold.split('PS')[1]
                        position = tmp[1]
			# Variant alleles
			ref_allele = tmp[3]
                        alt_allele = tmp[4]	

                        # Extract intersect data
                        match = tmp[-6:]
                        # The name of the contig
                        bopa = match[3]
			bopa_allele = match[4]

                        # Create the variant name
                        name = 'S' + str(chrom) + '_' + str(position)
                        # Format the GBS alleles
			gbs_allele = str(ref_allele) + '/' + str(alt_allele)                        

                        # Append to the toprint list
                        toprint.append(str(name))
                        toprint.append(str(pscaffold))
			toprint.append(str(gbs_allele))
                        toprint.append(str(bopa))
                        toprint.append(str(bopa_allele))

                        # Print the toprint list
                        handle.write('\t'.join(toprint) + '\n')

	print "File was written as " + filename
        handle.close()
##### End of function #####


#####
# Define the arguments
#####

# Description
DESC = """A Python program to intersect a VCF file using bedtools and/or 
convert a VCF file to a hapmap file in rrBLUP format {-1, 0, 1}. This tool is part of
the GBarleyS pipeline.\n\nMost options return the specified output independent of 
other options (i.e. providing the -r flag and the -b flag will return both a hapmap
file and a BOPA intersect report). The only exception is providing the -m flag and either
the -t flag or the -r flag. In this case, a reformatted hapmap file will be provided
in which the markers have been renamed to reflect their morex contig position. Note
that nothing is written to stdout; instead, files with the outfile name are written."""

# Argument parser
parser = argparse.ArgumentParser(description=DESC, add_help=True)

# Add arguments
# Input VCF file
parser.add_argument('-i',
                '--vcf_in',
                metavar = 'VCFIN',
                help = 'Input VCF file',
                required = True)
# Output file name
parser.add_argument('-o',
                '--outfile',
                metavar = 'OUTFILE',
                help = 'Output file basename (i.e. no extension)',
                required = True)
parser.add_argument('-f',
		'--hapmap',
		action = 'store_true',
		help = 'Boolean flag for whether a hapmap file should be exported',
		required = False)
parser.add_argument('-b',
		'--bedfile',
		metavar = 'BEDFILE',
		help = 'BED file for intersection',
		required = False)

# Parse the arguments
args = parser.parse_args()


#####
# Execute program functions
#####

# Print a statement if no flags are thrown
if not any((args.hapmap, args.bedfile)):
	print "\nNo options were given. What do you want to do?\n"
	parser.print_help()
	exit(1)

# Run bedtools intersect, if called
# Running bedtools using subprocess
# If the bedfile is defined, run the intersect function
if args.bedfile:
	bedfile = args.bedfile
	print "Intersecting with the " + str(bedfile) + " bedfile."
	# Set the arguments for the shell process
	cmd = ['/panfs/roc/itascasoft/bedtools/2.17.0/bin/bedtools', 'intersect', '-wo', '-header', '-a', args.vcf_in, '-b', bedfile]
	p = subprocess.Popen(cmd, stdout = subprocess.PIPE)
	out, err = p.communicate() # Store the standard out as an out variable

	# Split the string on newline characters
	inter = out.split('\n')

	# Output the file
	contig_int(inter)

# Basic conversion of the VCF to hapmap
if not args.bedfile: # If the -m flag is not thrown
	# If the haptype is not defined, don't convert the VCF
	if args.hapmap:
		with open(args.vcf_in, 'r') as VCF:
			basic_hapmap(VCF)

