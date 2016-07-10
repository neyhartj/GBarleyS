#!/usr/bin/env python

## This python script will take the output of the bedtools intersection
# of a VCF file with the bed file linking the physical barley reference "parts"
# positions with the full physical reference positions. This script requires the 
# bedtools intersect options to include the header (-header) and the '-wo' flag.

#####
# Define the arguments
#####

# Description
DESC = """A Python program to take the output of BEDtools intersect using the 
'150831_barley_pseudomolecules_parts_to_full_chromosomes.bed' file. The options for BEDtools
intersect should be -a [VCF] -b 150831_barley_pseudomolecules_parts_to_full_chromosomes.bed -header
-wo.\n"""

# Argument parser
parser = argparse.ArgumentParser(description=DESC, add_help=True)

# Add arguments
# Input VCF file
parser.add_argument('-i',
                '--intersect_in',
                metavar = 'INTERSECT',
                help = 'The output of BEDtools intersect run using the options described.',
                required = True)
# Output file name
parser.add_argument('-o',
                '--vcfout',
                metavar = 'VCFOUT',
                help = 'The output VCF file name without the .vcf extension',
                required = True)


# Parse the arguments
args = parser.parse_args()


### Perform python functions

# Define the output filename
outfilename = str(args.vcfout) + '.vcf'

# Open a file for writing
handle = open(outfilename, 'w')

# Read in the intersect file
with open(args.intersect_in, 'r') as vcf:

	# Interate over lines in the VCF file
	for line in vcf:

		# Simply write header lines to the handle
		if line.startswith('#'):
			handle.write(line)

        else:
        	# Split line on tabs
        	tmp = line.strip().split('\t')

        	# Subset the regular VCF info
        	vcf_info = tmp[:-7]
        	# Subset further the info to print
        	vcf_toprint = vcf_info[2:]
        	# Get the position in the chromosome part
        	part_pos = int(vcf_info[1])

        	# Subset the results from the intersection
        	intersect_info = tmp[-7:]
        	# Get the new chromosome name
        	new_chrom = intersect_info[3]
        	# Get the chromosome start position
        	chrom_start = int(intersect_info[4])

        	# Add the part_pos to the chrom_start
        	new_chrom_pos = chrom_start + part_pos

        	# Make the new toprint list
        	toprint = new_chrom + new_chrom_pos + vcf_toprint

        	# Print
        	handle.write('\t'.join(toprint) + '\n')

# Close the file
handle.close()
print "File written as " + outfilename