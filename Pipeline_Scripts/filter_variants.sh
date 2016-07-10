#!/bin/bash

# Description
# This submission runs a script that uses bcftools and vcftools to filter the .vcf
## file for high-quality SNPs. Be sure to have bcftools and vcftools installed
## for this script to work. Note that if you are using MSI, bcftools is an included
## module, however vcftools is not. You may download a copy of the vcftools 
## software from https://vcftools.github.io/index.html.
## This script also changes the identifier of each sample from "Flowcell_Lane_BC##"
## to the sample name found in the keyfile.

# Inputs
## .vcf file from the VariantRecalibration step
## The keyfile produced in the demultiplexing step (see below)
# Outputs
## A single, filtered .vcf file with only high-quality, biallelic SNP markers

##### Make changes below this line #####

# The complete filepath to the VCF file
VCFIN=

# The desired name of the VCF file (including the .vcf extension)
## (e.g. '2row_GBS_filtered_snps.vcf')
VCFOUT=

# The keyfile produced in the demultiplexing step. Note: this file should
## look like the following: Samples+barcode+ID_${PROJECT}.txt
## One may provide multiple keyfiles, separated by a comma (',').
KEYFILE=

# The file path to the Barley_VCP directory
## (i.e. /path/to/GBarleyS)
VCPWD=

# Should the script use bcftools to output stats on the vcf files?
# Acceptable inputs are 'true' or 'false'
STATS=''

# You may change the filtering parameters below:
# Note: by default, this filtering script will extract only biallelic 
## SNP variants.

# Default parameters are assigned already 
# Any variables left blank will be assigned default values 

# Minimum mapping quality score (phred scaled) # to include that variant
MinQ=40
# Minimum genotype quality score (phred scaled) to include that genotype
MinGQ=40
# The minimum average read depth for a variant to include that variant
MinMeanDP=10
# The minimum per-genotype read depth to include that genotype
MinDP=5
# The minimum minor allele frequency to retain a site
MinMAF=0
# The minimum non-missing data proportion to retain a site (where 0 allows 
## completely missing data and 1 retricts to no missing data)
MaxMISS=0

# Other parameters can be changed as well. You may edit the script after the line beginning
## with "bcftools reheader...". Please see https://vcftools.github.io/man_latest.html for
## a list of adjustable parameters.


####################################################
##### USE CAUTION WHEN EDITING BELOW THIS LINE #####
####################################################

set -e
set -u
set -o pipefail

# Error reporting
# Check if variables are empty
if [[ -z $VCFIN || -z $VCFOUT || -z $VCPWD || -z $KEYFILE || -z $STATS ]]; then
	echo -e "\nERROR: One or more variables was not specified. Please check the script and re-run.\n" && exit 1
fi
# Check the STAGE input
if [[ $STATS != "true" ]] && [[ $STATS != "false" ]]; then
        echo -e "\nERROR: An improper input was given to STATS. Acceptable inputs are 'true' or 'false'.\n" && exit 1
fi

# Save inputs as an array
INPUTS=( VCFIN:$VCFIN VCFOUT:$VCFOUT VCPWD:$VCPWD KEYFILE:$KEYFILE STATS:$STATS )

# Set empty variables to the default
if [[ -z $MinQ ]]; then MinQ=40; fi
if [[ -z $MinGQ ]]; then MinGQ=40; fi
if [[ -z $MinMeanDP ]]; then MinMeanDP=10; fi
if [[ -z $MinDP ]]; then MinDP=5; fi
if [[ -z $MinMAF ]]; then MinMAF=0; fi
if [[ -z $MaxMISS ]]; then MaxMISS=0; fi

# Check if the keyfile is correct


# Change directory
cd $VCPWD/Pipeline/VCF_Processing

# Copy the raw VCF file to the current directory
echo -e "\nCopying $(basename $VCFIN) to the directory $(pwd)."
cp -u $VCFIN $(pwd) && echo -e "\nVCF copied."

# Update the VCFIN variable to the local file
VCFIN=$(pwd)/$(basename $VCFIN)

# Run stats on the raw variants
if [[ $STATS ]]; then
	echo -e "\nUsing bcftools to output stats on ${VCFIN}."
	bcftools stats $VCFIN > $(basename $VCFIN ".vcf")_stats.out \
	&& echo -e "\nPre-filter stats complete."
fi

# This set of code will convert the names in the VCF from the barcode ID to the actual line name
# Notify
echo -e "\nUsing the provided KEYFILE to change sample names."

# Extract the sample names from the VCF file
NAMES=$(head -n 1000 $VCFIN | grep '#CHROM' | cut -f 10- | tr '\t' '\n')

# Clear the file if it exists
truncate -s 0 new_vcf_sample_names.txt

# Loop through the names in order
for name in $NAMES; do
        # Ensure that the name doesn't end with an underscore
        name=$(echo $name | cut -d'_' -f 1,2,3)
        # Use awk to find the flowcell-lane-barcode combination matching that name, then return the flowcell-lane-sample name and write it to a file
        awk -v name=$name '{ if ($1"_"$2"_"$5 == name) print $1"_"$2"_"$5"_"$4 }' $KEYFILE >> new_vcf_sample_names.txt
done

# Capture the BLANK samples in a variable
# THIS NEEDS TO BE CHANGED
BAD_SAMPLES=$(echo $(cat new_vcf_sample_names.txt | grep 'BLANK') | sed 's/ / --remove-indv /g')

# Use bcftools reheader to change the sample names, then pipe the output to vcftools for filtering
echo -e "\nStarting VCF filtering."
bcftools reheader --samples new_vcf_sample_names.txt $VCFIN | vcftools --vcf - \
	--remove-indels \
	--remove-filtered-all \
	--remove-indv ${BAD_SAMPLES} \
	--min-alleles 2 \
	--max-alleles 2 \
	--minQ $MinQ \
	--min-meanDP $MinMeanDP \
	--minDP $MinDP \
	--minGQ $MinGQ \
	--maf $MinMAF \
	--max-missing $MaxMISS \
	--recode \
	--recode-INFO-all \
	--out $VCFOUT \
&& echo -e "\nFiltering complete and output file created."

# Rename the file
mv $VCFOUT.recode.vcf $VCFOUT

# Obtain stats on the output
if [[ $STATS ]]; then
	echo -e "\nUsing bcftools to output stats on ${VCFOUT}."
	bcftools stats $VCFOUT > $(basename $VCFOUT ".vcf")_stats.out \
	&& echo -e "\nPost-filtration stats complete."
fi

# Remove the new sample names file
rm new_vcf_sample_names.txt

# Notify the user
echo -e "\nFiltering complete. The final VCF file can be found at: $(pwd)/$VCFOUT.\n\nStats prior to filtering can be found at: $(pwd)/$(basename $VCFIN ".vcf")_stats.out.\n\nStats after filtering can be found at: $(pwd)/$(basename $VCFOUT ".vcf")_stats.out."

YMD=$(date +%m%d%y-%H%M%S)

# Print a log for the user records
echo -e "
Basic Info
Script: ${0}
Date: $(date)
User: $(whoami)

Inputs
$(echo ${INPUTS[@]} | tr ' ' '\n')

VCF Filtering Parameters
Minimum variant quality: $MinQ
Minimum genotype quality: $MinGQ
Minimum mean depth across sites: $MinMeanDP
Minimum genotype depth: $MinDP
Minimum minor allele frequency: $MinMAF
Minimum call rate: $MaxMISS
" > $(basename ${0} ".sh")_log_${YMD}.out
