#!/bin/bash

# Description
# This submission runs a script that uses the GATK VariantRecalibrator program
## to adjust the variant quality scores in the .vcf file produced from
## the GenotypeGVCF step. Briefly, the GATK variant calling pipeline is designed
## to be fairly lenient to maximize sensitivity, so variant filtering and
## recalibration must be done to remove the false positives. GATK recommends
## using variant quality score recalibration (VQSR) to achieve this.
## Note that this script only applied recalibration to SNPs, not indels.
## The script also uses VCFtools to filter a VCF file. This script may be used
## without VQSR in order to filter a VCF file.

# Version Info
# GATK: 3.4-46
# VCFtools: 0.1.14

##### Make changes below this line #####

# The .vcf created from the GenotypeGVCF program.
## Note: this .vcf file should end in "raw_variants.vcf"
VCFIN=

# The name of the current project
## (e.g. "2row_GBS")
PROJECT=''

# The path to the GenomeAnalysisTK.jar file for calling GATK programs
# (e.g. /shared/software/GenomeAnalysisTK-3.4-46/GenomeAnalysisTK.jar)
GATK=

# Set settings for GATK IndelRealigner
# -drf DuplicateRead: disable the DuplicateRead filter (important for GBS)
# --maxGaussians 4: This parameter helps this program run effectively. This may be a step to debug if necessary
GATK_SETTINGS='-drf DuplicateRead'

# The reference genome for GATK processing. The file should end in a FASTA
# extension (i.e. .fasta, .fa)
REFERENCE=

# The keyfile produced in the demultiplexing step. Note: this file should
## look like the following: Samples+barcode+ID_${PROJECT}.txt
## One may provide multiple keyfiles, separated by a comma (',').
KEYFILE=

# VCF file(s) of known SNPS for model training
## File(s) should contain high confident and reliable SNPs.
## In the GATK VariantRecalibrator, these SNPs will be defined as truth training SNPs
## See https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_bqsr_BaseRecalibrator.php#--knownSites for accepted file types
## Multiple VCF files may be separated with a comma (',')
KNOWN_SNPS=

# VCF file(s) of false SNPs for model tranining
## File(s) should contain confident SNPs, but which have not been verified to a high degree of reliability
## In the GATK VariantRecalibrator, these SNPs will be defined as non-truth training SNPs
## Multiple VCF files may be separated with a comma (',')
TRAIN_SNPS=

# List of annotations to recalibrate
ANNOTATIONS='-an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP'

# The file path to the GBarleyS directory
## (i.e. /path/to/GBarleyS)
VCPWD=


#######################################
##### DO NOT EDIT BELOW THIS LINE #####
#######################################

set -e
set -u
set -o pipefail

module load bcftools/1.2
module load java


# Error reporting
# Check if variables are empty
if [[ -z $VCFIN || -z $VCPWD || -z $GATK || -z ${GATK_SETTINGS} \
	|| -z $REFERENCE ]]; then
	echo -e "\nERROR: One or more variables was not specified. Please check the script and re-run." \
	&& exit 1
fi

# Save inputs as an array
INPUTS=( VCFIN:$VCFIN VCPWD:$VCPWD GAKT:$GATK GATK_SET:${GATK_SETTINGS} REF:$REFERENCE \
KNOWN:${KNOWN_SNPS} TRAINING:${TRAIN_SNPS} )

# Change the working directory
cd $VCPWD/Pipeline/GATK_Pipeline/VariantRecalibration

# Copy the raw VCF file to the current directory
echo -e "\nCopying $(basename $VCFIN) to the directory $(pwd)."
cp -u $VCFIN $(pwd) && echo -e "\nVCF copied."

# Update the VCFIN variable to the local file
VCFIN=$(pwd)/$(basename $VCFIN)

# Set the name of the output VCF
## For GATK
OUTVCFG=${PROJECT}_recalibrated_snps.vcf

# Set the output directory
OUTDIR=$(pwd)

# Separate the VCF files with known SNPs into a string for GATK input
KS=$(echo "-resource:knownSNPs,known=false,training=true,truth=true,prior=15.0 $(echo ${KNOWN_SNPS} | sed 's/,/ -resource:highconfSNPs,known=false,training=true,truth=true,prior=15.0 /g')")
# Do the same thing for the lower confidence SNPs
TS=$(echo "-resource:highconfSNPs,known=false,training=true,truth=false,prior=10.0 $(echo ${TRAIN_SNPS} | sed 's/,/ -resource:lowconfSNPs,known=false,training=true,truth=false,prior=10.0 /g')")


# Set the YMD variable
YMD=$(date +%m%d%y-%H%M%S)

# First recode the sample names
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
BAD_SAMPLES=$(echo $(cat new_vcf_sample_names.txt | grep 'BLANK') | sed 's/ / --remove-indv /g')

# Use bcftools reheader to change the sample names, then pipe the output to vcftools
bcftools reheader --samples new_vcf_sample_names.txt $VCFIN | vcftools --vcf - \
        --remove-indv ${BAD_SAMPLES} \
        --recode-INFO-all \
	--recode \
        --out renamed_no_blank \

# Run the VariantRecalibrator
java -Xmx22g -jar $GATK \
-T VariantRecalibrator \
-R $REFERENCE \
${GATK_SETTINGS} \
-input renamed_no_blank.recode.vcf \
-recalFile $OUTDIR/${PROJECT}_output.recal \
-tranchesFile $OUTDIR/${PROJECT}_output.tranches ${GATK_SETTINGS} \
--maxGaussians 4 \
$KS $TS \
-tranche 100.0 \
-tranche 99.0 \
$ANNOTATIONS \
-mode SNP \
-rscriptFile $OUTDIR/${PROJECT}_output_plots.R

# Apply variant recalibration
java -Xmx22g -jar $GATK \
-T ApplyRecalibration \
-R $REFERENCE \
-input renamed_no_blank.recode.vcf \
-recalFile $OUTDIR/${PROJECT}_output.recal \
-tranchesFile $OUTDIR/${PROJECT}_output.tranches \
${GATK_SETTINGS} \
--ts_filter_level 99.5 \
-mode SNP \
-o recal_applied.vcf

# Keep only SNPs that passed the variant reclalibration
java -Xmx22g -jar $GATK \
-T SelectVariants \
-R $REFERENCE \
-V recal_applied.vcf \
-o $OUTVCFG \
--excludeFiltered \
--selectTypeToExclude INDEL

# Delete the intermediary files
rm renamed_no_blank.recode.vcf* renamed_entries.vcf* recal_applied.vcf* samples_to_exclude.txt new_vcf_sample_names.txt


# Print a log for the user records
echo -e "
Basic Info
Script: ${0}
Date: $(date)
User: $(whoami)

Inputs
$(echo ${INPUTS[@]} | tr ' ' '\n')

GATK Annotations
$ANNOTATIONS

" > $(basename ${0} ".sh")_log_${YMD}.out
