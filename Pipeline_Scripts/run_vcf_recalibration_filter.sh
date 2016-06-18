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
# extension (i.e. .fasta, .fa, .fas, etc)
REFERENCE=

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

# Boolean indicating whether recalibration should be performed.
## This will override any input for KNOWn_SNPS or TRAIN_SNPS
## Acceptable inputs are 'true' or 'false'
RECALIBRATE=''

# Settings for VCFtools filtering
## Boolean indicating whether bcftools should be used to output stats
## Acceptable inputs are 'true' or 'false'
STATS=''

# The keyfile produced in the demultiplexing step. Note: this file should
## look like the following: Samples+barcode+ID_${PROJECT}.txt
## Multiple keyfiles can be included by separating the file paths with a comma (',')
KEYFILE=

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

# The file path to the GBarleyS directory
## (i.e. /path/to/GBarleyS)
VCPWD=

# Settings for PBS job submission
# Queue settings
QUEUE_SETTINGS='-l walltime=03:00:00,mem=22g,nodes=1:ppn=1'

# Computing node
QUEUE=''


#######################################
##### DO NOT EDIT BELOW THIS LINE #####
#######################################

set -e
set -u
set -o pipefail

# Error reporting
# Check if variables are empty
if [[ -z $VCFIN || -z $VCPWD || -z $GATK || -z ${GATK_SETTINGS} \
	|| -z $REFERENCE || -z ${QUEUE_SETTINGS} || -z $NODE || -z $STATS || -z $RECALIBRATE \
	|| -z $KEYFILE ]]; then
	echo -e "\nERROR: One or more variables was not specified. Please check the script and re-run." \
	&& exit 1
fi

# Check the STATS input
if [[ $STATS != "true" ]] && [[ $STATS != "false" ]]; then
        echo -e "\nERROR: An improper input was given to STATS. Acceptable inputs are 'true' or 'false'.\n" && exit 1
fi

# Check the RECALIBRATE input
if [[ $RECALIBRATE != "true" ]] && [[ $RECALIBRATE != "false" ]]; then
        echo -e "\nERROR: An improper input was given to RECALIBRATE. Acceptable inputs are 'true' or 'false'.\n" && exit 1
fi

# Save inputs as an array
INPUTS=( VCFIN:$VCFIN VCPWD:$VCPWD GAKT:$GATK GATK_SET:${GATK_SETTINGS} REF:$REFERENCE \
KEYFILE:$KEYFILE KNOWN:${KNOWN_SNPS} TRAINING:${TRAIN_SNPS} STATS:$STATS RECALIBRATE:$RECALIBRATE \
QUEUE:${QUEUE_SETTINGS} )

# Set empty variables to the default
if [[ -z $MinQ ]]; then MinQ=40; fi
if [[ -z $MinGQ ]]; then MinGQ=40; fi
if [[ -z $MinMeanDP ]]; then MinMeanDP=10; fi
if [[ -z $MinDP ]]; then MinDP=5; fi
if [[ -z $MinMAF ]]; then MinMAF=0; fi
if [[ -z $MaxMISS ]]; then MaxMISS=0; fi

# Change the working directory
cd $VCPWD/Pipeline/GATK_Pipeline/VariantRecalibration

# Copy the raw VCF file to the current directory
echo -e "\nCopying $(basename $VCFIN) to the directory $(pwd)."
cp -u $VCFIN $(pwd) && echo -e "\nVCF copied."

# Update the VCFIN variable to the local file
VCFIN=$(pwd)/$(basename $VCFIN)

# Set the name of the output VCF
## For GATK
OUTVCFG=${PROJECT}_recalibrated_variants.vcf
## For VCFtools
OUTVCFV=${PROJECT}_recalibrated_filtered_snps.vcf

# Set the output directory
OUTDIR=$(pwd)

# Pre-job VCFtools code
# This set of code will convert the names in the VCF from the barcode ID to the actual line name
# If the keyfile is not included, filtering is done without the keyfile
if [[ $RECALIBRATE ]]; then
	# Notify the user
	echo -e "\nUsing the supplied KEYFILE to change the sample names in the VCF."
	# Combine multiple keyfiles
	KEYARRAY=( $(echo $KEYFILE | tr ',' ' ') )
	# Extract the sample names from the VCF file
	NAMES=$(head -n 1000 $VCFIN | grep '#CHROM' | cut -f 10- | tr '\t' '\n')
	# Clear the newnames file if it exists
	truncate -s 0 new_vcf_sample_names.txt
	# Loop through the names in order
	for name in $NAMES; do
        	# Ensure that the name doesn't end with an underscore
	        name=$(echo $name | cut -d'_' -f 1,2,3)
        	# Use awk to find the flowcell-lane-barcode combination matching that name, 
			# then return the flowcell-lane-sample name and write it to a file
	        awk -v name=$name '{ if ($1"_"$2"_"$5 == name) print $1"_"$2"_"$5"_"$4 }' ${KEYARRAY[@]} >> new_vcf_sample_names.txt
	done
	# Capture the BLANK samples in a variable
	BAD_SAMPLES=$(echo $(cat new_vcf_sample_names.txt | grep 'BLANK') | sed 's/ / --remove-indv /g')

fi

# Separate the VCF files with known SNPs into a string for GATK input
KS=$(echo "-resource:highconfSNPs,known=false,training=true,truth=true,prior=15.0 $(echo ${KNOWN_SNPS} | sed 's/,/ -resource:highconfSNPs,known=false,training=true,truth=true,prior=15.0 /g')")
# Do the same thing for the lower confidence SNPs
TS=$(echo "-resource:lowconfSNPs,known=false,training=true,truth=false,prior=10.0 $(echo ${TRAIN_SNPS} | sed 's/,/ -resource:lowconfSNPs,known=false,training=true,truth=false,prior=10.0 /g')")


# Set the YMD variable
YMD=$(date +%m%d%y-%H%M%S)
# Find the number of cores to use based on the QUEUE_SETTINGS
NCORES=$(echo ${QUEUE_SETTINGS} | grep -o "ppn=[0-9]*" | grep -o "[0-9]*")
# Find the memory specified in the Queue settings
MEM=$(echo ${QUEUE_SETTINGS} | grep -o "mem=[0-9]*" | grep -o "[0-9]*")

# Submit the job to run the variant recalibrator
# Notify the user
echo -e "Launching the job.\n"

# If the user specifies recalibration, proceed with it
if [[ $RECALIBRATE ]]; then

	# Run the VariantRecalibrator and VCFtools filter
	echo "cd $OUTDIR && \
module load java && \
module load bcftools/1.2 && \
if [[ $STATS ]]; then bcftools stats $VCFIN > $(basename $VCFIN ".vcf")_stats.out; fi; \
java -Xmx${MEM}g -jar $GATK -T VariantRecalibrator -R $REFERENCE \
-input $VCFIN -recalFile $OUTDIR/${PROJECT}_output.recal \
-tranchesFile $OUTDIR/${PROJECT}_output.tranches ${GATK_SETTINGS} --maxGaussians 4 \
$KS $TS \
-an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an InbreedingCoeff -mode SNP; \
java -Xmx${MEM}g -jar $GATK -T ApplyRecalibration -R $REFERENCE \
-input $VCFIN -recalFile $OUTDIR/${PROJECT}_output.recal -tranchesFile $OUTDIR/${PROJECT}_output.tranches \
${GATK_SETTINGS} --ts_filter_level 99.5 -mode SNP -o $OUTDIR/$OUTVCFG; \
VCFIN=$OUTDIR/$OUTVCFG; \
bcftools reheader --samples new_vcf_sample_names.txt "'"$VCFIN"'" | \
vcftools --vcf - --remove-indels --remove-filtered-all --remove-indv ${BAD_SAMPLES} --min-alleles 2 \
--max-alleles 2 --minQ $MinQ --min-meanDP $MinMeanDP --minDP $MinDP --minGQ $MinGQ --maf $MinMAF \
--max-missing $MaxMISS --recode --recode-INFO-all --out $OUTVCFV; \
mv $OUTVCFV.recode.vcf $OUTVCFV; \
if [[ $STATS ]]; then bcftools stats $OUTVCFV > $(basename $OUTVCFV ".vcf")_stats.out; fi; \
rm new_vcf_sample_names.txt" \
| qsub ${QUEUE_SETTINGS} -M $EMAIL -m abe -N $INFO -r n -q $NODE

# Otherwise just filter using VCFtools
else

		echo "cd $OUTDIR && \
if [[ $STATS ]]; then bcftools stats $VCFIN > $(basename $VCFIN ".vcf")_stats.out; fi; \
vcftools --vcf "'"$VCFIN"'" --remove-indels --remove-filtered-all --min-alleles 2 \
--max-alleles 2 --minQ $MinQ --min-meanDP $MinMeanDP --minDP $MinDP --minGQ $MinGQ --maf $MinMAF \
--max-missing $MaxMISS --recode --recode-INFO-all --out $OUTVCFV; \
mv $OUTVCFV.recode.vcf $OUTVCFV; \
if [[ $STATS ]]; then bcftools stats $OUTVCFV > $(basename $OUTVCFV ".vcf")_stats.out; fi" \
| qsub ${QUEUE_SETTINGS} -M $EMAIL -m abe -N $INFO -r n -q $NODE

fi && echo -e "\nJob away!"



# Print a log for the user records
echo -e "
Basic Info
Script: ${0}
Date: $(date)
User: $(whoami)

Inputs
$(echo ${INPUTS[@]} | tr ' ' '\n')

VCF Filtering Parameters:
Minimum variant quality: $MinQ
Minimum genotype quality: $MinGQ
Minimum mean depth across sites: $MinMeanDP
Minimum genotype depth: $MinDP
Minimum minor allele frequency: $MinMAF
Minimum call rate: $MaxMISS

" > $(basename ${0} ".sh")_log_${YMD}.out
