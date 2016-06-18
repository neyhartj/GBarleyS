#!/bin/bash

#PBS -l walltime=06:00:00,mem=15g,nodes=1:ppn=8
#PBS -N VariantRecalibration
#PBS -M email
#PBS -m abe
#PBS -r n
#PBS -q lab

# Description
# This submission runs a script that uses the GATK VariantRecalibrator program
## to adjust the variant quality scores in the .vcf file produced from
## the GenotypeGVCF step. Briefly, the GATK variant calling pipeline is designed
## to be fairly lenient to maximize sensitivity, so variant filtering and
## recalibration must be done to remove the false positives. GATK recommends
## using variant quality score recalibration (VQSR) to achieve this.
## Note that this script only applied recalibration to SNPs, not indels.

# Version Info
# GATK: 3.4-46

# Inputs
## .vcf file from the GenotypeGVCF step
## A file of known SNPs (the same file used for base recalibration)
# Outputs
## A single .vcf with recalibrated variant quality scores

##### Make changes below this line #####

# The .vcf created from the GenotypeGVCF program.
## Note: this .vcf file should end in "raw_variants.vcf"
VCFIN=

# Directory in which the output .vcf file will be placed.
## The output .vcf file will have a similar name as the input
## .vcf file.
OUTDIR=

# The name of the current project
## (e.g. "2row_GBS")
PROJECT=

# The path to the GenomeAnalysisTK.jar file for calling GATK programs
# (e.g. /shared/software/GenomeAnalysisTK-3.4-46/GenomeAnalysisTK.jar)
GATK=

# The reference genome for GATK processing. The file should end in a FASTA
# extension (i.e. .fasta, .fa, .fas, etc)
REFERENCE=

# The file of known snps
## See https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_bqsr_BaseRecalibrator.php#--knownSites for accepted file types
KNOWN_SNPS=

# The file path to the Barley_VCP directory
## (i.e. /path/to/Barley_VCP)
VCPWD=

# Remember to change the email address at the top of the script to receive
## notifications when the job is complete


#######################################
##### DO NOT EDIT BELOW THIS LINE #####
#######################################

set -e
set -u
set -o pipefail

# Load the modules
module load java

cd $VCPWD/Pipeline/GATK_Pipeline/VariantRecalibration 

# Make directories in the OUTDIR directory
mkdir -p $OUTDIR/$PROJECT $OUTDIR/$PROJECT
OUTDIR=$OUTDIR/$PROJECT

OUTVCF=$(basename $VCFIN "_raw_variants.vcf")_adj_variants.vcf

# Run the VariantRecalibrator
java -Xmx14g -jar $GATK \
	-T VariantRecalibrator \
	-R $REFERENCE \
	-input $VCFIN \
	-recalFile $OUTDIR/${PROJECT}_output.recal \
	-tranchesFile $OUTDIR/${PROJECT}_output.tranches \
	-drf DuplicateRead \
	--maxGaussians 4 \ # This parameter helps this program run effectively. This may be a step to debug if necessary
	-resource:2rowSNPs,known=true,training=true,truth=true,prior=10.0 ${KNOWN_SNPS} \
	-an MQ -an MQRankSum -an ReadPosRankSum -an FS -an InbreedingCoeff \ # Apply various annotations
	-mode SNP

# Apply the recalibration to the vcf file
java -Xmx14g -jar $GATK \
        -T ApplyRecalibration \
        -R $REFERENCE \
        -input $VCFIN \
        -recalFile $OUTDIR/${PROJECT}_output.recal \
        -tranchesFile $OUTDIR/${PROJECT}_output.tranches \
        -drf DuplicateRead \
        --ts_filter_level 99.5 \
        -mode SNP \
        -o $OUTDIR/$OUTVCF

