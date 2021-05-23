#!/bin/bash

#SBATCH -c 4
#SBATCH -t 0-01:00
#SBATCH -p short
#SBATCH --mem=10G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=<youremailaddress>

### WRITTEN:		3/30/21 by Mary Couvillion
###
### USE: 			This script runs featureCounts on mito mRNA mappers 
###					and reports reads per gene for both unique mappers only
###					(after removing PCR duplicates), and for all mappers including
###					multi mappers (after removing PCR duplicates)


### REQUIREMENTS:	./logs/ directory 
###					*_Aligned.Mito_mRNA.noDups.bam file
###					
###					Rsubread package installed in R environment
###					
###					./5_CountReadsOnFeatures/Scripts/
###					featureCounts_hMitoRP_allSizesAsiteOverap_CDS.R
###					
###					./Annotations/GRCh38_mitoCDSnoOverlapPlus6minus5.gtf


module load R/4.0.1

# Input from command line
LibName=$1
offset=$2 # distance to shift signal from 3' end
PWD=`pwd`

cd 5_CountReadsOnFeatures

############ COUNT READS ON FEATURES ###########


# Note: Here use CDS as GTF.featureType. This gtf is modified to remove first 6 codons and last 5 codons of CDS, as well as all overlapping regions


./Scripts/featureCounts_hMitoRP_allSizesAsiteOverlap_CDS.R ${LibName} ${PWD} $offset

# Add headers to files
sed  -i '1i\GeneID\t'${LibName} ${LibName}_mito_featureCounts_allSizesAsite_CDSignore1stlast_multi_noDups.txt
sed  -i '1i\GeneID\t'${LibName} ${LibName}_mito_featureCounts_allSizesAsite_CDSignore1stlast_unique_noDups.txt
sed  -i '1i\Entry\tLength' ${LibName}_mito_featureCounts_CDSignore1stlast_Length.txt



