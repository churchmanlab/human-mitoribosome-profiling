#!/bin/bash

#SBATCH -c 4
#SBATCH -t 0-02:00
#SBATCH -p short
#SBATCH --mem=50G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=<youremailaddress>

### WRITTEN:		3/30/21 by Mary Couvillion
###
### USE: 			This script runs featureCounts on rRNA-filtered 
###					and reports reads per gene for both unique mappers only
###					(after removing PCR duplicates), and for all mappers including
###					multi mappers (after removing PCR duplicates)

### REQUIREMENTS:	./logs/ directory 
###					*_Aligned_noDups.bam file in ./1_AlignData_CytoRP
###					
###					Rsubread package installed in R environment
###					
###					./5_CountReadsOnFeatures_CytoRP/Scripts/
###					featureCounts_hCytoRP_allSizesAsiteOverlap_CDS.R
###					
###					./Annotations/GRCh38_mitoCDSnoOverlapPlus6minus5.gtf

module load R/4.0.1

# Input from command line
LibName=$1
offset=$2 # distance to shift signal from 5' end
PWD=`pwd`

cd 5_CountReadsOnFeatures_CytoRP

############ COUNT READS ON FEATURES ###########

# This will make two .txt files for each library, one with multi mappers (after removing duplicates) and one with unique mappers (after removing duplicates). 

# Note: Here use CDS as GTF.featureType. Entire CDS is used here, instead of ignoring first and last 25 codons, as done in Taggart and Li, 2018).

# Reads are approximately transformed to A sites by shifting the read downstream before taking 5' read position for count (all read sizes are used)

./Scripts/featureCounts_hCytoRP_allSizesAsiteOverlap_CDS.R ${LibName} ${PWD} $offset

# Add headers to files
sed  -i '1i\GeneID\t'${LibName} ${LibName}_featureCounts_allSizesAsite_CDS_multi_noDups.txt
sed  -i '1i\GeneID\t'${LibName} ${LibName}_featureCounts_allSizesAsite_CDS_unique_noDups.txt
sed  -i '1i\Entry\tLength' ${LibName}_featureCounts_CDS_Length.txt

