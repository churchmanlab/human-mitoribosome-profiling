#!/bin/bash

#SBATCH -c 4
#SBATCH -t 0-2:00
#SBATCH -p short
#SBATCH --mem=5G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=<youremailaddress>


### WRITTEN:		06/10/2020 by Mary Couvillion

### USE: 			This script will first output sequence data and 
###					density data into CodonAnalysis directory, then use those files
###					to calculate mean codon density for each codon as well as individual
###					codon densities.
###					sbatch CodonOccupancy.sh libName sizeRange
###					
### REQUIREMENTS:	A site bedGraphs: *_Mito_mRNA.noDups.Asite_${sizeRange}_P.bedGraph,
###					*_Mito_mRNA.noDups.Asite_${sizeRange}_M.bedGraph
###					
###					./6_CodonOccupancy/Scripts/get_gene_sequence_densities_MTC_hMitoRP.py
###					./6_CodonOccupancy/Scripts/calculate_codon_density_MTC_hMitoRP.py
###
###					./Annotations/hMitoExonRanges.txt
###					./Annotations/hMitoGenome.fasta
###					./Annotations/hMitoGeneticCode_BattersbyMyoblast_tRNAabundance.txt
###

# Load required modules
module load gcc/6.2.0
module load python/2.7.12


libName=$1
sizeRange=$2


# Activate python virtual environment
source hMitoRP/bin/activate

# Go to directory
cd 6_CodonOccupancy

# Make input files
python ./Scripts/get_gene_sequence_densities_MTC_hMitoRP.py -f ../Annotations/hMitoExonRanges.txt -l ${libName} -g ../Annotations/hMitoGenome.fasta -p ../4_AsiteTransformation/${libName}_Mito_mRNA.noDups.Asite_${sizeRange}_P.bedGraph -m ../4_AsiteTransformation/${libName}_Mito_mRNA.noDups.Asite_${sizeRange}_M.bedGraph

# Calculate codon density
python ./Scripts/calculate_codon_density_MTC_hMitoRP.py -f ${libName}_${sizeRange}_MitoGenes_seq.txt -s ${sizeRange} -l ${libName} -c ../Annotations/hMitoGeneticCode_BattersbyMyoblast_tRNAabundance.txt

# Add header
sed -i '1i\codon\t'${libName}'_aveDensity\t'${libName}'_SD\tCodonFreq' ${libName}_${sizeRange}_CodonDensities.txt

# Clean up 
rm *MitoGenes_seq.txt
rm *AsiteFP.txt

# Deactivate python virtual environment
deactivate