These scripts use outputs from 5\_CountReadsOnFeatures(\_CytoRP/\_RNAseq) to produce barplots of RPK values across OXPHOS complexes (similar to Figure S4), cumulative distribution of relative synthesis rates of OXPHOS complexes (similar to Figure 4B), average synthesis of cyto- vs mito-translated OXPHOS complexes, and average RNA abundance of cyto- vs mito-translated OXPHOS complexes (similar to Figure 4C,D).  



# 1. Make directories for output plots and txt files
```bash
mkdir 7_OXPHOSstoichiometry/OXPHOSbarplots
mkdir 7_OXPHOSstoichiometry/CumFracPlots
mkdir 7_OXPHOSstoichiometry/ComplexAveCorrelations
```

# 2. Run StoichPlots_from_RPK  
Requires the following R packages: data.table, plotrix, rlist, scales  
Replace mito-, cyto-, and RNAseq- Experiment names below, as well as samplesOfinterest and soiName.  
Run code block for each set of replicates (samplesOfinterest), e.g. S1 and S2 may be two WT replicates and S3 and S4 may be two mutant replicates.  
For this script to run as-is, sample names have to be identical in all experiments (mito, cyto, and RNAseq)

```bash

mitoExperiment="hMitoRP_HEK"
cytoExperiment="CytoRP_HEK"
RNAseqExperiment="RNAseq_HEK"

samplesOfinterest="WT1,WT2"
soiName="WT"

mitoExperiment="MitoRP1"
cytoExperiment="CytoRP1"
RNAseqExperiment="RNAseq1"

samplesOfinterest="S1,S2"
soiName="WT"

datas="RP RNAseq"

for data in $datas
do
sbatch -p short -t 0-0:30 -e logs/StoichPlots.err -o logs/StoichPlots.log --wrap="Rscript ./7_OXPHOSstoichiometry/Scripts/StoichPlots_from_RPK.R $data $mitoExperiment $cytoExperiment $RNAseqExperiment $samplesOfinterest $soiName"
done
```

 
# 3. Outputs
  - Cyto and mito OXPHOS barplots for RP and RNAseq, and matching .txt files with values
  - Cyto and mito OXPHOS cumulative fraction of subunits synthesis rate relative to rate limiting (lowest synthesis) subunit
  - Cyto vs. mito average complex synthesis and RNA abundance
  