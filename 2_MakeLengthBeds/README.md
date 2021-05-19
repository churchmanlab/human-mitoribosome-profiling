Uses ${LibName}_Aligned.Mito_mRNA.noDups.bam file from 1_AlignData to make 5' and 3' plus and minus bed files to make Vplots (read length vs. genomic position)


# 1. Run MakeLengthBeds.sh
Replace S1, S2, etc with library (sample) names, and run code block.

```bash
Libs="S1 S2 S3 S4"
scriptpath="./2_MakeLengthBeds/Scripts/MakeLengthBeds.sh"

for lib in $Libs
do
sbatch -e logs/2_MakeLengthBeds_${lib}.err -o logs/2_MakeLengthBeds_${lib}.log $scriptpath ${lib}
done
```

# 2. Outputs  
  - ${LibName}_Mito_mRNA_lengths3p_M.bed  
  - ${LibName}_Mito_mRNA_lengths3p_P.bed  
  - ${LibName}_Mito_mRNA_lengths5p_M.bed  
  - ${LibName}_Mito_mRNA_lengths5p_P.bed  
   
# 3. Vplots
Open script ./2_AlignData/Scripts/Length_vs_Pos_Vplot.R with text editor  
Modify path, experiment name, etc as needed   
Run in R in interactive node: 
```bash
source('./2_MakeLengthBeds/Scripts/Length_vs_Pos_Vplot.R')
```
