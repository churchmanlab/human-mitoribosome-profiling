Uses \*\_Aligned.Mito\_mRNA.noDups.bam file from 1\_AlignData to make 5' and 3' plus and minus bed files to make Vplots (read length vs. genomic position). Soft-clipped reads are not counted toward length unless they are all As and at the 3' end of the read (so likely part of polyA tail)


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
  - \*\_Mito\_mRNA\_lengths3p\_P.bed  
  - \*\_Mito\_mRNA\_lengths3p\_M.bed  
  - \*\_Mito\_mRNA\_lengths5p\_P.bed  
  - \*\_Mito\_mRNA\_lengths5p\_M.bed  
   
# 3. Vplots
Requires the following R packages: scales
Modify commands below as needed and  
run: 
```bash
Experiment="MitoRP1"
libName="S1"
zoom="none" # none start stop

module load R/4.0.1
cd 2_MakeLengthBeds
Rscript ./Scripts/Length_vs_Pos_Vplot.R $Experiment $libName $zoom
mv *.png ../.
cd ../
```
Note: beware of ND6 in these plots. Reverse strand may not be coded correctly for some zoom levels
