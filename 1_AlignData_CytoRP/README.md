This pipeline produces a FASTQC report from adaptor-trimmed raw reads, then filters out human nuclear-encoded ncRNA (rRNA and tRNA), aligns to mito ncRNA, and oGAB control oligo to retrieve counts, then to the human genome. It removes PCR duplicates and if no UMI was used will need to be modified. It creates 5' end bedGraph (tdf) files for human ncRNA-filtered reads for visualization on IGV. See below (#3 Outputs) for the remaining files produced.  
Most outputs will be in the directory 1_AlignData\_CytoRP. Combined outputs from step 2 and plots will be in main directory.

# 1. Run ProcessFASTQ_hCytoRP.sh
Replace S1, S2, etc with library (sample) names, update UMI and number of samples (e.g. for 4 samples use {0..3}), and run code block.  
UMI choices are 3p6 (6 random nucleotides at 3' end), 3p6\_5p4 (6 at 3' end, 4 at 5' end), or 3p10\_5p4 (10 at 3' end, 4 at 5' end).
```bash
mkdir logs

Libs=("S1" "S2" "S3" "S4")
fastqs=("S1_CytoRP.fastq.gz" "S2_CytoRP.fastq.gz" "S3_CytoRP.fastq.gz" "S4_CytoRP.fastq.gz")
scriptpath="./1_AlignData_CytoRP/Scripts/ProcessFASTQ_hCytoRP.sh"
UMI="3p10_5p4"  # 3p6, 3p6_5p4, 3p10_5p4

for i in {0..3}
do
sbatch -e logs/1_AlignData_CytoRP_${Libs[i]}.err -o logs/1_AlignData_CytoRP_${Libs[i]}.log $scriptpath ${Libs[i]} ${fastqs[i]} $UMI
done
```

# 2. Combine outputs
**Do not do this until all jobs are complete**

```bash
Experiment="CytoRP1" # Replace "1" with descriptive name if desired
Libs=("S1" "S2" "S3" "S4")
```
- Read length distributions from FASTQC  
  - This will count soft-clipped nts as part of the length, not ideal for getting RPF lengths  
  - Be sure to modify paste command based on number of samples you have in your experiment
  ```bash
  cd 1_AlignData_CytoRP
  list='Aligned.Mito_mRNA.noDups.uniq Aligned.Nuc_mRNA.noDups.uniq hMito_tRNA_al hMito_rRNA_al'
  for type in $list
  do
  paste <(awk '{print $1"\t"$2}' ${Libs[0]}_${type}_LengthDist.txt) <(awk '{print $2}' ${Libs[1]}_${type}_LengthDist.txt) <(awk '{print $2}' ${Libs[2]}_${type}_LengthDist.txt) <(awk '{print $2}' ${Libs[3]}_${type}_LengthDist.txt) > ../${Experiment}_${type}_LengthDist.txt 
  done
  cd ../
  ```

- Read length distributions from samstats  
  - This does not count soft-clipped nts and is more accurate for RPF lengths
  - First make file with all possible lengths to join together:  
  ```bash
  cd 1_AlignData_CytoRP
  echo 'Length' > ${Experiment}_lengths.txt
  for x in {13..46}
  do
  echo $x >> ${Experiment}_lengths.txt
  done 
  ```

  - Join them  
    Notes:  
    -j 1: Match on column 1  
    -a1: Print all lines in file1 even if it doesn't have match in file2  
    -e 0: Fill missing with '0'  

  ```bash
  list='Aligned.Mito_mRNA.noDups Aligned.Nuc_mRNA.noDups'

  for type in $list
  do
  join -j 1 -a1 -e 0 -o 1.1,2.2 ${Experiment}_lengths.txt ${Libs[0]}_${type}_samstatsReadlength.noSoft.txt > ${Libs[0]}_${type}_samstatsReadlengthAll.noSoft.txt
  join -j 1 -a1 -e 0 -o 1.1,2.2 ${Experiment}_lengths.txt ${Libs[1]}_${type}_samstatsReadlength.noSoft.txt > ${Libs[1]}_${type}_samstatsReadlengthAll.noSoft.txt
  join -j 1 -a1 -e 0 -o 1.1,2.2 ${Experiment}_lengths.txt ${Libs[2]}_${type}_samstatsReadlength.noSoft.txt > ${Libs[2]}_${type}_samstatsReadlengthAll.noSoft.txt
  join -j 1 -a1 -e 0 -o 1.1,2.2 ${Experiment}_lengths.txt ${Libs[3]}_${type}_samstatsReadlength.noSoft.txt > ${Libs[3]}_${type}_samstatsReadlengthAll.noSoft.txt
  done
  ```

  - Combine samples into one text file
  ```bash
  for type in $list
  do
  paste <(awk '{print $1"\t"$2}' ${Libs[0]}_${type}_samstatsReadlengthAll.noSoft.txt) <(awk '{print $2}' ${Libs[1]}_${type}_samstatsReadlengthAll.noSoft.txt) <(awk '{print $2}'  ${Libs[2]}_${type}_samstatsReadlengthAll.noSoft.txt) <(awk '{print $2}' ${Libs[3]}_${type}_samstatsReadlengthAll.noSoft.txt) > ../${Experiment}_${type}_noSoftLengthDist.txt 
  done
  cd ../
  ```

- Counts files with read mapping stats
  ```bash
  cd 1_AlignData_CytoRP
  Samps="S1 S2 S3 S4"

  echo All Counts > ../${Experiment}_Counts.txt
  for Samp in $Samps
  do
  cat ${Samp}_Counts.txt >> ../${Experiment}_Counts.txt
  done
  cd ../
  ```
  
# 3. Outputs
  - FASTQC report on adapter-trimmed reads: \*\_Cleaned_fastqc.html
  - Read length distributions for different RNA species: ${Experiment}\_${type}\_LengthDist.txt
  - RPF length distributions: ${Experiment}\_${type}\_noSoftLengthDist.txt
  - Library compositions: ${Experiment}\_Counts.txt
  - 5' bedGraphs (tdfs): all mappers: \*\_noDups.5p.plus.tdf, \*\_noDups.5p.minus.tdf and unique mappers: \*\_noDups.unique.5p.plus.tdf \*\_noDups.unique.5p.minus.tdf

# 4. Read length distribution plots
Open script ./1\_AlignData/Scripts/ReadLengthPlot.R with text editor  
Modify experiment name, colors, etc as needed  
Requires the following R packages: data.table, pheatmap, RColorBrewer, stringr  
Run: 
```bash
module load R/4.0.1
Rscript ./1_AlignData/Scripts/ReadLengthPlot.R
```
