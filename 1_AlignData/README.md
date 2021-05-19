# 1. Run ProcessFASTQ_hMitoRP.sh
Replace S1, S2, etc with library (sample) names, update other variables and number of samples (e.g. for 4 samples use {0..3}), and run code block
```bash
Libs=("S1" "S2" "S3" "S4")
fastqs=("S1.fastq.gz" "S2.fastq.gz" "S3.fastq.gz" "S4.fastq.gz")
scriptpath="./1_AlignData/Scripts/ProcessFASTQ_hMitoRP.sh"
UMI="3p6"  # 3p6, 3p6_5p4, 3p10_5p4

for i in {0..3}
do
sbatch -e logs/1_AlignData_${Libs[i]}.err -o logs/1_AlignData_${Libs[i]}.log $scriptpath ${Libs[i]} ${fastqs[i]} $UMI
done
```

# 2. Combine outputs
**Do not do this until all jobs are complete**

```bash
Experiment="MitoRP1"
Libs=("S1" "S2" "S3" "S4")
```
- Read length distributions from FASTQC  
  - This will count soft-clipped nts as part of the length, not ideal for getting RPF lengths  
  - Be sure to modify paste command based on number of samples you have in your experiment
  ```bash
  list='Aligned.Mito_mRNA.noDups.uniq Aligned.Nuc_mRNA.noDups.uniq hMito_tRNA_al hMito_rRNA_al mMito_mRNAs_0mm_al'
  for type in $list
  do
  paste <(awk '{print $1"\t"$2}' ${Libs[0]}_${type}_LengthDist.txt) <(awk '{print $2}' ${Libs[1]}_${type}_LengthDist.txt) <(awk '{print $2}' ${Libs[2]}_${type}_LengthDist.txt) <(awk '{print $2}' ${Libs[3]}_${type}_LengthDist.txt) > ${Experiment}_${type}_LengthDist.txt 
  done
  ```

- Read length distributions from samstats  
  - This does not count soft-clipped nts and is more accurate for RPF lengths
  - First make file with all possible lengths to join together:  
  ```bash
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
  list='Aligned.Mito_mRNA.noDups'

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
  paste <(awk '{print $1"\t"$2}' ${Libs[0]}_${type}_samstatsReadlengthAll.noSoft.txt) <(awk '{print $2}' ${Libs[1]}_${type}_samstatsReadlengthAll.noSoft.txt) <(awk '{print $2}'  ${Libs[2]}_${type}_samstatsReadlengthAll.noSoft.txt) <(awk '{print $2}' ${Libs[3]}_${type}_samstatsReadlengthAll.noSoft.txt) > ${Experiment}_${type}_noSoftLengthDist.txt 
  done
  ```

- Counts files with read mapping stats
  ```bash
  Samps="S1 S2 S3 S4"

  echo All Counts > ${Experiment}_Counts.txt
  for Samp in $Samps
  do
  cat ${Samp}_Counts.txt >> ${Experiment}_Counts.txt
  done
  ```
  
# 3. Outputs
  - Read length distributions for different RNA species: ${Experiment}_${type}_LengthDist.txt
  - RPF length distributions: ${Experiment}_${type}_noSoftLengthDist.txt
  - Library compositions: ${Experiment}_Counts.txt
  - 5' bedGraphs: ${LibName}_Mito_mRNA.noDups.5p.plus.bedGraph, ${LibName}_Mito_mRNA.noDups.5p.plus.bedGraph

# 4. Read length distribution plots
Open script ./1_AlignData/Scripts/ReadLengthPlot.R with text editor  
Modify path, experiment name, colors, etc as needed   
Run in R in interactive node: source('./1_AlignData/Scripts/ReadLengthPlot.R')
