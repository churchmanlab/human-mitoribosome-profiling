This pipeline produces a FASTQC report from adaptor-trimmed raw reads, then filters out human nuclear-encoded ncRNA (rRNA and tRNA), aligns to mito ncRNA, oGAB control oligo, and mouse mito mRNA to retrieve counts, then to the concatenated human and mouse genomes. It removes PCR duplicates and if no UMI was used will need to be modified. It creates coverage bedGraph (tdf) files for human rRNA-filtered reads for visualization on IGV. See below (#3 Outputs) for the remaining files produced.  
Most outputs will be in the directory 1_AlignData\_RNAseq. Combined outputs from step 2 and plots will be in main directory.

# 1. Run ProcessFASTQ_hRNAseq.sh
Replace S1, S2, etc with library (sample) names, update UMI and number of samples (e.g. for 4 samples use {0..3}), and run code block.  
UMI choices are 3p6 (6 random nucleotides at 3' end), 3p6\_5p4 (6 at 3' end, 4 at 5' end), or 3p10\_5p4 (10 at 3' end, 4 at 5' end).
```bash
mkdir logs

Libs=("S1" "S2" "S3" "S4")
fastqs=("S1_RNAseq.fastq.gz" "S2_RNAseq.fastq.gz" "S3_RNAseq.fastq.gz" "S4_RNAseq.fastq.gz")
scriptpath="./1_AlignData_RNAseq/Scripts/ProcessFASTQ_hRNAseq.sh"
UMI="3p10_5p4"  # 3p6, 3p6_5p4, 3p10_5p4

for i in {0..3}
do
sbatch -e logs/1_AlignData_RNAseq_${Libs[i]}.err -o logs/1_AlignData_RNAseq_${Libs[i]}.log $scriptpath ${Libs[i]} ${fastqs[i]} $UMI
done
```

# 2. Combine outputs
**Do not do this until all jobs are complete**

```bash
Experiment="RNAseq1" # Replace "1" with descriptive name if desired
Libs=("S1" "S2" "S3" "S4")
```

- Counts files with read mapping stats
  ```bash
  cd 1_AlignData_RNAseq
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
  - Library compositions: ${Experiment}\_Counts.txt
  - Coverage tdfs (for loading on IGV): all mappers: \*\_Aligned\_noDups_CovPlus.tdf, \*\_Aligned\_noDups_CovMinus.tdf and unique mappers: \*\_uniqueAligned_noDups_CovPlus.tdf \*\_uniqueAligned_noDups_CovMinus.tdf  

