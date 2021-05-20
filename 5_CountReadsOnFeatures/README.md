This script runs Rsubread's featureCounts on mito mRNA mappers. Rsubread must be installed in your R environment. The gtf file used has the first 6 codons, last 5 codons, and all overlapping regions removed from CDS annotations, which is the feature used here by featureCounts.  
All read sizes are used here and the A site is approximated using a single offset. See comment in code block below for offset guidelines.  
Separate files are created to report unique mappers and all mappers ('multi'). Because of nuclear mitochondrial pseudogenes it is advised to use all ('multi') mappers. 

# 1. Run run_featureCounts.sh
Replace S1, S2, etc with library (sample) names, update offset and number of samples (e.g. for 4 samples use {0..3}), and run code block.
```bash
Libs=("S1" "S2" "S3" "S4")
offset=14  # 14 for RPF distribution peak at 31; 15 for peak at 32, 33, 34; 16 for peak at 35,36 
scriptpath="./5_CountReadsOnFeatures/Scripts/run_featureCounts.sh"

for i in {0..3}
do
sbatch -e logs/5_CountReadsOnFeatures_${Libs[i]}.err -o logs/5_CountReadsOnFeatures_${Libs[i]}.log $scriptpath ${Libs[i]} $offset
done
```

# 2. Combine outputs
**Do not do this until all jobs are complete**  

Feature length is also added here for calculating RPK in the next step.  
Be sure to modify paste command based on number of samples you have in your experiment but do not delete the final file ${Libs[0]}\_mito\_featureCounts\_CDSignore1stlast\_Length.txt

```bash
Experiment="MitoRP1"
Libs=("S1" "S2" "S3" "S4")

sizeRange='allSizes'

list='CDSignore1stlast_multi CDSignore1stlast_unique'

for type in $list
do
paste <(awk '{print $1"\t"$2}' ${Libs[0]}_mito_featureCounts_${sizeRange}Asite_${type}_noDups.txt) <(awk '{print $2}' ${Libs[1]}_mito_featureCounts_${sizeRange}Asite_${type}_noDups.txt) <(awk '{print $2}' ${Libs[2]}_mito_featureCounts_${sizeRange}Asite_${type}_noDups.txt) <(awk '{print $2}' ${Libs[3]}_mito_featureCounts_${sizeRange}Asite_${type}_noDups.txt) <(awk '{print $2}' ${Libs[0]}_mito_featureCounts_CDSignore1stlast_Length.txt) > ${Experiment}_featureCounts_${sizeRange}Asite_${type}_noDups.txt
done

```
  
# 3. Outputs
  - Unique read counts: ${Experiment}\_featureCounts\_allSizesAsite\_CDSignore1stlast\_unique\_noDups.txt
  - Multi read counts (meaning all reads including multi-mappers: ${Experiment}\_featureCounts\_allSizesAsite\_CDSignore1stlast\_multi\_noDups.txt

# 4. Add gene names and calculate RPK
Run R script AddGeneName_RPK.R
```bash
Experiment="MitoRP1"
Rscript ./5_CountReadsOnFeatures/Scripts/AddGeneName_RPK.R $Experiment
```
