This script runs Rsubread's featureCounts on rRNA-filtered mappers. Rsubread must be installed in your R environment. The gtf file needs to be unzipped first (do this only once: gunzip ./Annotations/GRCh38\_ncRNAs\_mM17\_merge_noMitoOverlap.gtf.gz).  
Here, overlaps are allowed in order to catch some special nuclear OXPHOS cases, but the overlapping regions are removed from mitochondrial annotations in order to not over-count those.
Separate files are created to report unique mappers and all mappers ('multi'). Because of pseudogenes it is advised to use all ('multi') mappers. 

# 1. Run run_featureCounts_hRNAseq.sh
Replace S1, S2, etc with library (sample) names, update number of samples (e.g. for 4 samples use {0..3}), and run code block.
```bash
Libs="S1 S2 S3 S4"
scriptpath="./5_CountReadsOnFeatures_RNAseq/Scripts/run_featureCounts_hRNAseq.sh"

for lib in $Libs
do
sbatch -e logs/5_CountReadsOnFeatures_RNAseq_$lib.err -o logs/5_CountReadsOnFeatures_RNAseq_$lib.log $scriptpath $lib $offset
done
```

# 2. Combine outputs
**Do not do this until all jobs are complete**  

Feature length is also added here for calculating RPK in the next step.  
Be sure to modify paste command based on number of samples you have in your experiment but do not delete the final file ${LibName}\_featureCounts\_CDS\_Length.txt

```bash
Experiment="RNAseq1"
Libs=("S1" "S2" "S3" "S4")

list='CDS_multi CDS_unique'

cd 5_CountReadsOnFeatures_RNAseq

for type in $list
do
paste <(awk '{print $1"\t"$2}' ${Libs[0]}_featureCounts_${type}_noDups.txt) <(awk '{print $2}' ${Libs[1]}_featureCounts_${type}_noDups.txt) <(awk '{print $2}' ${Libs[2]}_featureCounts_${type}_noDups.txt) <(awk '{print $2}' ${Libs[3]}_featureCounts_${type}_noDups.txt) <(awk '{print $2}' ${Libs[0]}_featureCounts_CDS_Length.txt) > ../${Experiment}_featureCounts_${type}_noDups.txt
done
cd ..
```
  
# 3. Outputs
  - Unique read counts: ${Experiment}\_featureCounts\_CDS\_unique\_noDups.txt
  - Multi read counts (meaning all reads including multi-mappers: ${Experiment}\_featureCounts\_CDS\_multi\_noDups.txt

# 4. Add gene names and calculate RPK
Requires the following R packages: data.table, rlist
Run: 
```bash
module load R/4.0.1
Experiment="RNAseq1"
type="multi" # multi unique
geneSets="all mito"
for geneSet in $geneSets
do
Rscript ./5_CountReadsOnFeatures_RNAseq/Scripts/AddGeneName_RPK_hRNAseq.R $Experiment $type $geneSet
done
```
This produces 4 new files: all features: \*\_all\_readCount.txt and \*\_all\_RPK.txt mito genes only: \*\_mito\_readCount.txt and \*\_mito\_RPK.txt