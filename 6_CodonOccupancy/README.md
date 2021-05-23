This script uses A site bedGraphs from step 4 and calculates mean and individual codon occupancy. 

# 1. Run CodonOccupancy.sh
This script calls python scripts that require biopython, numpy and scipy, so first need to set up a python virtual environment to install packages.  
Do this only once for each experiment (directory). It will create a new directory with the name of the virtual environment (hMitoRP).
```bash
module purge
module load gcc/6.2.0
module load python/3.7.4

# Set up virtual environment
virtualenv hMitoRP

# Activate virtual environment
source mrp/bin/activate
unset PYTHONPATH
pip3 install biopython
pip3 install numpy
pip3 install scipy

# Leave virtual environment
deactivate
```

Replace S1, S2, etc with library (sample) names, update number of samples (e.g. for 4 samples use {0..3}), and run code block.
```bash
Libs="S1 S2 S3 S4"
sizeRange="31to33"
scriptpath="./6_CodonOccupancy/Scripts/CodonOccupancy.sh"

for lib in $Libs
do
sbatch -e logs/6_CodonOccupancy_${lib}.err -o logs/6_CodonOccupancy_${lib}.log $scriptpath $lib $sizeRange
done
```

# 2. Combine outputs
**Do not do this until all jobs are complete**  

Be sure to modify paste command based on number of samples you have in your experiment. Only the last sample added should include the fourth column, which is codon frequency.

```bash
cd 6_CodonOccupancy
Experiment="MitoRP1"
Libs=("S1" "S2" "S3" "S4")
sizeRange="31to33"

paste <(awk '{print $1"\t"$2"\t"$3}' ${Libs[0]}_${sizeRange}_CodonDensities.txt) <(awk '{print $2"\t"$3}' ${Libs[1]}_${sizeRange}_CodonDensities.txt) <(awk '{print $2"\t"$3}' ${Libs[2]}_${sizeRange}_CodonDensities.txt) <(awk '{print $2"\t"$3"\t"$4}' ${Libs[3]}_${sizeRange}_CodonDensities.txt) > ../${Experiment}_${sizeRange}_CodonDensities.txt
cd ..
```
  
# 3. Outputs
  - Mean occupancy for each codon and standard deviation across all occurances of that codon: \*_CodonDensities.txt 
  - Individual occupancies for all occurances of each codon: \*_indCodonDensities.txt

# 4. Make plots
Requires the following R packages: Rfast
Run R script CodonDensities_tRNAabundance.R
```bash
module load R/4.0.1
Experiment="MitoRP1"
sizeRange="31to33"
Rscript ./6_CodonOccupancy/Scripts/CodonDensities_tRNAabundance.R $Experiment $sizeRange
```
This produces 3 plot files: 
- Barplot of mean codon occupancies with a bar for each sample
  - \*_CodonDens_barplot.pdf
- Boxplots of individual codon occupancies for Asp, Leu, Phe for each sample
  - \*_indCodonDens_boxplot.pdf
- Scatterplots of mean codon occupancy for each sample vs tRNA abundance from myoblasts (Richter et al, Nat Commun. 2018)
  - \*_CodonDens_vs_tRNAabundance.pdf