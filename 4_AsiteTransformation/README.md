This script produces A site bedGraphs for specified size range and offsets, as well as periodicity and coverage stats. 

# 1. Run AsiteAndCountFrame.sh
Replace S1, S2, etc with library (sample) names, choose size range based on outputs from steps (1), (2), and (3), and run code block
```bash
Exp="MitoRP1"
Libs="S1 S2 S3 S4"
scriptpath="./4_AsiteTransformation/Scripts/AsiteAndCountFrame.sh"
sizerange="31to33"
offsets="14,15,15"

for lib in ${Libs}
do
sbatch -e logs/4_AsiteTransform_${lib}.err -o logs/4_AsiteTransform_${lib}.log $scriptpath $lib $sizerange $offsets
done
```

# 2. Combine outputs
**Do not do this until all jobs are complete**

```bash
Experiment="MitoRP1"
Libs=("S1" "S2" "S3" "S4")
sizeRange="31to33"
```
- A site frame count (periodicity)
  - Be sure to modify paste command based on number of samples you have in your experiment
  ```bash
  paste <(awk '{print $2"\t"$3}' ${Libs[0]}_Asite_${sizeRange}_FrameCount_ignore1st6.txt) <(awk '{print $3}' ${Libs[1]}_Asite_${sizeRange}_FrameCount_ignore1st6.txt) <(awk '{print $3}' ${Libs[2]}_Asite_${sizeRange}_FrameCount_ignore1st6.txt) <(awk '{print $3}' ${Libs[3]}_Asite_${sizeRange}_FrameCount_ignore1st6.txt) > ${Experiment}_Asite_${sizeRange}_FrameCount_ignore1st6.txt
  ```

- Codon coverage  
  - Be sure to modify paste command based on number of samples you have in your experiment
  ```bash
  paste <(awk '{print $1"\t"$2}' ${Libs[0]}_Asite_${sizeRange}_CodonCoverage_ignore1st6.txt) <(awk '{print $2}' ${Libs[1]}_Asite_${sizeRange}_CodonCoverage_ignore1st6.txt) <(awk '{print $2}' ${Libs[2]}_Asite_${sizeRange}_CodonCoverage_ignore1st6.txt) <(awk '{print $2}' ${Libs[3]}_Asite_${sizeRange}_CodonCoverage_ignore1st6.txt) > ${Experiment}_Asite_${sizeRange}_CodonCoverage_ignore1st6.txt
  ```
  
- Combine all bedGraph files into one summed file (opt)
  - If samples are replicates and you want to increase coverage
  - Be sure to modify paste AND piped awk commands based on number of samples you have in your experiment
  ```bash
  list="P M"
  set="Mito_mRNA.noDups.Asite"
  for strand in $list
  do
  paste <(awk 'NR>1 {print $1"\t"$2"\t"$3"\t"$4}' ${Libs[0]}_${set}_${sizeRange}_${strand}All.txt) <(awk 'NR>1 {print $4}' ${Libs[1]}_${set}_${sizeRange}_${strand}All.txt) <(awk 'NR>1 {print $4}' ${Libs[2]}_${set}_${sizeRange}_${strand}All.txt) <(awk 'NR>1 {print $4}' ${Libs[3]}_${set}_${sizeRange}_${strand}All.txt) | awk '{print $1"\t"$2"\t"$3"\t"$4 + $5 + $6 + $7}' > ${Experiment}_${set}_${sizeRange}_${strand}All.bedGraph
  sed -i '1i\track type=bedGraph color=0,100,50' ${Experiment}_${set}_${sizeRange}_${strand}All.bedGraph
  done
  ```

# 3. Outputs
  - A site bedGraphs: \*\_Mito\_mRNA.noDups.Asite\_${sizeRange}\_P.bedGraph, ${LibName}\_Mito\_mRNA.noDups.Asite\_${sizeRange}\_M.bedGraph
  - A site bedGraphs with all positions filled, for optional combining: ${LibName}\_Mito\_mRNA.noDups.Asite\_${sizeRange}\_PAll.txt, ${LibName}\_Mito\_mRNA.noDups.Asite\_${sizeRange}\_MAll.txt
  - A site periodicity: ${Experiment}\_Asite\_${sizeRange}\_FrameCount\_ignore1st6.txt
  - Codon coverage: ${Experiment}\_Asite\_${sizeRange}\_CodonCoverage\_ignore1st6.txt
