Uses *_Aligned.Mito_mRNA.noDups.bam file from 1_AlignData  

# 1. Run CountFramePerLength.sh
Replace S1, S2, etc with library (sample) names, open script with text editor to modify RPF lengths analyzed (will work for any sizes 17nt-38nt), and run code block
```bash
Libs="S1 S2 S3 S4"
scriptpath="./3_CountFramePerLength/Scripts/CountFramePerLength.sh"

for lib in $Libs
do
sbatch -e logs/3_CountFramePerLength_${lib}.err -o logs/3_CountFramePerLength_${lib}.log $scriptpath $lib
done
```

# 2. Combine outputs
**Do not do this until all jobs are complete**  
Be sure to modify paste command based on number of samples you have in your experiment

```bash
Experiment="MitoRP1"
Libs=("S1" "S2" "S3" "S4")

ends="5 3"
for end in $ends
do
paste <(awk '{print $1"\t"$2}' ${Libs[0]}_${end}p_FrameCounts.txt) <(awk '{print $2}' ${Libs[1]}_${end}p_FrameCounts.txt) <(awk '{print $2}' ${Libs[2]}_${end}p_FrameCounts.txt) <(awk '{print $2}' ${Libs[3]}_${end}p_FrameCounts.txt) > ${Experiment}_${end}p_allFrameCounts.txt  
done
```

# 3. Outputs
  - *_5p_allFrameCounts.txt
  - *_3p_allFrameCounts.txt  
  
  Paste into excel or equivalent to make a quick bar plot and see how periodicity looks for each size
