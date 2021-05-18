# Muman mitoribosome profiling analysis

This repository includes the scripts and annotation files needed to analyze mitoribosome profiling data generated from human cells. The directories are listed in order to take raw fastq files through trimming, alignment, quality control, and many other library characteristics. There are README files for each analysis step explaining how scripts are run and what is generated with each. Optimized to run on HMS O2 computing cluster (SLURM job scheduler).

# Analysis steps

0. Create STAR index
1. Trim and align raw reads, remove PCR duplicates  >  get library compositions, RPF length distributions, 5' bedgraphs for viewing on IGV
2. Calculate periodicity on 5' and 3' ends PER RPF LENGTH  >  Needed for accurately determining A-site transformation, in combination with RPF length distibutions (step 1) and Vplots: RPFlength vs. genomic position (step 3)
3. Make bed files for Vplots  >  5' and 3' plus(P) and minus(M) files for input to ScatterPlotLengths.R to visualize read lengths along genes
4. A-site transformation  >  get periodicity, A-site bedgraphs for viewing on IGV
5. Count reads on features using featureCounts  >  get unique- and multi- (all-) aligned readcounts across genes. Use in featureCounts_addGeneName_RPK.R to get RPK values
6. Calculate codon coverage
7. Get codon occupancy

# Data availability and manuscript

Fastq files are deposited in the GEO database under the accession number GSE173283. The link to our full manuscript will be provided here upon publication.

![alt text](https://github.com/mtcouvi/human-mitoribosome-profiling/blob/main/Method.pdf?raw=true)

## 0_CreateSTARindex
Download fasta and gtf files from desired source (e.g. GENCODE) and follow instuctions.txt

## 1_AlignData

Steps for complete analysis. Please see each script for dependencies and details

### 1. ProcessFASTQ_hMitoRP.sh LibName inputFASTQ UMI(e.g. 3p6, 3p6_5p4, 3p10_5p4)
   
#### outputs:
   - ${LibName}_Counts.txt : alignment counts
   - ${LibName}_LengthDist.txt : mito mRNA read lengths (with soft-clipped nt counted in length
   - ${LibName}_Aligned.Mito_mRNA.noDups_samstatsReadlength.noSoft.txt : mito mRNA read lengths with soft-clipped nt not counted in length
   - ${LibName}_Mito_mRNA.noDups.5p.plus.bedGraph : 5' bedgraph for IGV visualization
   
   Use notes at the end of script to combine all samples into one file using unix commands (paste, awk, etc)
   
### 2. CountFramePerLength.sh LibName
   
   This is not required but is good to check for determining ideal Asite transformation, in combination with read length distributions. May also need to do (3) as well
   Use notes at the end of script to combine all samples into one file using unix commands (paste, awk, etc)
   
#### outputs:
   - ${Experiment}_5p_allFrameCounts.txt
   - ${Experiment}_3p_allFrameCounts.txt
   
   --> Plot frame frequency for each footprint length and choose which lengths to use in downstream analysis

### 3. MakeLengthBeds.sh LibName
   
   This is not required unless you want to make Vplots to visualize read profile
   
   #### outputs:
   - ${LibName}_Mito_mRNA_lengths3p_M.bed
   - ${LibName}_Mito_mRNA_lengths3p_P.bed                                                                                                                
   - ${LibName}_Mito_mRNA_lengths5p_M.bed                                                                                                                  
   - ${LibName}_Mito_mRNA_lengths5p_P.bed 
   
   --> Then put these through the script ScatterPlotLengths.R to make Vplots    

### 4. AsiteAndCountFrame.sh LibName Exp sizeRange
   
   This should only be run after determining A-site offset and read length sizes to be used. Calls script SAM2hMitoFBED_softClip_<Exp>.py, which provides experiment-specific offsets
   
   #### outputs:
   - ${Experiment}_Asite_${sizeRange}_FrameCount_ignore1st3.txt
   - ${Experiment}_${sizeRange}_Asite_RPK_inside.txt : Summing of readcounts. This is now replaced by featureCounts (step 5 below)
   - ${LibName}_Mito_mRNA.noDups.Asite_30to33_P/M.bedGraph : Asite bedGraphs
   - ${Experiment}_Mito_mRNA.noDups.Asite_31to33_PAll.bedGraph -> which goes into StackedFramesPlot_assignValues.R for ad hoc scoring
   
### 5. run_featureCounts_hMitoRP.sh LibName offset

### 6. CodonOccupancy.sh LibName sizeRange
To get codon occupancy, move Asite bedGraphs (${LibName}_Mito_mRNA.noDups.Asite_30to33_P/M.bedGraph) to personal computer, make new directory 'CodonFrequency' and run                                                                                                         
