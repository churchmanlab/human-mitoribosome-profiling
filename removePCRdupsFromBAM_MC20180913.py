#!/usr/bin/env python
"""

Date : March 18, 2016

Author : Heather Landry

Remove reads in bam file resulting from PCR duplicates. This script records all barcodes and coordinates
at a specific position. For every bam line, if the barcode and coordinate has not been seen previously, it
will print; if the barcode and position has been seen previously, it will not print to a new file.

use : python removePCRdupsFromBAM.py    iBAM (input BAM file with only unique alignments) [1]
                                        oBAM (output BAM file containing only non duplicated reads) [2]
                                        
Update April 25, 2017  - Mary Couvillion:
Skip lines where target id (tid) is <0
                                            
"""
import sys, pysam, os, numpy, re

iBAM = pysam.Samfile(sys.argv[1], 'rb')
oBAM = pysam.Samfile(sys.argv[2], 'wb', template=iBAM)

MB = set()

# read through starting bam file
for read in iBAM:
    mb = read.qname.split('_MolecularBarcode:')[1]
    if read.tid < 0:
        continue
    if read.tid >= 0:
        chrom = iBAM.getrname(read.tid)
        cigar = read.cigarstring
        
    # selecting the 5' position for pos strand
    if not read.is_reverse:
        start = read.reference_start
        std='pos'
    
    # selecting the 5' position for neg strand
    if  read.is_reverse:
        start = read.reference_end
        std='neg'

    key = str(chrom)+"_"+str(start)+"_"+str(std)+"_"+str(mb)+"_"+cigar
    
    # output 1 read per molecular barcode
    if key not in MB:
        MB.add(key)
        oBAM.write(read)

iBAM.close()
oBAM.close()

# module load dev/python/2.7.10
# pip freeze  # to see all modules sub-loaded

