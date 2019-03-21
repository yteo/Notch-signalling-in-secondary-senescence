#!/usr/bin/env python
import argparse,subprocess,os,sys
from datetime import datetime
from subprocess import call
import shlex
from subprocess import Popen
from subprocess import PIPE
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import csv
import timeit
from Bio import SeqIO

parser = argparse.ArgumentParser(description='Single cell analysis -- Trimming, HISAT2, HTSeq and counting')
parser.add_argument('-o',action='store',help='Output directory')
parser.add_argument('-t',action='store',help='Number of threads [1]',type=int,default=1)
parser.add_argument('-infile1',action='store',help='Path to gzipped fastq 1')
parser.add_argument('-infile2',action='store',help='Path to gzipped fastq 2')
parser.add_argument('-outputprefix',action='store',help='Output prefix')
#parser.add_argument('-inputtype',action='store',help='Input file type: trimmed or untrimmed')
parser.add_argument('-hisat2Idx',action='store',help='HISAT2 index',default='/data/hg19_scRNA')
parser.add_argument('-htseqann',action='store',help='htseq annotation',default='/data/hg19_ERCC.annotation.gtf ')
args = parser.parse_args()


######################################################################
# Defining arguments
OUTPUTDIR = args.o
N = args.t
trimmedfile1 = args.infile1
trimmedfile2 = args.infile2
OUTPUTPREFIX = args.outputprefix
#NONHUMANLIST = args.nonhumanlist
hisat2Idx = args.hisat2Idx
#inputtype = args.inputtype
htseqanno = args.htseqann

#######################################################################
# Check if output directory exists
def ensure_dir(OUTPUTDIR):
    if not os.path.exists(OUTPUTDIR):
        os.makedirs(OUTPUTDIR)


#######################################################################
# HISAT2 alignment of trimmed fastq files
print datetime.now(), "..... Aligning using HISAT2"
hisat2_align= subprocess.Popen(("hisat2 -q -x  %s -1 %s -2 %s -S %s -p %s --dta-cufflinks" % (hisat2Idx,trimmedfile1 , trimmedfile2, OUTPUTDIR+"/"+OUTPUTPREFIX+".sam", N)), shell = True)
hisat2_align.communicate()

# sort sam files
call("samtools view -bS %s|samtools sort -n - -o %s" %(OUTPUTDIR+"/"+OUTPUTPREFIX+".sam", OUTPUTDIR+"/"+OUTPUTPREFIX+"_sorted"))
#sorting=subprocess.Popen("samtools view -bS %s |samtools sort -n - -o %s" %(OUTPUTDIR+"/"+OUTPUTPREFIX+".sam", OUTPUTDIR+"/"+OUTPUTPREFIX+"_sorted"))
#sorting.communicate()
# remove unsorted sam files
os.remove(OUTPUTDIR+"/"+OUTPUTPREFIX+".sam")


#######################################################################
# HTSeq count
print datetime.now(), "..... Running HTSeq-Count"
htseq=subprocess.Popen("python htseq-count -f bam --stranded=no %s -r name %s > %s" %(OUTPUTDIR+"/"+OUTPUTPREFIX+"_sorted.bam", htseqanno, OUTPUTDIR+"/"+OUTPUTPREFIX+"_htseq.txt"))

#######################################################################
# Counting Hg19 mapped reads or ERCC mapped reads
print datetime.now(), "..... Counting mapped reads"
countERCCreads=subprocess.Popen("samtools view -f2 -F256 %s |grep ERCC- |wc -l" %(OUTPUTDIR+"/"+OUTPUTPREFIX+"_sorted.bam"))
outERCCR=countERCCreads.communicate()
countHg19reads=subprocess.Popen("samtools view -f2 -F256 %s |grep -v ERCC- |wc -l" %(OUTPUTDIR+"/"+OUTPUTPREFIX+"_sorted.bam"))
outHg19R=countHg19reads.communicate()

# out put total mapped reads, Hg19 mapped reads and ERCC mapped reads to a txt file
fmap = open(OUTPUTDIR+"/mappingreadsstat.txt",'a')
fmap.write(OUTPUTPREFIX+"\t"+str(float(outHg19R)/2+float(outERCCR)/2)+"\t"+str(float(outHg19R)/2)+"\t"+str(float(outERCCR)/2)+"\n")
fmap.close()

