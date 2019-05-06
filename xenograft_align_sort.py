"""
Author: xdo
Aim: A simple Snakemake workflow to map paired-end RNAseq with STAR and then run Disambiguate to take out xenograft mouse reads
Pre-Run:    
            source activate py3.5
            ./dl_JHU_synapse.py (in 01/scripts)
Run: snakemake -s xenograft_align_sort.py --cores 8  

Latest modification: 
    - 
"""
import sys
import subprocess 
import re
import os

import synapseclient
import synapseutils
from synapseclient import Project, File, Folder, Activity
from synapseclient import Schema, Column, Table, Row, RowSet, as_table_columns


##-----------------------------------------------##
## A set of functions
##-----------------------------------------------##
def message(mes):
  sys.stderr.write("\n " + str(mes) + "\n")

##-----------------------------------------------##
## set wd 
##-----------------------------------------------##
BASE_DIR = "/home/xdoan/MPNST/"

##-----------------------------------------------##
## declare variables
##-----------------------------------------------##
REF_DIR = BASE_DIR + "reference/"
hg38_FA = BASE_DIR + "reference/hg38/hg38.fa" 
hg38_GTF = BASE_DIR + "reference/hg38/hg38.gtf" 
mm10_FA = BASE_DIR + "reference/mm10/mm10.fa" 
mm10_GTF = BASE_DIR + "reference/mm10/mm10.gtf" 

# INDEX = BASE_DIR + "/index/hg38.fa" #hg38 index 
# REF_BED = "/home/xdo/Data/bwa/reference/hg38.ALR.bed"
# MKV_OUT = ["tis", "ois", "suf", "bwt", "bck", "lcp", "skp"]
# GENOME = "/home/xdo/Data/bwa/wgs/Esoph/eso-1060/hg38.genome"

##-----------------------------------------------##
## The list of samples to be processed
##-----------------------------------------------##
WC = glob_wildcards(BASE_DIR + "01_dl_files/output/{smp}_1.fastq.gz")
SAMPLES = set(WC.smp)
# for smp in SAMPLES:
#   message(smp)

WC2 = glob_wildcards(BASE_DIR + "01_dl_files/output/{ind}_CCDYYANXX_{etc}.fastq.gz")
INDIVIDUALS = set(WC2.ind)
for ind in INDIVIDUALS:
  message(ind)
##-----------------------------------------------##
## all rule
##-----------------------------------------------##
rule all:
  input: 
    BASE_DIR + "reference/mm10/chrName.txt",
    BASE_DIR + "reference/hg38/chrName.txt",
    expand(BASE_DIR + "02_run_STAR/output/hg38/{smp}_Aligned.sortedByCoord.out.bam", smp=SAMPLES),
    expand(BASE_DIR + "03_merge/output/hg38/{ind}_Aligned.sortedByCoord.out.bam", ind=INDIVIDUALS),
    expand(BASE_DIR + "02_run_STAR/output/mm10/{smp}_Aligned.sortedByCoord.out.bam", smp=SAMPLES),
    expand(BASE_DIR + "03_merge/output/mm10/{ind}_Aligned.sortedByCoord.out.bam", ind=INDIVIDUALS),
    expand(BASE_DIR + "04_disambiguate/{ind}.disambiguatedSpeciesA.bam", ind = INDIVIDUALS)

## STAR alignment
##-----------------------------------------------##
rule starIndex_hg38:
  input: 
    ref = hg38_FA,
    gtf = hg38_GTF,
    starref = REF_DIR +"hg38/", 
    # tmpdir = TMP_DIR
  output: 
    BASE_DIR + "reference/hg38/chrName.txt"
  threads: 8
  shell: 
    "STAR --runThreadN {threads} --runMode genomeGenerate --genomeDir {input.starref} --genomeFastaFiles {input.ref} --sjdbGTFfile {input.gtf} --sjdbOverhang 100 --genomeSAindexNbases 4 --limitGenomeGenerateRAM 30000000000 --genomeSAsparseD 2"# --outTmpDir {input.tmpdir}"

rule starIndex_mm10:
  input:
    ref=mm10_FA,
    gtf=mm10_GTF,
    starref=REF_DIR +"mm10/",
    # tmpdir = TMP_DIR 
  output:
    BASE_DIR + "reference/mm10/chrName.txt"
  threads: 8
  shell: 
    "STAR --runThreadN {threads} --runMode genomeGenerate --genomeDir {input.starref} --genomeFastaFiles {input.ref} --sjdbGTFfile {input.gtf} --sjdbOverhang 100 --genomeSAindexNbases 4 --limitGenomeGenerateRAM 30000000000 --genomeSAsparseD 2"# --outTmpDir {input.tmpdir}"

rule star_pe_hg38:
    input:
        fq1 = BASE_DIR + "01_dl_files/output/{smp}_1.fastq.gz",
        fq2 = BASE_DIR + "01_dl_files/output/{smp}_2.fastq.gz" 
    output:
        BASE_DIR + "02_run_STAR/output/hg38/{smp}_Aligned.sortedByCoord.out.bam"
    log:
        BASE_DIR + "00logs/star/hg38/{smp}.log"
    params:
      starref = REF_DIR +"hg38",
      smp = BASE_DIR + "02_run_STAR/output/hg38/{smp}_", 
      logdir = BASE_DIR + "00logs/star/hg38/"
    threads: 8
    shell:
      """ STAR --runThreadN {threads} --genomeDir {params.starref} \
      --readFilesIn {input.fq1} {input.fq2} --outSAMtype BAM SortedByCoordinate \
      --outFileNamePrefix {params.smp} --outStd Log {log} --readFilesCommand zcat """

# Function for retrieving list of files with same ind
def bamlist(individual):
    paths = []
    ind = str(individual)
    for root, dirs, files in os.walk(r"/home/xdoan/MPNST/02_run_STAR/output/hg38", topdown=True):
        for name in files:
            if name.startswith(ind) & name.endswith(".bam"):
                path = os.path.join(root, name)
                paths.append(path)
    return(paths)

rule merge_hg38_star_output:
    input:
      # i = BASE_DIR + "02_run_STAR/output/hg38/{smp}_Aligned.sortedByCoord.out.bam",
      bamlist = bamlist
    output:
      BASE_DIR + "03_merge/output/hg38/{ind}_Aligned.sortedByCoord.out.bam"
    shell:
      "samtools merge {output} {input.bamlist} "

rule star_pe_mm10:
    input:
        fq1 = BASE_DIR + "01_dl_files/output/{smp}_1.fastq.gz",
        fq2 = BASE_DIR + "01_dl_files/output/{smp}_2.fastq.gz" 
    output:
        BASE_DIR + "02_run_STAR/output/mm10/{smp}_Aligned.sortedByCoord.out.bam"
    log:
        BASE_DIR + "00logs/star/mm10/{smp}.log"
    params:
      starref = REF_DIR +"mm10",
      smp = BASE_DIR + "02_run_STAR/output/mm10/{smp}_", 
      logdir = BASE_DIR + "00logs/star/mm10/"
    threads: 8
    shell:
      """ STAR --runThreadN {threads} --genomeDir {params.starref} \
      --readFilesIn {input.fq1} {input.fq2} --outSAMtype BAM SortedByCoordinate\
      --outFileNamePrefix {params.smp} --outStd Log {log} --readFilesCommand zcat """
  
rule merge_mm10_star_output:
    input:
      # i = BASE_DIR + "02_run_STAR/output/hg38/{smp}_Aligned.sortedByCoord.out.bam",
      bamlist = bamlist
    output:
      BASE_DIR + "03_merge/output/mm10/{ind}_Aligned.sortedByCoord.out.bam"
    shell:
      "samtools merge {output} {input.bamlist} "

rule disambiguate:
    input:
        human= BASE_DIR + "03_merge/output/hg38/{ind}_Aligned.sortedByCoord.out.bam",
        mouse= BASE_DIR + "03_merge/output/mm10/{ind}_Aligned.sortedByCoord.out.bam",
    output:
        BASE_DIR + "04_disambiguate/{ind}.disambiguatedSpeciesA.bam"
    params:
        prefix = "{ind}",
        outdir = BASE_DIR + "04_disambiguate/"
    shell:
        "ngs_disambiguate {input.human} {input.mouse} -s {params.prefix} -a star -o {params.outdir}"
      