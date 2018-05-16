# NUI_pipeline

## Softwares used in the pipeline:
samtools v1.2 <br>
sambamba v0.5.9 <br>
bedtools v2.17.0 <br>
seqtk subseq v1.0 <br>
FASTX Toolkit v0.0.14 <br>
samblaster v0.1.24 (Source code was modified to redirect the fastq output to stdout) <br>
bwa v0.7.15

**R packages:** <br>
GenomicRanges v1.26.4 <br>
stringr v1.2.0 <br>

## Pipeline overview:
![pipeline](https://user-images.githubusercontent.com/22200237/40133546-f892daa8-58f4-11e8-9f16-2b68354be019.jpg)
* An asterisk indicates that the steps are repeated twice, one for each pseudo-haplotype.
