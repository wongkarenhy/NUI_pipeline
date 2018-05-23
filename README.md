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

## System requirements
16-core (or greater) Intel processor <br>
At least 128GB RAM <br>

## Pipeline overview:
![pipeline](https://user-images.githubusercontent.com/22200237/40133546-f892daa8-58f4-11e8-9f16-2b68354be019.jpg)
* An asterisk indicates that the steps are repeated twice, one for each pseudo-haplotype

## Notes:
1. Place all the source codes, along with hg38_primary_header_meaning.txt, segdups.bedpe, and sv_blacklist.bed into the same directory (segdups.bedpe and sv_blacklist.bed are provided by 10x Genomics) <br>
2. NUI_config.sh must be edited manually each time <br>
3. After all samples are processed individually, results can be merged together to generate a unified, non-redundant list of NUIs by running combine_metaNUI.R <br> 
4. Translocated sequences are removed at this point
