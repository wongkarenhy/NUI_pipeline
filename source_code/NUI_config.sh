# Configuration file for the NUI pipeline
# Place this script together with NUI_pipeline.sh 

#number of CPUs to be used
CORES=32

#amount of memory to be used (mostly for sorting)
MEM=80G

#sample to be analyzed
SAMPLEID="HG00250"

#### File paths for input data
BASEDIR="/path"
WORKDIR="$BASEDIR/NUI" 

#temp directory for sorting
TEMPDIR="$BASEDIR/temp"

#FASTA files corresponding to the two pseudohap generated with 10xG Supernova
#FASTA file names should be "$SAMPLEID"_pseudohap2.1.fasta and "$SAMPLEID"_pseudohap2.2.fasta
FASTADIR="$BASEDIR/reference/10X_pseudohap"

#BAM file corresponding to the 10xG Long Ranger alignment output
BAMLOC="/path/to/HG00250_longranger/outs/phased_possorted_bam.bam"

#FASTQ directory should contain gzipped FASTQ files
#FASTQ filenames format used in the pipeline: *R*.gz
FASTQDIR="/path/to/HG00250_mkfastq/H5VHGALXX_34/HG00250"

#path to the core hg38 reference FASTA
HG38="/path/to/reference/hg38/hg38_core_chrs.fa"

#path to the hg38.p11 reference FASTA
HG38ALL="/path/to/reference/hg38/hg38.p11.fa"
