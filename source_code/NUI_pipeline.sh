#### Input to this pipeline:
# Long Ranger alignment BAM file and Supernova pseudohaplotypes in FASTA format (2 pseudohaplotypes per sample)
#
#------------------------------------------- Pipeline descriptions:
# 1. Unmapped and poorly mapped reads are extracted from the alignment BAM file
# 2. This collection of reads are removed if less than 70% of the bases have a quality score of 30 or above
# 3. Align these reads to the two pseudohaplotypes and filter for good alignments
# 4. Extract the corresponding pseudohaplotype sequences if the reads form clusters
# 5. Extend the pseudohaplotype sequences by 7kb on both ends
# 6. Align these sequences back to hg38 reference genome to identify breakpoints
# 7. Compute the breakpoints for each alignmnet
# 8. Filter alignment results to meet NUI definitions
# 
#!/bin/bash 
usage() {
    NAME=$(basename $0)
    cat <<EOF
Usage:
  ${NAME} 
Please edit NUI_config.sh accordingly before each run. The config file, 
along with the R scripts, hg38_primary_header_meaning.txt, segdups.bedpe, 
and sv_blacklist.bed must be placed together with this script.

EOF
}

## load variables
if [[ ! -f ./NUI_config.sh ]]; then
 usage
 echo "Error: Missing configuration file (NUI_config.sh) !"
 exit 1
fi

source ./NUI_config.sh

if [ ! -f $BAMLOC ]; then
    echo "ERROR: Alignment BAM file not found!"
    exit 1
fi 
if [ ! -f $HG38 ]; then
    echo "ERROR: HG38 core reference file not found!"
    exit 1
fi 
if [ ! -f $HG38ALL ]; then
    echo "ERROR: HG38 all reference file not found!"
    exit 1
fi 
if [ ! -d $FASTQDIR ]; then
    echo "ERROR: $FASTQDIR doesn't exist!"
    exit 1
fi 
if [ ! -d $FASTADIR ]; then
    echo "ERROR: $FASTADIR doesn't exist!"
    exit 1
fi 

LOGFILE=./"$SAMPLEID"_run.log

# exit immediately if there's an error in the main pipeline
set -e

## main pipeline
pipeline(){
echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> START: " $0 

score=30
percent=70

cd "$WORKDIR"
    
if [ ! -d "$SAMPLEID" ]; then
    mkdir "$SAMPLEID"
fi
    
cd "$SAMPLEID"

echo [`date +"%Y-%m-%d %H:%M:%S"`] "   * Extracting unaligned reads from $SAMPLEID"
# extract reads that have more than 40bp clipped off
samtools view -h "$BAMLOC" | samblaster --ignoreUnmated --minClipSize 40 -r | \
    fastq_quality_filter -Q33 -q "$score" -p "$percent" -i - | paste - - - - | cut -c2- | \
    cut -f1 -d"_" > "$SAMPLEID"_unaligned_1.txt &

# get the pid to check the progress
pid=$!
    
# extract reads that mapped poorly
(samtools view "$BAMLOC" | \
    awk '$2==81 || $2==161 || $2==97 || $2==145 || $2==65 || $2==129 || $2==113 || $2==177 || $2==67 || $2==131 || $2==115 || $2==177'; \
    sambamba view -t "$CORES" -F "[AS] <= -80" "$BAMLOC") | \
    cat | awk '{print $1}' > "$SAMPLEID"_unaligned_2.txt

# need to make sure both jobs are finished before moving on
# this is necessary if sample_unaligned_2.txt finishes writing before sample_unaligned_1.txt
wait $pid 
    
echo [`date +"%Y-%m-%d %H:%M:%S"`] "   * Sorting unaligned reads for $SAMPLEID"
cat "$SAMPLEID"_unaligned_1.txt "$SAMPLEID"_unaligned_2.txt | sort -u -S "$MEM" -T "$TEMPDIR" > "$SAMPLEID"_unaligned.txt

# remove intermediate files
rm "$SAMPLEID"_unaligned_1.txt
rm "$SAMPLEID"_unaligned_2.txt

# extract fastq based on the read names
echo [`date +"%Y-%m-%d %H:%M:%S"`] "   * Extracting fastq"
eval zcat "$FASTQDIR"/*R* | seqtk subseq - "$SAMPLEID"_unaligned.txt > unaligned.fq 

echo [`date +"%Y-%m-%d %H:%M:%S"`] "   * Filtering and trimming fastq"
fastq_quality_filter -Q33 -q "$score" -p "$percent" -i unaligned.fq | paste - - - - | \
    awk 'length($3)>150' | sort -u -k1,1 -k2,2 -S "$MEM" -T "$TEMPDIR" > temp.fq

# discard reads if not in pairs
awk 'BEGIN { FS=" " } { c[$1]++; l[$1,c[$1]]=$0 } END { for (i in c) { if (c[i] == 2) for (j = 1; j <= c[i]; j++) print l[i,j] } }' temp.fq | \
    awk '{print $1, $3, $4, $5}' - | \
    
    # trim read1.first 16+7bp of R1
    # R1 contains the 16bp 10x barcode + 7bp of low accuracy sequence from an N-mer oligo.
    awk '(NR % 2) {$2 = substr($2,24)}1 !NR % 2 {print}' - | awk '(NR % 2) {$4 = substr($4,24)}1 !NR % 2 {print}' - | \
    tr ' ' '\n' | tr '\t' '\n'> "$SAMPLEID"_final.fq 

rm unaligned.fq
rm temp.fq

# BWA alignment
echo [`date +"%Y-%m-%d %H:%M:%S"`] "   * Starting BWA-MEM"
    
# align the fastq to the hg38 core sequences to spot good alignments and remove them later
bwa mem -t "$CORES" -p "$HG38" "$SAMPLEID"_final.fq > aln-se_hg38.sam
sambamba view -t "$CORES" -F "(sequence_length <= 128 and [AS] >= 90) or (sequence_length >= 151 and [AS] >= 113)" \
    -S aln-se_hg38.sam | awk '($2==99||$2==147||$2==83||$2==163) && $5>=30 {print $1}' | uniq -d \
    > aln-se_hg38_goodAln_ID_pairs.txt

for i in 2.1 2.2; do
    
    # align to the pseudohaplotypes 
    bwa mem -t "$CORES" -p "$FASTADIR"/"$SAMPLEID"_pseudohap"$i".fasta "$SAMPLEID"_final.fq > aln-se_pseudohap"$i".sam
    
    # good alignments spotted earlier are filtered out here
    sambamba view -t "$CORES" -F "(sequence_length <= 128 and [AS] >= 90) or (sequence_length >= 151 and [AS] >= 113)" \
        -S aln-se_pseudohap"$i".sam | awk '($2==99||$2==147||$2==83||$2==163) && $5>=30' | \
        awk 'FNR==NR{a[$1];next};!($1 in a)' aln-se_hg38_goodAln_ID_pairs.txt -  > aln-se_pseudohap"$i"_goodAln.sam

done

# compress the fastq and alignment sam files
echo [`date +"%Y-%m-%d %H:%M:%S"`] "   * Compressing FASTQ and SAM files"
gzip "$SAMPLEID"_final.fq &
samtools view -Sbh aln-se_hg38.sam > aln-se_hg38.bam; rm aln-se_hg38.sam &
samtools view -Sbh aln-se_pseudohap2.1.sam > aln-se_pseudohap2.1.bam; rm aln-se_pseudohap2.1.sam 
samtools view -Sbh aln-se_pseudohap2.2.sam > aln-se_pseudohap2.2.bam; rm aln-se_pseudohap2.2.sam 

# identify read cluster
echo [`date +"%Y-%m-%d %H:%M:%S"`] "   * Calculating coverage"
for i in 2.1 2.2; do

    #Calculate coverage genomewide based on aln-se_pseudohap1.bam and extract regions with coverages between 8-100x
    samtools view -H aln-se_pseudohap"$i".bam | cat - aln-se_pseudohap"$i"_goodAln.sam | samtools view -Sbh - | \
    samtools sort -@ 4 -m 16G -o - - | bedtools genomecov -ibam - -bg | awk '$4>8 && $4<100' | \
    bedtools merge -i - -d 200 > aln-se_pseudohap"$i"_genomecov_new.bed 
    
    #Make a fasta file and extend 7kb on each end
    awk '{$2=$2-7000; $3=$3+7000; print}' aln-se_pseudohap"$i"_genomecov_new.bed | awk '$2<0 {$2=1}1' | \
    awk '{print $1":"$2"-"$3}' > aln-se_pseudohap"$i"_genomecov_new_7000.bed
    
    bedtools getfasta -fi "$FASTADIR"/"$SAMPLEID"_pseudohap"$i".fasta \
        -bed aln-se_pseudohap"$i"_genomecov_new.bed -fo pseudohap"$i"_ROI_new.fa
    
    #Remove any sequence with N inside the ROI and make a fasta file ready for lastz alignment
    #This fasta file should have a 7kb anchor on each side
    cat pseudohap"$i"_ROI_new.fa | paste - - |grep -v N - | awk '{print $1}' | \
    awk -F':|-' '{print $1, $2-7000, $3+7000}' | awk '$2<0 {$2=1}1' | awk '{print $1":"$2"-"$3}' | cut -c2- | \
    grep -Fwf - aln-se_pseudohap"$i"_genomecov_new_7000.bed | awk -F':|-' -v OFS='\t' '{print $1, $2, $3}' | \
    bedtools merge -i - | bedtools getfasta -fi "$FASTADIR"/"$SAMPLEID"_pseudohap"$i".fasta -bed - \
    -fo lastz_ready_pseudohap"$i".fa

done 

if [ ! -d lastz ]; then
    mkdir lastz
fi
    
cd lastz

# NUI anchoring step
echo [`date +"%Y-%m-%d %H:%M:%S"`] "   * Anchoring NUI"
for i in 2.1 2.2; do

    (lastz_32 "$HG38ALL"[multiple][unmask] ../lastz_ready_pseudohap"$i".fa[nameparse=darkspace] \
        --format=general:number,name1,start1,end1,length1,size1,name2,start2+,end2+,length2,size2,strand2,score,nmatch,nmismatch,cigar,identity,cov% \
        --step=20 --seed=match15 --notraNUItion --exact=400 --identity=95 \
        --match=1,5 --ambiguous=iupac > unmask_hg38.p11vcontig_ROI"$i".txt;
    
    # If the query is covered more than 99%, remove them. They are good alignments.
    awk '$19>99 || $19=="100.0%" {print $7}' unmask_hg38.p11vcontig_ROI"$i".txt | \
        grep -vf - unmask_hg38.p11vcontig_ROI"$i".txt | \
        awk 'substr($2,1,3)=="CM0" {print}' > unmask_hg38.p11vcontig_ROI"$i"_CLEAN.txt;
    
    # If an NUI aligns to less than 5 loci, append all alignment results to the output, which will be further processed 
    awk '{print $7}' unmask_hg38.p11vcontig_ROI"$i"_CLEAN.txt | uniq -c | awk '$1<5 {print $2}'| \
        grep -Fwf - unmask_hg38.p11vcontig_ROI"$i"_CLEAN.txt > unmask_hg38.p11vcontig_ROI"$i"_FINAL_7000_merged.txt;
    
    # If an NUI aligns to multiple locations, keep the alignments that are at least 3.5kb long AND fewer than 5 aln results
    awk '{print $7}' unmask_hg38.p11vcontig_ROI"$i"_CLEAN.txt | uniq -c | awk '$1>=5 {print $2}'| \
        grep -Fwf - unmask_hg38.p11vcontig_ROI"$i"_CLEAN.txt | awk '$10 >=3500 {print $7}' | \
        uniq -c | awk '$1<5 {print $2}' | grep -Fwf - unmask_hg38.p11vcontig_ROI"$i"_CLEAN.txt | \
        awk '$10 >=3500 {print}' >> unmask_hg38.p11vcontig_ROI"$i"_FINAL_7000_merged.txt;
     
    # identify insertions 
    # ensure insertions don't overlap with segmental duplication/blacklist
    echo [`date +"%Y-%m-%d %H:%M:%S"`] "   * Identifying NUI from pseudohap$i";
    Rscript "$WORKDIR"/process_NUI_aln.R -s "$SAMPLEID" -v "$i";
     
    # remove insertions if the sequences (including 50 flanking bases on each side) contain more than 10 N's
    awk -F'\t' '{print $10}' SV_final_"$i".txt | awk -F':|-' -v OFS='\t' '{print $1, $2-15, $3+15}' | \
    bedtools getfasta -fi "$FASTADIR"/"$SAMPLEID"_pseudohap"$i".fasta -bed - -fo SV_final_"$i".fa;
    grep -B 1 NNNNNNNNNN SV_final_"$i".fa | grep ">" | cut -c 2- | awk -F':|-' '{print $1":"$2+15"-"$3-15}' |\
    grep -Fvf - SV_final_"$i".txt > SV_final_"$i"_CLEAN.txt;
    
    # remove intermediate files
    rm SV_final_"$i".txt) &
    
done
wait
sleep 1h
    
# Combine the two pseudohap results 
Rscript "$WORKDIR"/combine_SV_calculate_genotype.R 

# Filter for sequences that meet the definitions of an NUI
mkdir ../NUI
cd ../NUI
    
# Extract the FASTA sequences for RepeatMasker and dustmasker
awk -F'\t' '$11=="2.1" && $7=="Insertion" {print $10}' ../lastz/SV_final_combined.txt | \
    awk -F':|-' '{print $1, $2, sprintf("%9d", $3)}' | awk -v OFS='\t' '{print $1, $2, $3}' | \
    bedtools getfasta -fi "$FASTADIR"/"$SAMPLEID"_pseudohap2.1.fasta -bed /dev/stdin -fo temp2.1.fa 

awk -F'\t' '$11=="2.2" && $7=="Insertion" {print $10}' ../lastz/SV_final_combined.txt | \
    awk -F':|-' '{print $1, $2, sprintf("%9d", $3)}' | awk -v OFS='\t' '{print $1, $2, $3}' | \
    bedtools getfasta -fi "$FASTADIR"/"$SAMPLEID"_pseudohap2.2.fasta -bed /dev/stdin -fo temp2.2.fa 

cat temp2.1.fa temp2.2.fa | paste - - | awk 'length($2)>=50 {print $0}'|tr '\t' '\n' > NUI.fa
rm temp2.1.fa temp2.2.fa

# Mask the NUIs
RepeatMasker --species human -pa "$CORES" NUI.fa 

# Remove sequences with low complexity using dustmasker
# Dustmasker is more stringent that repeatMasker
dustmasker -outfmt fasta -in NUI.fa.masked -out NUI.dustmasked.fa

# Check which sequence has less than 50 unique bases
# Unique meaning non-masked
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < NUI.dustmasked.fa | \
    grep .... | paste - - | awk '{print gsub(/A|T|G|C/, "") "\t" $0}' | awk '$1<50 {print $2}' | \
    cut -c 2- > repeatMasker_discard_ID.txt 

# Make sure NUIs don't align well with something already in the reference
# Add 50bp to each end
awk -F'\t' '$11=="2.1" && ($7=="Insertion") {print $10}' ../lastz/SV_final_combined.txt | \
    awk -F':|-' '{print $1, $2, sprintf("%9d", $3)}' | awk -v OFS='\t' '{print $1, $2-50, $3+50}' | \
    bedtools getfasta -fi "$FASTADIR"/"$SAMPLEID"_pseudohap2.1.fasta -bed /dev/stdin -fo temp2.1.fa 
awk -F'\t' '$11=="2.2" && ($7=="Insertion") {print $10}' ../lastz/SV_final_combined.txt | \
    awk -F':|-' '{print $1, $2, sprintf("%9d", $3)}' | awk -v OFS='\t' '{print $1, $2-50, $3+50}' | \
    bedtools getfasta -fi "$FASTADIR"/"$SAMPLEID"_pseudohap2.2.fasta -bed /dev/stdin -fo temp2.2.fa 

cat temp2.1.fa temp2.2.fa | paste - - | awk 'length($2)>=50 {print $0}'|tr '\t' '\n' > NUI_50.fa
rm temp2.1.fa temp2.2.fa

# Use BLAST to realign all the sequences (+50bp flanking sequences) to confirm they don't align to hg38
blastn -query NUI_50.fa \
        -db /media/KwokRaid02/karen/reference/hg38/hg38.p11.fa \
        -task megablast \
        -dust no \
        -outfmt "7 qseqid sseqid evalue bitscore qlen pident length salltitles qstart qend sstart send nident mismatch gapopen gaps qcovs qcovhsp" \
        -max_target_seqs 1 \
        -max_hsps 1 \
        -out NUI_50_output.blast \
        -num_threads 16

awk -F'\t' '$6>=95 && $17==100 {print $1}' NUI_50_output.blast | sort -u | \
    awk -F':|-' '{print $1":"$2+50"-"$3-50}' > blast_discard_ID.txt

cat blast_discard_ID.txt repeatMasker_discard_ID.txt | sort -u | grep -Fvf - ../lastz/SV_final_combined.txt | \
    awk -F'\t' '($7=="Insertion") && ($9>=50) {print $0}' > SV_final_combined_NR_FILTERED.txt
    
} #end of pipeline

pipeline 2>&1 | tee $LOGFILE

