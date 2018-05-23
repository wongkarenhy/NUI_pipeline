# This script identifies insertions and filter out NUIs that are in the blacklist/segmental duplication regions

#!/usr/lib64/R/bin/Rscript
######################################################################################################################
#command line options
library(optparse)

option_list = list(
  
  make_option(c("-s", "--sample"), action="store", default=NA, type='character',
              help="sample ID"),
  make_option(c("-v", "--version"), action="store", default=NA, type='character',
              help="pseudohap version either 2.1 or 2.2")

)
opt = parse_args(OptionParser(option_list=option_list))
################################################################################################################

suppressMessages(library(stringr))
suppressMessages(library(GenomicRanges))

# read files
lastz = read.table(paste0("unmask_hg38.p11vcontig_ROI", opt$version,"_FINAL_7000_merged.txt"), stringsAsFactors = F)
colnames(lastz) = c("number","name1","start1","end1","length1","size1","name2","start2+","end2+","length2","size2","strand2","score","nmatch","nmismatch","cigar","identity","idPct","cov%")
definition = read.table("../../hg38_primary_header_meaning.txt", stringsAsFactors = F)
blackList = read.table("../../sv_blacklist.bed",stringsAsFactors = F)[,1:3]
colnames(blackList) = c("chr", "start", "end")
blackList = blackList[nchar(blackList$chr)<=5,]
blackList$chr = substr(blackList$chr,4,5)
segDup = read.table("../../segdups.bedpe", stringsAsFactors = F)[,1:3]
colnames(segDup) = c("chr", "start", "end")
segDup = segDup[nchar(segDup$chr)<=5,]
segDup$chr = substr(segDup$chr,4,5)

#Change the chr codes to chr numbers. 
chr = sapply(as.character(lastz$name1), function(x) definition[,2][match(x, definition[,1], nomatch = NA)]) 
chr = substring(chr, 4,nchar(as.character(chr)))
lastz = cbind(chr, lastz[,3:(ncol(lastz))])

lastz$chr = as.character(lastz$chr)
lastz$name2 = as.character(lastz$name2)

#Rearrange the lastz file
lastz = lastz[order(lastz$name2, lastz$start2),]

#########################################################################################################################################################
#Actual work flow
name = as.character(unique(lastz$name2))

DF_analysis = NULL
for (i in 1:length(name)){
  df = NULL
  df = lastz[lastz$name2==name[i],]
  
  # lots of special cases that we need to consider
  # if df only has 1 line (anchor only on one side), discard
  # if df aligns to >= 3 chromosomes, discard
  if ((nrow(df) == 1) | (length(unique(df$chr)) >= 3)) next
  else{
    
    #Ranges: if one range completely overlaps another, discard the inner alignment
    discard = NULL
    if (nrow(df) > 2){
      
      firstPos = which(df$start2==min(df$start2))
      lastPos = which(df$end2==max(df$end2))
      
      #If df has no alignment that starts at 1 or ends at the last position, discard
      if ((length(firstPos)==0) & (length(lastPos)==0)) {
        discard = df[1:nrow(df),7:8]
      }

      if ((length(firstPos)==1) & (length(lastPos)==1)){
        left_anchor = df[firstPos,7:8]
        right_anchor = df[lastPos, 7:8]
        test_ranges = df[-c(firstPos,lastPos),7:8]
        
        #define the ranges using IRanges function
        rangesA = IRanges(left_anchor$start2, left_anchor$end2)
        rangesB = IRanges(right_anchor$start2, right_anchor$end2)
        for (r in 1:nrow(test_ranges)){
          rangesC = IRanges(test_ranges$start2[r], test_ranges$end2[r])
          ov = countOverlaps(rangesC, c(rangesA,rangesB), type="within")
          if (ov >= 1){
            discard = rbind(discard,test_ranges[r,])
          }
        }
      }
      #If two sequences have the same starting position, choose the one that is the longest
      else if ((length(firstPos)>1) & nrow(df)==3){
        left_anchor = df[firstPos[1],7:8]
        test_anchor = df[firstPos[2], 7:8]
        if (left_anchor$end2 > test_anchor$end2){
          discard = test_anchor
        } else discard = left_anchor
      }
      #If two sequences have the same ending position, choose the one that is the longest
      else if ((length(lastPos)>1) & nrow(df)==3){
        right_anchor = df[lastPos[1], 7:8]
        test_anchor = df[lastPos[2],7:8]
        if (right_anchor$start2 > test_anchor$start2){
          discard = right_anchor
        } else discard = test_anchor
      }
    }
    
    #Remove the discard alignments from df
    if (!is.null(discard)){
      df = df[which(!(df$start2 %in% discard$start2 & df$end2 %in% discard$end2)),]
    }
    
    #Move on to the next one if df is now an empty dataframe
    if (nrow(df)==0 | nrow(df)==4) next
    
    #If one alignment is + and the other is -
    if (length(unique(df$strand2))==2 & nrow(df)==2){
        next
    } else if (length(unique(df$strand2))==1 & nrow(df)==2){ #Most of the alignment should be in this category
        DF_analysis = rbind.data.frame(DF_analysis, df)
    } else if ((nrow(df)==3 & length(unique(df$chr))==2) | (nrow(df)==3 & length(unique(df$strand2))==2)){
        if (df$start2[3]-df$end2[1] >0 & df$start2[3]-df$end2[1]-df$length2[2]>=(-500)){
            DF_analysis = rbind.data.frame(DF_analysis,df)
        }else{
            #Remove the middle alignment if the coordinates don't make sense
            DF_analysis = rbind.data.frame(DF_analysis, df[c(1,3),]) 
        }
    } else {
        df_1 = df
        df_2 = NULL
        #Split the 3 liners into 2 under the following conditions: all three alignments align to the same chr and same strand
        if ((nrow(df)==3) & (length(unique(df$chr))==1) & (length(unique(df$strand2))==1)){
          df = df_1[1:2,]
          df_2 = df_1[2:3,]
          #df_1 and df_2 name2 need to be changed
          df$name2 = paste0(df$name2, "_1")
          df_2$name2 = paste0(df_2$name2, "_2")
          DF_analysis = rbind.data.frame(DF_analysis, df, df_2)
        }
    }
  }
}

DF_analysis[,c(2:5,7:9)] <- as.data.frame(lapply(DF_analysis[,c(2:5,7:9)], as.numeric))
analysis_name = unique(DF_analysis$name2)
DF_SV = NULL
for (i in 1:length(analysis_name)){
    ID = analysis_name[i]
    df = NULL
    df = DF_analysis[DF_analysis$name2==ID,]

    if (nrow(df)==2){
        query_start = df$end2[1]
        query_end = df$start2[2]
        query_gap_size = query_end - query_start
        ref_chr = df$chr[1]
        if (unique(df$strand2)=="+"){
            ref_start = df$end1[1]
            ref_end = df$start1[2]
            strand = "+"
        } else{
            ref_start = df$end1[2]
            ref_end = df$start1[1]
            strand = "-"
        }
        ref_gap_size = ref_end - ref_start
        
        if (query_gap_size==0){
            query_gap_size_ratio = 1
        }else query_gap_size_ratio = query_gap_size
        if (ref_gap_size==0){
            ref_gap_size_ratio = 1
        }else ref_gap_size_ratio = ref_gap_size
        
        size = abs(ref_gap_size - query_gap_size)
        type = NULL
        
        #Calculate ratio of query_gap_size/ref_gap_size
        if (ref_gap_size_ratio > 0 & query_gap_size_ratio > 0) {
            ratio = query_gap_size_ratio/ref_gap_size_ratio
            if (ratio >=1.5) {
              type = "Insertion"
            } else {next}
        } else if (ref_gap_size_ratio < 0 & query_gap_size_ratio > 0) {
            type = "Insertion"
        } else {next}
        
    } else {next} #Probably an inversion 
    
    #Need to make sure ref_start is smaller than ref_end when reporting SVs
    #Same with ref coordinates
    start = min(ref_start, ref_end)
    end = max(ref_start, ref_end)
    ref_start = start
    ref_end = end
    
    start = min(query_start, query_end)
    end = max(query_start, query_end)
    query_start = start
    query_end = end
    
    #Also report the precise query breakpoint
    query_contig = as.numeric(str_split_fixed(df$name2[1], ':|-', 3)[,1])
    start_point = as.numeric(str_split_fixed(df$name2[1], ':|-', 3)[,2])
    first_breakpoint = start_point + query_start
    second_breakpoint = start_point + query_end
    query_breakpoint_coord = paste0(query_contig, ":", first_breakpoint, "-", second_breakpoint)
    
    #Insertion with over 800bp overlap causes error 
    if (type == "Insertion" & ref_gap_size <= -800) next
    else if (size >= 100000) next
    
    newline = c(ref_chr, ref_start, ref_end, ID, size, strand, type, ref_gap_size, query_gap_size, query_breakpoint_coord, opt$version, opt$sample)
    DF_SV = rbind(DF_SV, newline)
    
}

DF_SV = as.data.frame(DF_SV, stringsAsFactors = F)
colnames(DF_SV) = c("ref_chr", "ref_start", "ref_end", "ID", "size", "strand", "type", "ref_gap_size", "query_gap_size", "query_breakpoint_coord", "pseudohap", "sample")

# make sure there are no duplicates
# remove NUI if it anchors to chrY
DF_SV = DF_SV[(!(duplicated(DF_SV[,c('ref_start','ref_end')])) & DF_SV$ref_chr != "Y"),]

#Put blackList and segDup in GRange objects
blackList.gr = GRanges(blackList)
segDup.gr = GRanges(segDup)
SV.gr = GRanges(seqnames = Rle(as.character(DF_SV$ref_chr)), IRanges(start = as.numeric(DF_SV$ref_start), end = as.numeric(DF_SV$ref_end)))
elementMetadata(SV.gr) = DF_SV$ID  

#Remove SVs in the blackList and the segDup list
ov1 = findOverlaps(SV.gr,blackList.gr, type = "any")
ov2 = findOverlaps(SV.gr,segDup.gr, type = "any")
discard = unique(sort(c(ov1@from, ov2@from)))
keepID = SV.gr[-discard]@elementMetadata$X

#Output final SV list  
DF_SV = DF_SV[DF_SV$ID %in% keepID,]

write.table(DF_SV, paste0("SV_final_",opt$version,".txt"), col.names = F, row.names = F, quote = F, sep = '\t')

