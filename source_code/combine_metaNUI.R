library(GenomicRanges)
library(stringr)

#sample_list = c("NA19068", "NA19440", "HG00250", "HG00851", "HG03115", "HG00353", "HG03838", "NA20587", "HG01971", "HG02623", "NA19789", "NA21125", "NA18552", "NA19921", "HG00512", "HG00733", "NA19240")
sample_list = c("HG00250", "HG00353")

SV_meta = NULL
total_SV_per_individual = NULL
for (i in 1:length(sample_list)){
    sample = sample_list[i]
    SV = read.table(paste0(sample,"/NUI/SV_final_combined_NR_FILTERED.txt"), header = F, stringsAsFactors = F, sep = '\t')
    total_SV_per_individual = rbind.data.frame(total_SV_per_individual, cbind(nrow(SV), sum(SV[,9])), stringsAsFactors = F)
    SV_meta = rbind.data.frame(SV_meta, SV)
}

colnames(SV_meta) = c("ref_chr", "ref_start", "ref_end", "ID", "size", "strand", "type", "ref_gap_size", "query_gap_size", "query_breakpoint_coord", "pseudohap", "sample", "geno")
SV_meta = SV_meta[order(SV_meta$ref_chr, SV_meta$ref_start, SV_meta$ref_end),]
SV_meta$sample = paste0(SV_meta$sample,":", SV_meta$pseudohap, ":",SV_meta$geno)

#-------Get a list of breakpoints----------------
uniq_breakpoints = SV_meta[!duplicated(SV_meta[,c("ref_start")]),]
uniq_breakpoints = uniq_breakpoints[!duplicated(uniq_breakpoints[,c("ref_end")]),]
uniq_breakpoints$ref_start_diff = c(NA,diff(uniq_breakpoints$ref_start))
uniq_breakpoints$ref_end_diff = c(NA,diff(uniq_breakpoints$ref_end))
discard = which(abs(uniq_breakpoints$ref_start_diff)< 50 & abs(uniq_breakpoints$ref_end_diff)<50)
uniq_breakpoints = uniq_breakpoints[-discard,1:13]

#Get rid of duplicated sampels
#--------Determine who has what breakpoints--------------
count = NULL
for (i in 1:nrow(uniq_breakpoints)){
    occuring_sample = SV_meta[SV_meta$ref_chr==uniq_breakpoints[i,"ref_chr"] & (SV_meta$ref_start==uniq_breakpoints[i,"ref_start"] | SV_meta$ref_end==uniq_breakpoints[i,"ref_end"] | (abs(SV_meta$ref_start-uniq_breakpoints[i,"ref_start"]) < 50 & abs(SV_meta$ref_end - uniq_breakpoints[i,"ref_end"]) < 50 )),"sample"]
    counts = str_split_fixed(occuring_sample, ":", 4)[,1]
    counts = unique(counts)
    count[i] = length(unique(counts))
    uniq_breakpoints[i,"sample_combined"] = do.call(paste, c(as.list(occuring_sample), sep=";"))
}
uniq_breakpoints$counts = count

#--------Rearrange the dataframe----------
data = uniq_breakpoints[,c(1:10,14,15)]

write.table(data_final, "NUI_combined_all.txt", row.names = F, col.names = T, quote = F, sep = '\t')

