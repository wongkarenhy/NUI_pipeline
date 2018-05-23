# This script combines results from the two pseudohaplotype 

inNameSV2.1 = "SV_final_2.1_CLEAN.txt"
inNameSV2.2 = "SV_final_2.2_CLEAN.txt"

SV2.1 = read.table(inNameSV2.1, header = F, stringsAsFactors = F, sep = '\t')
SV2.2 = read.table(inNameSV2.2, header = F, stringsAsFactors = F, sep = '\t')

SV_combined = rbind.data.frame(SV2.1, SV2.2)
colnames(SV_combined) = c("ref_chr","ref_start","ref_end","ID","size","strand","type","ref_gap_size","query_gap_size","query_breakpoint_coord", "pseudohap", "sample")
SV_combined = SV_combined[order(SV_combined$ref_chr, SV_combined$ref_start, SV_combined$ref_end),]

SV_combined$ref_start_diff = c(NA,diff(SV_combined$ref_start))
SV_combined$ref_end_diff = c(NA,diff(SV_combined$ref_end))

#If both breakpoints are within 10bp, combine the two and make genotype = 2
rows_to_be_combined = NULL
for (i in 2:nrow(SV_combined)){
    if ((SV_combined$ref_start_diff[i] <= 10 ) & (SV_combined$ref_end_diff[i] <= 10) & (SV_combined$pseudohap[i]!=SV_combined$pseudohap[(i-1)])){
        if (grepl("_",SV_combined$ID[i]) & (SV_combined$size[i] - SV_combined$size[(i-1)] < 10)){
            rows_to_be_combined = c(rows_to_be_combined, i)
        }else if (!grepl("_",SV_combined$ID[i])){
            rows_to_be_combined = c(rows_to_be_combined, i)
        }
    }
}
rows_to_be_combined_one_above = rows_to_be_combined - 1 #Get genotype 2

#-------Assign genotype-----------
SV_combined$genotype = 1
SV_combined[rows_to_be_combined_one_above, "genotype"] = 2

#Remove rows that are duplicated (appears in both psuedohaps)
SV_combined = SV_combined[-rows_to_be_combined,]

#Remove the diff columns
SV_combined = SV_combined[,c(1:12,15)]
write.table(SV_combined, "SV_final_combined.txt", row.names = F, col.names = T, quote = F, sep = '\t')



