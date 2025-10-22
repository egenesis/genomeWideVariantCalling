# Load necessary libraries
library(data.table)
library(stringr)

setwd("~/Library/CloudStorage/OneDrive-eGenesisBio/Computational Biology Team/Sanaz/GenomeWideVariantCalling/Rscripts")

table_cl <- fread("Results/Clair3_MCB_MasterTable_simplified.csv") 
table_dv <- fread("Results/DeepVariant_MCB_MasterTable_simplified.csv") 


##### QUAL metric filtered #####
MCB_qualmetric_filtered_Clair3 <- table_cl[MCB_filter == TRUE, .(variant_id)]
MCB_qualmetric_filtered_DV <- table_dv[MCB_filter == TRUE, .(variant_id)]

##### QUAL metric filtered + MCB-specific #####
MCB_only_qualmetric_filtered_Clair3 <- table_cl[MCB_specific == TRUE & MCB_filter == TRUE, .(variant_id) ]
MCB_only_qualmetric_filtered_DV <- table_dv[MCB_specific == TRUE & MCB_filter == TRUE, .(variant_id) ]
                                            

##### + Sequence feature filters #####
MCB_only_qualmetric_seqfeat_filtered_Clair3 <- table_cl[
  MCB_specific == TRUE & 
  MCB_filter == TRUE &
  in_any_feature == FALSE,
  .(variant_id)  
]
MCB_only_qualmetric_seqfeat_filtered_DV <- table_dv[
  MCB_specific == TRUE &
  MCB_filter == TRUE &
    in_any_feature == FALSE,
  .(variant_id)  
]


##### + Protein-coding and intron/exon filters #####

### Exon or UTR ###
MCB_only_qualmetric_seqfeat_filtered_in_exon_UTR_PC_Clair3 <- table_cl[
  MCB_specific == TRUE &
    MCB_filter == TRUE &
    in_any_feature == FALSE &
    in_exon_or_utr_pc == TRUE,
  .(variant_id)  
]
MCB_only_qualmetric_seqfeat_filtered_in_exon_UTR_PC_DV <- table_dv[
  MCB_specific == TRUE &
    MCB_filter == TRUE &
    in_any_feature == FALSE &
    in_exon_or_utr_pc == TRUE,
  .(variant_id)  
]

### Intron protein-coding ###
MCB_only_qualmetric_seqfeat_filtered_in_intron_PC_Clair3 <- table_cl[
  MCB_specific == TRUE &
    MCB_filter == TRUE &
    in_any_feature == FALSE &
    in_nonexonic_pc == TRUE,
  .(variant_id)  
]
MCB_only_qualmetric_seqfeat_filtered_in_intron_PC_DV <- table_dv[
  MCB_specific == TRUE &
    MCB_filter == TRUE &
    in_any_feature == FALSE &
    in_nonexonic_pc == TRUE,
  .(variant_id)  
]

#~~~~~~~~~~~~~~#
#### Union #####
#~~~~~~~~~~~~~~#

step1 <- union(table_cl$variant_id, table_dv$variant_id)
step2 <- union(MCB_qualmetric_filtered_Clair3$variant_id, MCB_qualmetric_filtered_DV$variant_id)
step3 <- union(MCB_only_qualmetric_filtered_Clair3$variant_id, MCB_only_qualmetric_filtered_DV$variant_id)
step4 <- union(MCB_only_qualmetric_seqfeat_filtered_Clair3$variant_id, MCB_only_qualmetric_seqfeat_filtered_DV$variant_id)
step5 <- Reduce(union, list(MCB_only_qualmetric_seqfeat_filtered_in_exon_UTR_PC_Clair3$variant_id, MCB_only_qualmetric_seqfeat_filtered_in_exon_UTR_PC_DV$variant_id,
               MCB_only_qualmetric_seqfeat_filtered_in_intron_PC_Clair3$variant_id, MCB_only_qualmetric_seqfeat_filtered_in_intron_PC_DV$variant_id))
step6 <- union(MCB_only_qualmetric_seqfeat_filtered_in_exon_UTR_PC_Clair3$variant_id, MCB_only_qualmetric_seqfeat_filtered_in_exon_UTR_PC_DV$variant_id)
step7 <- union(MCB_only_qualmetric_seqfeat_filtered_in_intron_PC_Clair3$variant_id, MCB_only_qualmetric_seqfeat_filtered_in_intron_PC_DV$variant_id)

message('Union results MCB - Clair3 + DeepVariant')
message("#Variants Raw, unfiltered: ", length(step1))
message("#Variants after Quality filtered: ", length(step2))
message("#Variants MCB-specific, Quality filtered: ", length(step3))
message("#Variants MCB-specific, Quality filtered, Region filtered: ", length(step4))
message("#Variants MCB-specific, Quality filtered, Region filtered, PC genes: ", length(step5))
message("#Variants MCB-specific, Quality filtered, Region filtered, PC exon/UTR: ", length(step6))
message("#Variants MCB-specific, Quality filtered, Region filtered, PC intron: ", length(step7))

#~~~~~~~~~~~~~~#
#### Prepare VEP ready inputs #####
#~~~~~~~~~~~~~~#

### Exonic ###

# Parse the variant IDs
exons_parsed <- tstrsplit(MCB_only_qualmetric_seqfeat_filtered_in_exon_UTR_PC_Union, "_", fixed = TRUE)

# Assign to columns
chrom <- exons_parsed[[1]]
ref <- exons_parsed[[3]]
alt <- exons_parsed[[4]]
start <- as.integer(exons_parsed[[2]])

end <- ifelse(nchar(ref) > 1, start + nchar(ref) - 1, start)
# Adjust start if it's an insertion (i.e., ref is 1 base and alt is longer)
# start <- ifelse(nchar(ref) == 1 & nchar(alt) > 1, end + 1, start)


# Construct allele and strand
allele <- paste0(ref, "/", alt)
strand <- "+"

# Combine into a data.frame
vep_df_exon <- data.frame(
  chromosome = chrom,
  start = start,
  end = end,
  allele = allele,
  strand = strand,
  identifier = MCB_only_qualmetric_seqfeat_filtered_in_exon_UTR_PC_Union,
  stringsAsFactors = FALSE
)

fwrite(vep_df_exon, file = "MCB-union_exon-UTR_PC_vep_input.txt", sep = " ", col.names = FALSE)


### Intronic ###

# Parse the variant IDs
introns_parsed <- tstrsplit(MCB_only_qualmetric_seqfeat_filtered_in_intron_PC_Union, "_", fixed = TRUE)

# Assign to columns
chrom <- introns_parsed[[1]]
ref <- introns_parsed[[3]]
alt <- introns_parsed[[4]]
start <- as.integer(introns_parsed[[2]])

end <- ifelse(nchar(ref) > 1, start + nchar(ref) - 1, start)
# Adjust start if it's an insertion (i.e., ref is 1 base and alt is longer)
# start <- ifelse(nchar(ref) == 1 & nchar(alt) > 1, end + 1, start)

# Construct allele and strand
allele <- paste0(ref, "/", alt)
strand <- "+"

# Combine into a data.frame
vep_df_intron <- data.frame(
  chromosome = chrom,
  start = start,
  end = end,
  allele = allele,
  strand = strand,
  identifier = MCB_only_qualmetric_seqfeat_filtered_in_intron_PC_Union,
  stringsAsFactors = FALSE
)

fwrite(vep_df_intron, file = "MCB-union_intron_PC_vep_input.txt", sep = " ", col.names = FALSE)



