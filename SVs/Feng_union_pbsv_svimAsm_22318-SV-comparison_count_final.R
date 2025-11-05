#### Load libraries ####
# updated by Feng on Oct 16 2025.
# union of pbsv and svim-asm

#install.packages("vcfR")
#BiocManager::install("rtracklayer")
#BiocManager::install("txdbmaker")
#BiocManager::install("GenomicFeatures")
#install.packages("tidyr")

library(data.table)
#library(tidyverse)
library(stringr)
library(vcfR)
library(purrr)
library(dplyr)
library(GenomicRanges)
library(rtracklayer)
library(GenomicFeatures)
library(rlang)
library(ggplot2)
library(ggsci)

library(tidyr)


setwd("/home/rstudio/genome-wide-SV/")
#-----------------------------------#
##### Load & preprocess VCF(s) ######
#-----------------------------------#

##### Read in data #####
VCF ="/data5/sv/sam_SV/Duroc/union_pbsv_svim-asm/Yuc104F_22318_to_Ssc11.1_pbsv_svimAsm_Survivor1100_union.vcf"
CSV = "/data5/sv/sam_SV/Duroc/union_pbsv_svim-asm/Yuc104F_22318_to_Ssc11.1_pbsv_svimAsm_Survivor1100_union_VCF_data.tsv"

SV_vcf <- read.vcfR(VCF)
fixed <- as.data.frame(getFIX(SV_vcf)) # Fixed fields: CHROM, POS, ID, REF, ALT, QUAL, FILTER
info <- extract_info_tidy(SV_vcf) # INFO fields: parsed into a matrix
df <- cbind(fixed, info) %>% dplyr::select(-SUPP_VEC)# Merge

#### pre-process ####
df$POS <- as.integer(df$POS)
df$END <- as.integer(df$END)
df$SVLEN <- as.integer(df$SVLEN)
df$SUPP <- as.integer(df$SUPP)

# Data rows only
SV_vcf_data <- fread(CSV) %>%
  dplyr::select(CHROM, POS, ID, REF, ALT, FORMAT, Yuc104F, `22318`)


### Widen data to include sample-specific columns
format_fields <- strsplit(SV_vcf_data$FORMAT[1], ":")[[1]] # Split the FORMAT field into a vector

# Function to extract sample-specific data and for FORMAT fields
expand_sample_column <- function(sample_col, format_fields, prefix) {
  sample_col_clean <- as.character(sample_col)
  placeholder <- paste(rep(".", length(format_fields)), collapse = ":")   # Ensure no NA or unexpected blanks
  sample_col_clean[sample_col_clean == "" | is.na(sample_col_clean)] <- placeholder
  split_list <- strsplit(sample_col_clean, ":")   # Split and bind
  split_matrix <- do.call(rbind, split_list)
  if (ncol(split_matrix) != length(format_fields)) {   # Confirm success before converting
    stop("Mismatch between split fields and FORMAT field count")}
  as.data.frame(split_matrix, stringsAsFactors = FALSE) %>%   # Convert to dataframe
    setNames(paste0(prefix, "_", format_fields))
}

# Expand both sample columns & recombine
yuc_expanded <- expand_sample_column(SV_vcf_data$Yuc104F, format_fields, "Yuc104F")
S22318_expanded <- expand_sample_column(SV_vcf_data$`22318`, format_fields, "22318")
SV_vcf_wide <- cbind(SV_vcf_data, yuc_expanded, S22318_expanded)  %>%   dplyr::select(-Yuc104F, -`22318`,-FORMAT)

# combine with vcfR-parsed DF & remove unneeded columns
SVs <- inner_join(df,SV_vcf_wide, by = c("CHROM","POS","ID","REF","ALT")) %>% dplyr::select(-RE, -MAPQ, -PRECISE, -IMPRECISE)

# If 4 callers were passed to SURVIVOR
#SVs$Yuc104F_PSV[SVs$Yuc104F_PSV=="NaN"] <- 0000
#SVs$22318_PSV[SVs$22318_PSV=="NaN"] <- 0000

# If 2 callers were passed to SURVIVOR
SVs$Yuc104F_PSV[SVs$Yuc104F_PSV=="NaN"] <- 00
SVs$`22318_PSV`[SVs$`22318_PSV`=="NaN"] <- 00

# Function to parse SV support and widen into binary values per caller and sample
parse_supp_vec <- function(vec, caller_names, prefix) {
  do.call(rbind, lapply(vec, function(x) {
    bits <- strsplit(x, "")[[1]]
    bits <- c(bits, rep("0", length(caller_names) - length(bits)))
    setNames(as.integer(bits), paste0(prefix, "_", caller_names))
  }))
}

# Define caller names (Be sure order is correct)
caller_names <- c("pbsv","svim-asm")
support_Yuc <- parse_supp_vec(SVs$Yuc104F_PSV, caller_names, "Yuc104F")
support_22318 <- parse_supp_vec(SVs$`22318_PSV`, caller_names, "22318")

SVs <- cbind(SVs, support_Yuc, support_22318)
# remove chrY and MT variants
SVs <- SVs %>% filter(CHROM != "Y") %>%  filter(CHROM != "MT")

# Set variable for whether the SV was only supported in 22318, as in was wild-type or unknown for Yuc104F and non-WT for 22318
WT_GT_values <- c("./.", "0/0", "0/.", "./0") # WT values

# Make as variable to keep info in main data frame
SVs$`22318_specific` <- FALSE 
SVs$`22318_specific`[SVs$Yuc104F_GT %in% WT_GT_values & !(SVs$`22318_GT` %in%  WT_GT_values)] <- TRUE 

# Represent SV coordinates as Genomic Ranges objects
gr <- GRanges(
  seqnames = SVs$CHROM,
  ranges = IRanges(start = pmin(SVs$POS,SVs$END), end = pmax(SVs$POS,SVs$END)),
  SVTYPE = SVs$SVTYPE,
  SVLEN = SVs$SVLEN
)

##### Load BED files for sequence region filters #####
perv_bed <- fread("/data5/sv/sam_SV/BED_files/Ssc11.1_HC_41PERV.bed", col.names = c("CHROM", "START", "END","pType", "type"))    # high confidence PERV loci on Sus scrofa v11.1 (Duroc)
centro_bed <- fread("/data5/sv/sam_SV/BED_files/Ssc11.centromere_and_like_regions.bed", col.names = c("CHROM", "START", "END","feature"))    # centromere regions in Sus scrofa v11.1 (Duroc)
trf_bed <- fread("/data5/sv/sam_SV/BED_files/susScr11.trf_mod_swaped.bed", col.names = c("CHROM", "START", "END","feature"))     # modified from UCSC trf track "ftp://hgdownload.soe.ucsc.edu/goldenPath/susScr11/bigZips/susScr11.trf.bed.gz"
repeats_bed <- fread("/data5/sv/sam_SV/BED_files/nestedRepeats_swaped.bed", col.names = c("CHROM", "START", "END","name"))    # modified from UCSC repeat track containing LINE, SINE etc. https://hgdownload.soe.ucsc.edu/goldenPath/susScr11/database/. This is used in manual curation in IGV to understand the complex regions
#lowMQ10_bed_Yuc104F <- fread(file = "../projects/egenesis/sammodlin/egenesis/sammodlin/sammodlin-sandbox/manual_curation/guides/Yuc104F_HiFito_Ssc11.1.lowMQ10.bed", col.names = c("CHROM", "START", "END", "threshold"))
#lowMQ10_bed_22318 <- fread(file = "../projects/egenesis/sammodlin/egenesis/sammodlin/sammodlin-sandbox/manual_curation/guides/22318_HiFito_Ssc11.1.lowMQ10.bed", col.names = c("CHROM", "START", "END", "threshold"))
#lowMQ10_bed <- rbind(lowMQ10_bed_22318, lowMQ10_bed_Yuc104F) 

# Adjust to 1-based coordinates (for START only) to match VCF representation 
adjust_bed <- function(bed_df) {
  bed_df %>% mutate(START = START + 1) 
}

perv_bed <- adjust_bed(perv_bed)
centro_bed <- adjust_bed(centro_bed)
trf_bed <- adjust_bed(trf_bed)
repeats_bed <- adjust_bed(repeats_bed)
#lowMQ10_bed <- adjust_bed(lowMQ10_bed)

##### Annotate overlaps #####
# Define overlap detection function
annotate_overlap <- function(gr, bed_df) {
  bed_gr <- GRanges(seqnames = bed_df$CHROM,
                    ranges = IRanges(start = bed_df$START, end = bed_df$END))
  overlapsAny(gr, bed_gr)
}

# Annotate with logical column for each feature type
SVs$in_perv <- annotate_overlap(gr, perv_bed)
SVs$in_centromere  <- annotate_overlap(gr, centro_bed)
SVs$in_trf  <- annotate_overlap(gr, trf_bed)
SVs$in_repeat      <- annotate_overlap(gr, repeats_bed)
#SVs$in_lowMQ10        <- annotate_overlap(gr, lowMQ10_bed)

feature_cols <- c("in_centromere", "in_repeat", "in_trf", "in_perv")
filt_feature_cols <- c("in_centromere", "in_perv")

# Create "total" column that == TRUE if any feature is TRUE
SVs <- SVs %>%
  mutate(in_any_feature = if_any(all_of(feature_cols), ~ .x))

SVs <- SVs %>%
  mutate(in_filt_feature = if_any(all_of(filt_feature_cols), ~ .x))

#---------------------------------#
##### Genomic feature overlap #####
#---------------------------------#

#### Read in GTF file ####

gtf <- import("/data5/sv/sam_SV/BED_files/Sus_scrofa.Sscrofa11.1.106.gtf.gz")

# extract all Exons and UTRs from GTF
exons <- gtf[gtf$type == "exon"]
utrs  <- gtf[gtf$type %in% c("five_prime_utr", "three_prime_utr")]

# Subset to protein-coding
exons_pc <- exons[mcols(exons)$gene_biotype == "protein_coding"]
utrs_pc  <- utrs[mcols(utrs)$gene_biotype == "protein_coding"]

# Define genomic feature columns
gen_feature_cols <- c("in_exon", "in_utr")
gen_pc_feature_cols <- c("in_pc_exon", "in_pc_utr")

# Find overlaps
SVs$in_exon     <- overlapsAny(gr, exons)
SVs$in_utr      <- overlapsAny(gr, utrs)
SVs$in_pc_exon  <- overlapsAny(gr, exons_pc)
SVs$in_pc_utr   <- overlapsAny(gr, utrs_pc)

# Create "total" column that == TRUE if any feature is TRUE
SVs <- SVs %>%
  mutate(in_exon_or_utr = if_any(all_of(gen_feature_cols), ~ .x))

SVs <- SVs %>%
  mutate(in_exon_or_utr_pc = if_any(all_of(gen_pc_feature_cols), ~ .x))

#--------------------------------------------#
##### Genomic feature overlap: Introns #####
#--------------------------------------------#
# Create TxDb from GTF
#txdb <- makeTxDbFromGFF("/data5/sv/sam_SV/BED_files/Sus_scrofa.Sscrofa11.1.106.gtf.gz", format = "gtf")
# Build a TxDb from the GTF using GenomicFeatures
txdb <- makeTxDbFromGRanges(gtf)

# Extract introns by transcript
introns_by_tx <- intronsByTranscript(txdb, use.names = TRUE)

# Get transcript annotations with biotype info
tx2gene <- unique(gtf[gtf$type == "transcript", c("transcript_id", "gene_biotype")])

# Keep only introns from protein-coding genes
tx_ids_pc <- tx2gene$transcript_id[tx2gene$gene_biotype == "protein_coding"]
introns_pc <- introns_by_tx[names(introns_by_tx) %in% tx_ids_pc]
introns_pc_all <- reduce(unlist(introns_pc, use.names = FALSE))

# Add overlap columns
SVs$in_intron_pc <- overlapsAny(gr, introns_pc_all)

# Create a combined "intronic or exonic/UTR protein-coding" feature column
SVs <- SVs %>%
  mutate(in_exonic_or_intronic_pc = in_pc_exon | in_pc_utr | in_intron_pc)

# Write out results to file
fwrite(SVs, "./22318/Feng_Processed_union_pbsv_svimAsm_Yuc104F_var_22318_SVs.csv")

##############################################################
# Section 2: Filtering steps
##############################################################
#~~~~~~~~~~~~~~~~~~~~~#
#### Read in files ####
#~~~~~~~~~~~~~~~~~~~~~#

SVs <- fread("./22318/Feng_Processed_union_pbsv_svimAsm_Yuc104F_var_22318_SVs.csv") 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#### Generate counts & summaries ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
##### Summaries #####
### SV Type ###

svtypes <- c("DEL","DUP","INS","INV","TRA") # define character vector of SV types
# Make a list of successively stricter filter criteria
filters <- list(	
  "No filter" = expr(TRUE),							# unfiltered
  "22318-specific" = expr(`22318_specific`),					# 22318-specific
  "22318-specific + region" = expr(`22318_specific` & !in_filt_feature), 	# + sequence pattern
  "22318-specific + region + size (20–1000 bp)" =				# + SV Length (20-1000 bp)
    expr(`22318_specific` & !in_filt_feature & ((abs(SVLEN) >= 20 & abs(SVLEN) <= 1000) | SVTYPE == "TRA")),
  "22318-specific + region + size (20–500 bp)" = 				# + SV Length (20-500 bp)
    expr(`22318_specific` & !in_filt_feature & ((abs(SVLEN) >= 20 & abs(SVLEN) <= 500)  | SVTYPE == "TRA")),
  "22318-specific + region + PC (any size)" =   # + exon/intron
    expr(`22318_specific` & !in_filt_feature & in_exonic_or_intronic_pc),
  "22318-specific + region + size (20–1000 bp) + PC" = 			# + exon/intron (20-1000 bp)
    expr(`22318_specific` & !in_filt_feature & ((abs(SVLEN) >= 20 & abs(SVLEN) <= 1000) | SVTYPE == "TRA") & in_exonic_or_intronic_pc),
  "22318-specific + region + size (20–500 bp) + PC" =			# + exon/intron (20-500 bp)
    expr(`22318_specific` & !in_filt_feature & ((abs(SVLEN) >= 20 & abs(SVLEN) <= 500)  | SVTYPE == "TRA") & in_exonic_or_intronic_pc)
)

sv_summary <-
  imap_dfr(filters, ~ {
    SVs %>%
      filter(!! .x) %>%                        # apply the filter expression
      count(SVTYPE, name = "n") %>%
      complete(SVTYPE = svtypes, fill = list(n = 0)) %>%
      pivot_wider(names_from = SVTYPE, values_from = n, values_fill = 0) %>%
      mutate(Filter = .y) %>%
      dplyr::select(Filter, all_of(svtypes)) %>%
      mutate(Total = rowSums(across(all_of(svtypes))))
  }) %>%
  
  # Order rows to match the list order above
  mutate(Filter = factor(Filter, levels = names(filters))) %>%
  arrange(Filter)

fwrite(sv_summary, "./22318/Feng-22318-unionof2-SV-filtering-summary.csv")


#-------------------------------------#
##### Prepare for manual curation ##### 
#-------------------------------------#

#### Output VCFs for manual curation ####
# Filtered from union of 2 callers by everything but size (either of 2 callers, so 1 out 2)

VCF_full_filt_no_size_1_2 <- SVs %>% 
  filter(`22318_specific` == TRUE) %>%
  filter(in_filt_feature == FALSE) %>%
  filter(in_exonic_or_intronic_pc == TRUE) 

# For VCF subsetting
# 270 rows for 22318
VCF_full_filt_no_size_1_2_small <- VCF_full_filt_no_size_1_2 %>% as.data.table() %>% dplyr::select(CHROM,POS,ID)
fwrite(VCF_full_filt_no_size_1_2_small, "./22318/Feng_22318_unionOf2_full-filtered_no_size.csv")

# For manual curation excel
VCF_full_filt_no_size_1_2_small_csv <- VCF_full_filt_no_size_1_2 %>% as.data.table() %>% dplyr::select(CHROM, POS, ID, REF, ALT, SVLEN, SVTYPE, QUAL, CHR2, STRANDS, `22318_GT`, Yuc104F_GT, `22318_pbsv`, Yuc104F_pbsv, `22318_svim-asm`, `Yuc104F_svim-asm`, `22318_LN`, Yuc104F_LN, `22318_QV`, Yuc104F_QV, `22318_CO`, Yuc104F_CO, `22318_ID`, Yuc104F_ID,  `22318_TY`, Yuc104F_TY, `22318_AAL`, Yuc104F_AAL )
fwrite(VCF_full_filt_no_size_1_2_small_csv, "./22318/Feng_22318_unionOf2_full-filtered_no_size_man_cur.csv")

######################################################################
### Distribution of feature overlaps with 22318-specific SVs ###
######################################################################
summary_SVs_feat_overlap_22318_spec_all <- SVs %>% filter(`22318_specific` == TRUE) %>%
  group_by(SVTYPE) %>%
  summarise(across(all_of(c(feature_cols, "in_any_feature")), ~ sum(.x, na.rm = TRUE), .names = "{.col}")) %>%
  ungroup()
summary_SVs_feat_overlap_22318_spec_all

summary_SVs_feat_overlap_22318_spec <- SVs %>% filter(`22318_specific` == TRUE) %>%
  group_by(SVTYPE) %>%
  summarise(across(all_of(c(filt_feature_cols, "in_filt_feature")), ~ sum(.x, na.rm = TRUE), .names = "{.col}")) %>%
  ungroup()
summary_SVs_feat_overlap_22318_spec

summary_SVs_feat_overlap_22318_spec_sizefilt_all <- SVs %>% filter(`22318_specific` == TRUE) %>%
  filter((abs(SVLEN)<=1000 & abs(SVLEN)>=20) | SVTYPE == "TRA") %>%
  group_by(SVTYPE) %>%
  summarise(across(all_of(c(feature_cols, "in_any_feature")), ~ sum(.x, na.rm = TRUE), .names = "{.col}")) %>%
  ungroup()
summary_SVs_feat_overlap_22318_spec_sizefilt_all

summary_SVs_feat_overlap_22318_spec_sizefilt <- SVs %>% filter(`22318_specific` == TRUE) %>%
  filter((abs(SVLEN)<=1000 & abs(SVLEN)>=20) | SVTYPE == "TRA") %>%
  group_by(SVTYPE) %>%
  summarise(across(all_of(c(filt_feature_cols, "in_filt_feature")), ~ sum(.x, na.rm = TRUE), .names = "{.col}")) %>%
  ungroup()
summary_SVs_feat_overlap_22318_spec_sizefilt

#---------------#
##### Plots #####
#---------------#


# Keepo colors consistent:
# Get a palette with enough distinct colors
sv_palette <- pal_npg("nrc")(9)   # up to 9 colors from NPG palette

# Map them consistently to SV types
sv_palette <- pal_npg("nrc")(9)   # up to 9 colors from NPG palette
sv_colors <- c(
  DEL = sv_palette[1],
  INS = sv_palette[2],
  DUP = sv_palette[3],
  INV = sv_palette[4],
  TRA = sv_palette[5]
)


# Distribution facetted by SURVIVOR
p1 <- SVs %>% ggplot(aes(x = log10(abs(SVLEN)+1), fill = SVTYPE)) +
  geom_histogram(bins = 50, alpha = 0.9, position = "stack") +
  scale_fill_manual(values = sv_colors, drop = FALSE) +
  theme_minimal(base_size = 14) +
  labs(x = "SV Length (log10)", y = "Count") +
  theme(legend.position = "top", strip.background = element_blank(), strip.text = element_text(face = "bold")
  )


p2 <- SVs %>% filter(`22318_specific` == TRUE) %>% ggplot(aes(x = log10(abs(SVLEN)+1), fill = SVTYPE)) +
  geom_histogram(bins = 50, alpha = 0.9, position = "stack") +
  scale_fill_manual(values = sv_colors, drop = FALSE) +  theme_minimal(base_size = 14) +
  labs(x = "SV Length (log10)", y = "Count") +
  theme(legend.position = "top", strip.background = element_blank(), strip.text = element_text(face = "bold")
  )

p3 <- SVs %>% filter(`22318_specific` == TRUE) %>%
  filter(in_filt_feature == FALSE) %>%
  ggplot(aes(x = log10(abs(SVLEN)+1), fill = SVTYPE)) +
  geom_histogram(bins = 50, alpha = 0.9, position = "stack") +
  scale_fill_manual(values = sv_colors, drop = FALSE) +
  theme_minimal(base_size = 14) +
  labs(x = "SV Length (log10)", y = "Count") +
  theme(legend.position = "none", strip.background = element_blank(), strip.text = element_text(face = "bold")
  )

p4 <- SVs %>% filter(`22318_specific` == TRUE) %>%
  filter(in_filt_feature == FALSE) %>%
  filter((abs(SVLEN)<=1000 & abs(SVLEN)>=20) | SVTYPE == "TRA") %>%
  ggplot(aes(x = abs(SVLEN), fill = SVTYPE)) +
  geom_histogram(bins = 50, alpha = 0.9, position = "stack") +
  scale_fill_manual(values = sv_colors, drop = FALSE) +
  theme_minimal(base_size = 14) +
  labs(x = "SV Length", y = "Count") +
  theme(legend.position = "none", strip.background = element_blank(), strip.text = element_text(face = "bold")
  )

p5 <- SVs %>% filter(`22318_specific` == TRUE) %>%
  filter(in_filt_feature == FALSE) %>%
  filter((abs(SVLEN)<=500 & abs(SVLEN)>=20) | SVTYPE == "TRA") %>%
  ggplot(aes(x = abs(SVLEN), fill = SVTYPE)) +
  geom_histogram(bins = 50, alpha = 0.9, position = "stack") +
  scale_fill_manual(values = sv_colors, drop = FALSE) +
  theme_minimal(base_size = 14) +
  labs(x = "SV Length", y = "Count") +
  theme(legend.position = "none", strip.background = element_blank(), strip.text = element_text(face = "bold")
  )

p6 <- SVs %>% filter(`22318_specific` == TRUE) %>%
  filter(in_filt_feature == FALSE) %>%
  filter(in_exonic_or_intronic_pc == TRUE) %>%
  ggplot(aes(x = log10(abs(SVLEN)+1), fill = SVTYPE)) +
  geom_histogram(bins = 50, alpha = 0.9, position = "stack") +
  scale_fill_manual(values = sv_colors, drop = FALSE) +
  theme_minimal(base_size = 14) +
  labs(x = "SV Length (log10)", y = "Count") +
  theme(legend.position = "none", strip.background = element_blank(), strip.text = element_text(face = "bold")
  )

p7 <- SVs %>% filter(`22318_specific` == TRUE) %>%
  filter(in_filt_feature == FALSE) %>%
  filter((abs(SVLEN)<=500 & abs(SVLEN)>=20) | SVTYPE == "TRA") %>%
  filter(in_exonic_or_intronic_pc == TRUE) %>%
  ggplot(aes(x = abs(SVLEN), fill = SVTYPE)) +
  geom_histogram(bins = 50, alpha = 0.9, position = "stack") +
  scale_fill_manual(values = sv_colors, drop = FALSE) +
  theme_minimal(base_size = 14) +
  labs(x = "SV Length", y = "Count") +
  theme(legend.position = "none", strip.background = element_blank(), strip.text = element_text(face = "bold")
  )

p8 <- SVs %>% filter(`22318_specific` == TRUE) %>%
  filter(in_filt_feature == FALSE) %>%
  filter((abs(SVLEN)<=1000 & abs(SVLEN)>=20) | SVTYPE == "TRA") %>%
  filter(in_exonic_or_intronic_pc == TRUE) %>%
  ggplot(aes(x = abs(SVLEN), fill = SVTYPE)) +
  geom_histogram(bins = 50, alpha = 0.9, position = "stack") +
  scale_fill_manual(values = sv_colors, drop = FALSE) +
  theme_minimal(base_size = 14) +
  labs(x = "SV Length", y = "Count") +
  theme(legend.position = "none", strip.background = element_blank(), strip.text = element_text(face = "bold")
  )

# Print plots
p1
p2
p3
p4
p5
p6
p7



# 
# GT_summary <- SVs %>% 
#   count(Yuc104F_GT, name = "count") #%>%
# #  pivot_wider(names_from = SVTYPE, values_from = count, values_fill = 0)
# 
# S22318_GT_summary <- SVs %>% 
#   count(`22318_GT`, name = "count") #%>%
# #  pivot_wider(names_from = SVTYPE, values_from = count, values_fill = 0)
# 
# SUPP_summary <- SVs %>% 
#   count(SUPP, name = "count") #%>%
# #  pivot_wider(names_from = SVTYPE, values_from = count, values_fill = 0)
