#### Load libraries ####
library(data.table)
library(tidyverse)
# install.packages("vcfR")
library("vcfR")
library(purrr)
library(dplyr)
library(tidyr)


setwd("~/Library/CloudStorage/OneDrive-eGenesisBio/Computational Biology Team/Sanaz/GenomeWideVariantCalling/Rscripts")


clair3_vcf <- read.vcfR("Inputfiles/deepvariant.cohort.biallelic.norm.vcf.gz")
# clair3_vcf <- clair3_vcf[1:20000, ]

# ensure all IDs are unique
clair3_vcf@fix[, "ID"] <- paste0("var", seq_len(nrow(clair3_vcf@fix)))

# Convert to tidy data format
tidy_clair3_vcf <- vcfR2tidy(
  clair3_vcf,
  single_frame = TRUE,
  format_fields = c("GT", "DP", "AD"),
  info_fields = c("AQ"),
  dot_is_NA = FALSE # preserve "./." in gt_GT -- default behavior is to convert them to NA
)
remove(clair3_vcf)

# de-convert NAs back to "./." -- vcf2tidy() does this by default and as far as I can tell there is not an option to suppress this behavior
tidy_clair3_vcf$dat <- tidy_clair3_vcf$dat %>% mutate(gt_GT = ifelse(is.na(gt_GT), "./.", gt_GT))

# Add unique variant ID variable
tidy_clair3_vcf$dat <- tidy_clair3_vcf$dat %>% mutate(variant_id = paste(CHROM, POS, REF, ALT, sep = "_"))

# Reshape to wide format
tidy_clair3_vcf$dat$Indiv[tidy_clair3_vcf$dat$Indiv == "GC64-23_008_MCB_FSDC"] <- "MCB" # Rename for simplicity
tidy_clair3_vcf$dat$Indiv <- factor(tidy_clair3_vcf$dat$Indiv) # pre-factor for speed

wide <- tidy_clair3_vcf$dat %>%
  dplyr::select(variant_id, CHROM, POS, REF, ALT, Indiv, gt_GT, gt_DP, gt_AD, QUAL) %>%
  tidyr::pivot_wider(names_from = Indiv, values_from = c(gt_GT, gt_DP, gt_AD, QUAL))

remove(tidy_clair3_vcf)

# Separate allele depths into REF and ALT
wide <- wide %>% 
  separate(gt_AD_MCB, into = c("AD_REF_MCB", "AD_ALT_MCB"), sep = ",", convert = TRUE) %>%
  separate(gt_AD_Yuc104F, into = c("AD_REF_Yuc104F", "AD_ALT_Yuc104F"), sep = ",", convert = TRUE)

## add var_caller & ref_genome
# wide$ref_genome <- "Duroc"
# wide$var_caller <- "Clair3"
# wide$var_caller <- "DeepVariant"
wide$caller_ref <- "DeepVariant_Duroc"

#---------------------------------------------------------------------------#
##### Identify & catalog problematic rows & enforce correct data types #####
#---------------------------------------------------------------------------#

# Record rows where either GT col is multivalued
wide_bad <- wide %>% 
  mutate(row = row_number()) %>%
  filter(map_int(gt_GT_MCB, length) > 1 |
           map_int(gt_GT_Yuc104F, length) > 1) %>%
  mutate(culprit = case_when(
    map_int(gt_GT_MCB, length) > 1 ~ "gt_GT_MCB",
    map_int(gt_GT_Yuc104F, length) > 1 ~ "gt_GT_Yuc104F"
  ))

message("Number of rows with multiple genotype entries: ", nrow(wide_bad))
# Proceed with assuming nrow(wide_bad) > 1

### Helper Functions #####

# Helper functions to pull first element or NA, then coerce if needed
first_or_na <- function(x) {
  if (length(x) >= 1) x[[1]] else NA_character_
}

first_num <- function(x) {
  v <- first_or_na(x)
  # if it’s already NA, ".", or empty, return an NA_real_ (NA but as a numeric rather than character)
  if (is.na(v) || v == ".") return(NA_real_)
  # otherwise coerce safely
  as.numeric(v)
}

## convert "." to NA with correct data formats
wide <- wide %>%         
  mutate(
    # GT to character
    across(starts_with("gt_GT_"), ~ map_chr(.x, first_or_na)),
    
    # DP & GQ to numeric
    across(matches("^gt_(DP|GQ)_"), ~ map_dbl(.x, first_num)),
    
    # AD_REF & AD_ALT to numeric
    across(matches("^AD_(REF|ALT)_"), ~ map_dbl(.x, first_num)),
    
    # QUAL to numeric
    across(starts_with("QUAL_"), ~ map_dbl(.x, first_num))
  ) 

  # %>% mutate(
  #   alt_ratio_MCB = AD_ALT_MCB /
  #     (AD_ALT_MCB + AD_REF_MCB),
  #   alt_ratio_Yuc104F = AD_ALT_Yuc104F /
  #     (AD_ALT_Yuc104F + AD_REF_Yuc104F)
  # )

#-----------------------#
##### Clean dataset ##### 
#-----------------------#


# Clean up chromosome names, add variant type, & make single variable for VC & ref genome combo
get_variant_type <- function(ref, alt) {
  if (nchar(ref) == 1 && nchar(alt) == 1) {
    return("SNP")
  } else if (nchar(ref) < nchar(alt)) {
    return("Insertion")
  } else if (nchar(ref) > nchar(alt)) {
    return("Deletion")
  } else {
    return("Complex")
  }
}

wide <- wide %>% 
  mutate(
    # CHROM_clean = gsub("^chr(\\d+)_.*", "\\1", CHROM),     [no need for adding another chrom column]
    variant_type = mapply(get_variant_type, REF, ALT),
  )

# remove variants in chrY and MT
wide <- wide %>%
  filter(!CHROM %in% c("Y", "MT"))

# Add values of 0 to NA for columns that should be zero when NA
setDT(wide)
na_zero_cols <- grep("^(AD_REF_|AD_ALT_)", names(wide), value = TRUE) # grab relevant columns 
wide[ , (na_zero_cols) := lapply(.SD, function(x) fifelse(is.na(x), 0, x)),
      .SDcols = na_zero_cols ]

#----------------------------------------------------------------------------------#
##### Classify variants by their similarity/difference between Yuc104F & MCB #####
#----------------------------------------------------------------------------------#

# add MCB_info, Yuc104F_info, call_category columns
classify_variant_pair <- function(df, gt_cols) {
  # gt_cols: character vector of length 2, e.g., c("gt_GT_MCB","gt_GT_Yuc104F")
  DT <- as.data.table(df)
  
  # extract suffixes & build the info column names
  suffs     <- sub("^gt_GT_", "", gt_cols)
  info_cols <- paste0(suffs, "_info")
  
  # compute info-status ("none"/"partial"/"full")
  DT[, (info_cols) := lapply(.SD, function(gt) {
    gt <- as.character(gt)   # ensure NAs are real R NAs, not factor levels
    fifelse(
      is.na(gt) | gt == "./.",      "none",     # NA _or_ literal “./.” 
      fifelse(grepl("\\.", gt),     "partial",  # any single “.”
              "full")     # no “.” and not NA
    )
  }),
  .SDcols = gt_cols
  ]
  
  # classification into a new column
  # alias the two genotype columns
  col1 <- gt_cols[1]
  col2 <- gt_cols[2]
  
  DT[, call_category := fcase(
    # exact match if both are fully called → shared
    (
      get(col1) == get(col2) &
        get(info_cols[1]) == "full" &
        get(info_cols[2]) == "full"
    ) |
      # pure het vs pure het → shared
      (
        get(col1) %in% c("0/1","1/0") &
          get(col2) %in% c("0/1","1/0")
      ),
    "shared",
    
    # missing/partial vs alt → different_call
    (
      (is.na(get(col1)) | get(col1) %in% c("./.","./0","0/."))
      & grepl("1", get(col2))
    ) |
      (
        (is.na(get(col2)) | get(col2) %in% c("./.","./0","0/."))
        & grepl("1", get(col1))
      ),
    "different_call",
    
    # ref vs any containing "1" → different_call
    (get(col1) == "0/0" & grepl("1", get(col2))) |
      (get(col2) == "0/0" & grepl("1", get(col1))),
    "different_call",
    
    # het vs hom-alt → different_call
    (get(col1) %in% c("0/1","1/0") & get(col2) == "1/1") |
      (get(col2) %in% c("0/1","1/0") & get(col1) == "1/1"),
    "different_call",
    
    # fallback → ambiguous
    default = "ambiguous"
  )]
  
  # return table with new columns
  return(DT[])
}

wide <- as.data.table(wide)
wide <- classify_variant_pair(wide, c("gt_GT_MCB", "gt_GT_Yuc104F"))

#--------------------------------------#
##### Annotate additional features ##### 
#--------------------------------------#


### Define sample-specific variants ###
wide <- wide %>% 
  mutate(sample_specificity = case_when(
    call_category == "shared" ~ "shared",
    call_category == "ambiguous" ~ "ambiguous",
    call_category == "different_call" & str_detect(gt_GT_Yuc104F, "1") & !str_detect(gt_GT_MCB, "1") ~ "Yuc104F_specific",
    call_category == "different_call" & str_detect(gt_GT_MCB, "1") & !str_detect(gt_GT_Yuc104F, "1") ~ "MCB_specific",
    call_category == "different_call" & str_detect(gt_GT_Yuc104F, "1") & str_detect(gt_GT_MCB, "1") ~ "complex"
  ))


########## Annotate Zygosity ###########

# Annotate per-sample zygosities
zygosity_func <- function(gt) {
  if (is.na(gt)) return("other")
  if (gt %in% c("0/1", "1/0")) return("het")
  if (gt == "0/0") return("hom_ref")
  if (gt == "1/1") return("hom_alt")
  return("other")
}

wide <- wide %>% 
  mutate(
    zygosity_MCB = sapply(gt_GT_MCB, zygosity_func),
    zygosity_Yuc104F = sapply(gt_GT_Yuc104F, zygosity_func)
  )

# Annotate overall zygosity
wide <- wide %>% 
  mutate(overall_zygosity = case_when(
    call_category == "shared" ~ zygosity_MCB,  # could be Yuc104F as well
    call_category == "different_call" & zygosity_MCB == "hom_ref" ~ zygosity_Yuc104F,
    call_category == "different_call" & zygosity_Yuc104F == "hom_ref" ~ zygosity_MCB,
    call_category == "different_call" ~ if_else(zygosity_MCB != "hom_ref", zygosity_MCB, zygosity_Yuc104F),
    call_category == "ambiguous" ~ "ambiguous"
  ))




######################

### Script2. MCB-small-var-filter.R

wide <- as.data.table(wide)

gt_cols <- grep("^gt_GT_", names(wide), value = TRUE) # automatically find every column that starts with "gt_GT_"
suffs <- sub("^gt_GT_", "", gt_cols) # extract the sample names ("MCB", "Yuc104F", etc.)
metrics <- c("DP","AD_ALT","QUAL") # set filtering metrics

#### Set filtering thesholds ####
dp_threshold     <- 10
ad_threshold     <- 3
qual_thresholds <- list( SNP  = 30, INDEL = 10 )

#### Add Per-metric PASS flags & per-sample filter  ####
for(s in suffs) {
  # define metric‐column names
  dp_col   <- paste0("gt_DP_", s)
  ad_col   <- paste0("AD_ALT_",  s)
  qual_col <- paste0("QUAL_",    s)
  
  # PASS flags
  wide[, paste0("DP_pass_",  s) := get(dp_col)   >= dp_threshold]
  wide[, paste0("AD_pass_",  s) := get(ad_col)   >= ad_threshold]
  wide[, paste0("Q_pass_", s) := get(qual_col) >= fifelse(variant_type == "SNP", qual_thresholds$SNP, qual_thresholds$INDEL)]
  
  # combined sample-level filter PASS/FAIL flag [TRUE only if all filters are PASS]
  wide[, paste0(s, "_filter") := 
         get(paste0("DP_pass_",s)) &
         get(paste0("AD_pass_",s)) &
         get(paste0("Q_pass_",s))
  ]
}


gt_cols <- grep("^gt_GT_", names(wide), value = TRUE) # automatically find every column that starts with "gt_GT_"
suffs <- sub("^gt_GT_", "", gt_cols) # extract the sample names ("MCB", "Yuc104F", etc.)
metrics <- c("DP","AD_ALT","QUAL") # set filtering metrics

#### Mask GTs that fail filtering to reclassify into filt_call_category #### 
# masked GT columns: if filter=FALSE then "./.", else keep original genotype (gt_GT_<sample>_mod)
for(i in seq_along(gt_cols)) {
  col_gt    <- gt_cols[i]
  col_filt  <- paste0(suffs[i], "_filter")
  col_mod   <- paste0(col_gt, "_mod")
  
  wide[, (col_mod) := fifelse( get(col_filt), get(col_gt), "./.")]
}

# compute filtered info-status (<sample>_info_filt)
for(i in seq_along(suffs)) {
  suff       <- suffs[i]
  gtmod_col  <- paste0(gt_cols[i], "_mod")
  info_filt  <- paste0(suff, "_info_filt")
  
  wide[, (info_filt) := fifelse(
    get(gtmod_col) %in% c(NA, "./."), "none",
    fifelse(grepl("\\.", get(gtmod_col)), "partial", "full")
  )]
}


# re‐run classification logic on the masked GTs
col1_mod <- paste0(gt_cols[1], "_mod")
col2_mod <- paste0(gt_cols[2], "_mod")
info1_f  <- paste0(suffs[1], "_info_filt")
info2_f  <- paste0(suffs[2], "_info_filt")

wide[, filt_call_category := fcase(
  # shared: full/full identical or pure het
  (
    (get(col1_mod) == get(col2_mod) &
       get(info1_f) == "full" & get(info2_f) == "full")
  ) | (
    get(col1_mod) %in% c("0/1","1/0") & get(col2_mod) %in% c("0/1","1/0")
  ), "shared",
  
  # missing/partial vs alt
  (
    (get(col1_mod) %in% c(NA, "./.","./0","0/.") & grepl("1", get(col2_mod))) |
      (get(col2_mod) %in% c(NA, "./.","./0","0/.") & grepl("1", get(col1_mod)))
  ), "different_call",
  
  # ref vs any 1
  (
    (get(col1_mod) == "0/0" & grepl("1", get(col2_mod))) |
      (get(col2_mod) == "0/0" & grepl("1", get(col1_mod)))
  ), "different_call",
  
  # het vs hom-alt
  (
    (get(col1_mod) %in% c("0/1","1/0") & get(col2_mod) == "1/1") |
      (get(col2_mod) %in% c("0/1","1/0") & get(col1_mod) == "1/1")
  ), "different_call",
  
  # fallback
  default = "ambiguous"
)]


# Remove all columns with *pass* in them
pass_cols <- grep("_pass_", names(wide), value = TRUE)
wide[, (pass_cols) := NULL]
#thr_cols <- c("dp_thr","gq_thr","ar_min_thr","ar_max_thr","ad_thr","qual_thr")
#wide[, (thr_cols) := NULL]

### Define updated sample-specific variants ###
wide <- wide %>% 
  mutate(sample_specificity_mod = case_when(
    filt_call_category == "shared" ~ "shared",
    filt_call_category == "ambiguous" ~ "ambiguous",
    filt_call_category == "different_call" & str_detect(gt_GT_Yuc104F_mod, "1") & !str_detect(gt_GT_MCB_mod, "1") ~ "Yuc104F_specific",
    filt_call_category == "different_call" & !str_detect(gt_GT_Yuc104F_mod, "1") & str_detect(gt_GT_MCB_mod, "1") ~ "MCB_specific",
    filt_call_category == "different_call" & str_detect(gt_GT_Yuc104F_mod, "1") & str_detect(gt_GT_MCB_mod, "1") ~ "complex"
  ))

#### MCB-specific variant counts ####
# Select rows where Yuc104F is “0/0” or “./.” AND MCB is neither “0/0” nor “./.”

# set variable for MCB-specific post-filtering
wide <- wide %>%
  mutate(MCB_specific = gt_GT_Yuc104F_mod %in% c("0/0","./.", "0/.", "./0") &
           !gt_GT_MCB_mod %in% c("0/0","./.","0/.", "./0") &
           MCB_filter == TRUE)

sel_MCB_only <- wide[
  gt_GT_Yuc104F_mod %in% c("0/0","./.", "0/.", "./0") & 
    !gt_GT_MCB_mod  %in% c("0/0","./.","0/.", "./0") &
    `MCB_filter` == TRUE,
  .(variant_id)  
]

message("Number of variants present only in MCB: ", nrow(sel_MCB_only))

######################

### Script3. MCB-region-annotator.R

library(stringr)
# install.packages("BiocManager")
# BiocManager::install("GenomicRanges")
library(GenomicRanges)

### Read in BED files with genomic regions - Duroc
perv_bed <- fread("BEDfiles/Ssc11.1_HC_41PERV.bed", col.names = c("CHROM", "START", "END","pType", "type"))    
centro_bed <- fread("BEDfiles/Ssc11.centromere_and_like_regions.bed", col.names = c("CHROM", "START", "END","feature"))    
homopolymer_bed <- fread("BEDfiles/homopolymers_Sus_scrofa.Sscrofa11.1.dna.toplevel_min5bp_with_3bpflanks.bed", col.names = c("CHROM", "START", "END","hp_ID",  "V5", "V6"))[, .(CHROM, START, END, hp_ID)]   
homopolymer_bed[, CHROM := tstrsplit(CHROM, " ")[[1]]]
trf_bed <- fread("BEDfiles/susScr11.trf_mod_swaped.bed", col.names = c("CHROM", "START", "END","feature"))     
gaps_bed <- fread("BEDfiles/Ns_all_chromosomes_SusScrofa11.bed", col.names = c("CHROM", "START", "END","gap_ID", "V5", "V6"))[, .(CHROM, START, END, gap_ID)]
# repeats_bed <- fread("BEDfiles/nestedRepeats_swaped.bed", col.names = c("CHROM", "START", "END","name"))  

# Adjust to 1-based coordinates (for START only) to match VCF representation 
adjust_bed <- function(bed_df) {
  bed_df %>% mutate(START = START + 1) 
}

# Adjust BEDs before GRanges conversion
homopolymer_bed <- adjust_bed(homopolymer_bed)
perv_bed <- adjust_bed(perv_bed)
centro_bed <- adjust_bed(centro_bed)
trf_bed <- adjust_bed(trf_bed)
gaps_bed <- adjust_bed(gaps_bed)


#---------------------------------------------------#
##### Define Genomic Ranges for BED overlaps #####
#---------------------------------------------------#

# Convert variants into Genomic Ranges
variant_gr <- GRanges(
  seqnames = wide$CHROM,
  ranges = IRanges(
    start = wide$POS,
    end   = wide$POS + pmax(nchar(wide$REF), 1) - 1
  )
)


# Define overlap detection function
annotate_overlap <- function(gr_variants, bed_df) {
  bed_gr <- GRanges(seqnames = bed_df$CHROM, ranges = IRanges(start = bed_df$START, end = bed_df$END))
  overlapsAny(gr_variants, bed_gr)
}


# Annotate with logical column for each feature type
wide$in_homopolymer <- annotate_overlap(variant_gr, homopolymer_bed)
wide$in_perv <- annotate_overlap(variant_gr, perv_bed)
wide$in_centromere  <- annotate_overlap(variant_gr, centro_bed)
wide$in_trf  <- annotate_overlap(variant_gr, trf_bed)
# wide$in_repeat      <- annotate_overlap(variant_gr, repeats_bed)
wide$in_gap         <- annotate_overlap(variant_gr, gaps_bed)


#------------------------#
##### Summary Counts #####
#------------------------#

feature_cols <- c("in_homopolymer", "in_centromere", "in_gap", "in_trf", "in_perv")  # "in_repeat", 

# Create "total" column that == TRUE if any feature is TRUE
wide <- wide %>%
  mutate(in_any_feature = if_any(all_of(feature_cols), ~ .x))



#---------------------------------#
##### Genomic feature overlap : Exons
#---------------------------------#

# BiocManager::install("rtracklayer")
library(rtracklayer)

#### Read in GTF file ####
gtf <- import("BEDfiles/Sus_scrofa.Sscrofa11.1.106.gtf.gz")

# extract all Exons and UTRs from GTF, in protein coding (pc) regions
exons <- gtf[gtf$type == "exon"]
utrs  <- gtf[gtf$type %in% c("five_prime_utr", "three_prime_utr")]
exons_pc <- exons[mcols(exons)$gene_biotype == "protein_coding"]
utrs_pc  <- utrs[mcols(utrs)$gene_biotype == "protein_coding"]

# Find overlaps
wide$in_exon     <- overlapsAny(variant_gr, exons)
wide$in_utr      <- overlapsAny(variant_gr, utrs)
wide$in_pc_exon  <- overlapsAny(variant_gr, exons_pc)
wide$in_pc_utr   <- overlapsAny(variant_gr, utrs_pc)

# Define genomic feature columns
gen_feature_cols <- c("in_exon", "in_utr")
gen_pc_feature_cols <- c("in_pc_exon", "in_pc_utr")

# Create "total" column that == TRUE if any feature is TRUE
wide <- wide %>%
  mutate(in_exon_or_utr = if_any(all_of(gen_feature_cols), ~ .x))
wide <- wide %>%
  mutate(in_exon_or_utr_pc = if_any(all_of(gen_pc_feature_cols), ~ .x))

#### GTF-informed MCB-specific variant counts ####
sel_MCB_only <- wide[
  gt_GT_Yuc104F_mod %in% c("0/0","./.", "0/.", "./0") & 
    !gt_GT_MCB_mod  %in% c("0/0","./.","0/.", "./0") &
    MCB_specific == TRUE &
    in_exon_or_utr == TRUE &
    in_any_feature == FALSE,
  .(variant_id)  
]
message("Number of variants present only in MCB, overlapping exon/UTR regions, filtered by qual and region filters: ", nrow(sel_MCB_only))

sel_MCB_only <- wide[
  gt_GT_Yuc104F_mod %in% c("0/0","./.", "0/.", "./0") & 
    !gt_GT_MCB_mod  %in% c("0/0","./.","0/.", "./0") &
    MCB_specific == TRUE &
    in_exon_or_utr_pc == TRUE &
    in_any_feature == FALSE,
  .(variant_id)  
]

message("Number of variants present only in MCB, overlapping exon/UTR regions of protein coding genes, filtered by qual and region filters: ", nrow(sel_MCB_only))

#---------------------------------#
##### Genomic feature overlap : Introns
#---------------------------------#

BiocManager::install("GenomicFeatures")
library(GenomicFeatures)
BiocManager::install("txdbmaker")
library(txdbmaker)

txdb <- makeTxDbFromGFF("BEDfiles/Sus_scrofa.Sscrofa11.1.106.gtf.gz", format = "gtf")

# Extract introns by transcript, and then reduce to unique ranges
introns_by_tx <- intronsByTranscript(txdb, use.names = TRUE)
introns_all   <- reduce(unlist(introns_by_tx, use.names = FALSE))

# Get transcript annotations with biotype info
tx2gene <- unique(gtf[gtf$type == "transcript", c("transcript_id", "gene_biotype")])

# Keep only introns from protein-coding genes
tx_ids_pc <- tx2gene$transcript_id[tx2gene$gene_biotype == "protein_coding"]
introns_pc <- introns_by_tx[names(introns_by_tx) %in% tx_ids_pc]
introns_pc_all <- reduce(unlist(introns_pc, use.names = FALSE))

# Add overlap columns
wide$in_intron_pc <- overlapsAny(variant_gr, introns_pc_all)

# Create a combined "non-exonic" feature column
wide <- wide %>%
  mutate(in_nonexonic_pc = !in_pc_exon & !in_pc_utr & in_intron_pc)


fwrite(wide, file = "DeepVariant_MCB_MasterTable.csv") 

