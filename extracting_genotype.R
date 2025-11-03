library(data.table)
library(dplyr)


setwd("~/Library/CloudStorage/OneDrive-eGenesisBio/Computational Biology Team/Sanaz/GenomeWideVariantCalling/Rscripts")

clair3 <- fread("Results/Clair3_MCB_MasterTable_simplified.csv")
deep <- fread("Results/DeepVariant_MCB_MasterTable_simplified.csv")

varfile <- fread("Appendix_Table_for_small_variant_and_SV_manual_curation.csv", header=TRUE)

result_gt <- data.frame()
result_gt$Uploaded_ID <- file$Uploaded_ID


matches <- rbind(
  clair3[variant_id %in% varfile$Uploaded_ID, 
         .(Uploaded_ID = variant_id, gt_GT_MCB_mod, gt_GT_Yuc104F_mod)],
  
  deep[variant_id %in% varfile$Uploaded_ID, 
       .(Uploaded_ID = variant_id, gt_GT_MCB_mod, gt_GT_Yuc104F_mod)]
)

matches_filtered <- matches[grepl("1", gt_GT_MCB_mod) & !grepl("1", gt_GT_Yuc104F_mod)]
matches_unique <- matches_filtered[!duplicated(Uploaded_ID)]

matches_ordered <- matches_unique[match(varfile$Uploaded_ID, Uploaded_ID)]

fwrite(matches_ordered, "genotype_MCB_small_variants.csv")