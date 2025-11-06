# Oct 16, 2025
# Generate multi-caller multi-sample VCF (union of 2 callers, pbsv and SVIM-asm)
############################################################################################################
# STEP 1: merge all vcfs for the same sample into 1 vcf per sample with SURVIVOR merge. Called by at least 2 out of 2 callers. MULTI-CALLER
# IP42
cd /data5/sv/sam_SV/Duroc_MCB/

vim Yuc104F_2callers.vcf.list
pbsv/Yuc104F_toSsc11.1_pbsv.vcf
svim_asm/Yuc104F_toSsc11.1_svim_asm.vcf

vim MCB_GC64-23_008_2callers.vcf.list
pbsv/GC64-23_008_MCB_FSDC_toSsc11.1_pbsv.vcf
svim_asm/GC64-23_008_MCB_FSDC_toSsc11.1_svim_asm.vcf

mamba activate survivor
# Callset: same SV type, no strand required, at least called by 2 out of 3 callers per sample, sawfish is not considered
cd /data5/sv/sam_SV/Duroc_MCB

mkdir union_pbsv_svim-asm

SURVIVOR merge Yuc104F_2callers.vcf.list 1000 1 1 0 0 20 union_pbsv_svim-asm/Yuc104F_to_Ssc11.1_pbsv_svimAsm_Survivor1100.vcf
merging entries: 106616
merging entries: 53145

SURVIVOR merge MCB_GC64-23_008_2callers.vcf.list 1000 1 1 0 0 20 union_pbsv_svim-asm/MCB_to_Ssc11.1_pbsv_svimAsm_Survivor1100.vcf
merging entries: 106616
merging entries: 53145

############################################################################################################
#STEP 2: combine the 1 vcf per sample for all samples with SURVIVOR merge. MULTI-CALLER and MULTI-SAMPLE
# at least one sample calls it report it (Minimum number of supporting caller: 1)
cd /data5/sv/sam_SV/Duroc_MCB/union_pbsv_svim-asm

vim Yuc104F_MCB_union_pbsv_svim-asm.list
Yuc104F_to_Ssc11.1_pbsv_svimAsm_Survivor1100.vcf
MCB_to_Ssc11.1_pbsv_svimAsm_Survivor1100.vcf


SURVIVOR merge Yuc104F_MCB_union_pbsv_svim-asm.list 1000 1 1 0 0 20 \
Yuc104F_MCB_to_Ssc11.1_pbsv_svimAsm_Survivor1100_union.vcf
merging entries: 100290
merging entries: 100290

#Check if known SVs are called in MCB from this merged set
# download MCB BED from SBG: https://igor.sbgenomics.com/u/egenesis/pacbio-iga-aberrantinsertionofpayload/files/6706dc1770d8a41dd6e4902a/
BED="/data5/sv/sam_SV/Duroc/all_OGM_sv_GC64-23_rmCentro_050524.bed"
VCF="/data5/sv/sam_SV/Duroc_MCB/union_pbsv_svim-asm/Yuc104F_MCB_to_Ssc11.1_pbsv_svimAsm_Survivor1100_union.vcf"


cd /data5/sv/sam_SV/Duroc_MCB/union_pbsv_svim-asm
conda activate ont_ass
bedtools intersect -a $VCF -b $BED -header > Yuc104F_MCB_all_OGM_sv_GC64-23_rmCentro_050524_pbsv_svimAsm_union.vcf

bcftools query -f'%CHROM %POS %SVTYPE %ID %SVLEN %END[ %GT]\n' Yuc104F_MCB_all_OGM_sv_GC64-23_rmCentro_050524_pbsv_svimAsm_union.vcf | awk '{if ($7 == $8) next} 1' | awk '{if ($7 == "0|1" && $8 == "1|0") next} 1' | awk '{if ($8 == "0|1" && $7 == "1|0") next} 1' | column -t > Yuc104F_MCB_all_OGM_sv_GC64-23_rmCentro_050524_pbsv_svimAsm_union.filtered.diffgeno.tsv

