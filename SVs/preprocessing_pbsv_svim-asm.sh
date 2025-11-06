# Process svim-asm and pbsv files to be used as SURVIVOR input
## Note this was originally processed for 4 SV callers, pbsv, svim-asm, sniffles and sawfish. But can just keep pbsv and svim-asm for the final result
# MCB pbsv result: https://igor.sbgenomics.com/u/egenesis/pacbio-iga-aberrantinsertionofpayload/tasks/a0b312d9-098d-471b-825e-a3a278147479/
# Convert each SV caller result to single-sample VCF per caller
(base) ubuntu@ip-192-168-10-42:/data5/sv/sam_SV/Duroc_MCB$ 

cat Yuc104F_MCB_GC64-23_008_multi-sample_perTool.list
pbsv/GC64-23_008_MCB_FSDC_Yuc104F_pbsv_on_winnowmap2.vcf
sawfish/Yuc104F_GC64-23_008_MCB_joint_call_sawfish.vcf
sniffles/Yuc104F_GC64-23_008_MCB_map2_Ssc11_v260.snf.vcf
svim_asm/Yuc104F_GC64-23_008_MCB_FSDC_to_Ssc11.1_svim_asm.vcf

# split joint-calls per tool into per caller per vcf per sample

conda activate ont_ass
cd /data5/sv/sam_SV/Duroc_MCB

for VCF in `cat Yuc104F_MCB_GC64-23_008_multi-sample_perTool.list`
do
    for SAMPLE in $(bcftools query -l $VCF)
    do
        bcftools view -s $SAMPLE -o $(dirname $VCF)/${SAMPLE}_toSsc11.1_$(dirname $VCF).vcf -O v $VCF
    done
done

# check results
(ont_ass) ubuntu@ip-192-168-10-42:/data5/sv/sam_SV/Duroc_MCB$ 
 ls -1 */Yuc104F_toSsc11.1*.vcf > Yuc104F_4callers.vcf.list
 
cat Yuc104F_4callers.vcf.list
pbsv/Yuc104F_toSsc11.1_pbsv.vcf
sawfish/Yuc104F_toSsc11.1_sawfish.vcf
sniffles/Yuc104F_toSsc11.1_sniffles.vcf
svim_asm/Yuc104F_toSsc11.1_svim_asm.vcf

 ls -1 */GC64-23_008_MCB_FSDC_toSsc11.1*.vcf > MCB_GC64-23_008_4callers.vcf.list
 pbsv/GC64-23_008_MCB_FSDC_toSsc11.1_pbsv.vcf
sawfish/GC64-23_008_MCB_FSDC_toSsc11.1_sawfish.vcf
svim_asm/GC64-23_008_MCB_FSDC_toSsc11.1_svim_asm.vcf
sniffles/GC64-23_008_MCB_toSsc11.1_sniffles.vcf

