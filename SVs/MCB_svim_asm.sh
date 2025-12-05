# genome assemblies of MCB are completed on SBG
#Transfer to /data5/hifi_assembly/GC64/GC64-23_008_MCB: https://igor.sbgenomics.com/u/egenesis/pacbio-hifi-assembly/tasks/e55bd492-6f15-4301-b1d3-5d0049f83b2f/

# Run SVIM-asm for GC64-23_008_MCB_FSDC hifiasm hap1 and hap2 as well as Ssc11.1
REF=/data5/sv/ref/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa
HAP1="/data5/hifi_assembly/GC64/GC64-23_008_MCB/GC64-23_008_MCB_FSDC.asm.bp.hap1.p_ctg.fa"
HAP2="/data5/hifi_assembly/GC64/GC64-23_008_MCB/GC64-23_008_MCB_FSDC.asm.bp.hap2.p_ctg.fa"
NAME1=$(basename $HAP1 .fa)
NAME2=$(basename $HAP2 .fa)
THREADS=32

mkdir /data5/sv/GC64-23_008_MCB/svim_asm
working_dir="/data5/sv/GC64-23_008_MCB/svim_asm"
prefix=GC64-23_008_MCB_FSDC_hap1_hap2_to_Ssc11
cd ${working_dir}

# HAP1
mamba activate svimasm_env
minimap2 -a -x asm5 --cs -r2k -t $THREADS $REF $HAP1 > ${NAME1}_to_Ssc11_hap1.sam

conda activate ont
samtools sort -m4G -@4 -o ${NAME1}_to_Ssc11_hap1.sorted.bam ${NAME1}_to_Ssc11_hap1.sam
samtools index ${NAME1}_to_Ssc11_hap1.sorted.bam


# HAP2
mamba activate svimasm_env
minimap2 -a -x asm5 --cs -r2k -t $THREADS $REF $HAP2 > ${NAME2}_to_Ssc11_hap2.sam
conda activate ont

conda activate ont
samtools sort -m4G -@4 -o ${NAME2}_to_Ssc11_hap2.sorted.bam ${NAME2}_to_Ssc11_hap2.sam
samtools index ${NAME2}_to_Ssc11_hap2.sorted.bam

#--> alignment of Hap2 and Hap1 takes about 1hr to finish each

mamba activate svimasm_env
time svim-asm diploid  ${working_dir} ${NAME1}_to_Ssc11_hap1.sorted.bam ${NAME2}_to_Ssc11_hap2.sorted.bam $REF --sample GC64-23_008_MCB_FSDC

# rename vcf
prefix=GC64-23_008_MCB_FSDC_hap1_hap2_to_Ssc11
mv variants.vcf ${prefix}.vcf

#!!!! This following step is not necessary if the result will be merged with another sample. Can stay as single sample VCF if it's going to be multi-caller, multi-sample VCF.
# Multi-sample VCF between Yuc104F and GC64-23_008_MCB_FSDC
conda activate ont_ass

cd /data5/sv/svim_asm
# multi-sample VCF
VCF1="/data5/sv/Yuc104F.30X/svim_asm/Yuc104F_hap1_hap2_to_Ssc11.vcf"
VCF2="/data5/sv/GC64-23_008_MCB/svim_asm/GC64-23_008_MCB_FSDC_hap1_hap2_to_Ssc11.vcf"
prefix="Yuc104F_GC64-23_008_MCB_FSDC_to_Ssc11.1_svim_asm"

# min SV size 40bp because svim-asm minSV = 40bp
surpyvor merge --variants $VCF1 $VCF2 \
 -o ${prefix}.vcf \
 -d 1000 -l 40 -c 1 -s

#  format vcf into simple format to filter based on diff genotypes
conda activate ont_ass 
cd /data5/sv/svim_asm
prefix="Yuc104F_GC64-23_008_MCB_FSDC_to_Ssc11.1_svim_asm"

# modify VCF to tsv to filter out by Yuc104F and MCB genotype
# bcftools query -f'%CHROM %POS %SVTYPE %ID %SVLEN %END[ %GT]\n' ${prefix}.vcf |  awk '{if ($NF == $(NF-1)) next} 1' - |  awk '{if (($NF == "0/1" && $(NF-1) == "1/0") || ($NF == "1/0" && $(NF-1) == "0/1")) next} 1' - >  ${prefix}.diffgenotype.tsv


# the chr6 26kb del in OGM, two 20kb inv in PBSV. In svim-asm, one inv is called
# 6 77074997 INV svim_asm.INV.44 20452 77095449 ./. 1/0


