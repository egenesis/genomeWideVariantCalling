# Perform Yuc104F svim-asm SV calling against Ssc11.1

# install necessary tools for svim-asm
mamba create -n svimasm_env --channel bioconda svim-asm minimap2

 mamba activate svimasm_env
mamba install samtools

# call SV using assembly-based approach.
# INPUT:
# Yuc104F hifiasm hap1 and hap2 as well as Ssc11.1
mamba activate svimasm_env
REF=/data5/sv/ref/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa
# haplotype assemblies are downloaded from SBG: https://igor.sbgenomics.com/u/egenesis/pacbio-hifi-assembly/tasks/fca105dd-4e74-442d-8f54-4c3e4218ef79/
HAP1="/data5/hifi_assembly/GC64/30X_assembly/Yuc104F.asm.bp.hap1.p_ctg.fa"
HAP2="/data5/hifi_assembly/GC64/30X_assembly/Yuc104F.asm.bp.hap2.p_ctg.fa"
NAME1=$(basename $HAP1 .fa)
NAME2=$(basename $HAP2 .fa)
THREADS=32
mkdir /data5/sv/Yuc104F.30X/svim_asm
working_dir="/data5/sv/Yuc104F.30X/svim_asm"
prefix=Yuc014F_hap1_hap2_to_Ssc11

cd ${working_dir}

minimap2 -a -x asm5 --cs -r2k -t $THREADS $REF $HAP1 > ${NAME1}_to_Ssc11_hap1.sam
minimap2 -a -x asm5 --cs -r2k -t $THREADS $REF $HAP2 > ${NAME2}_to_Ssc11_hap2.sam

conda activate ont
samtools sort -m4G -@4 -o ${NAME1}_to_Ssc11_hap1.sorted.bam ${NAME1}_to_Ssc11_hap1.sam
samtools index ${NAME1}_to_Ssc11_hap1.sorted.bam

samtools sort -m4G -@4 -o ${NAME2}_to_Ssc11_hap2.sorted.bam ${NAME2}_to_Ssc11_hap2.sam
samtools index ${NAME2}_to_Ssc11_hap2.sorted.bam

mamba activate svimasm_env
time svim-asm diploid ${working_dir} ${NAME1}_to_Ssc11_hap1.sorted.bam ${NAME2}_to_Ssc11_hap2.sorted.bam $REF  --sample Yuc104F

# rename vcf
prefix=Yuc104F_hap1_hap2_to_Ssc11
mv variants.vcf ${prefix}.vcf
