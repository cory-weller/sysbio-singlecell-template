#!/usr/bin/env bash

OUTDIR=${1}

chrs=( chr{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X} )
pre='/fdb/1000genomes/20201028_3202_raw_GT_with_annot/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_'
post='.recalibrated_variants.annotated.vcf.gz'

module purge
module load samtools/1.23

cd $TMPDIR

for chr in ${chrs[@]}; do
    vcf_all="${pre}${chr}${post}"
    echo $vcf_all
    coding_variants=${vcf_all%.vcf.gz}.coding.txt
    echo $coding_variants

    ln -s $vcf_all ${chr}.vcf.gz
    if [ ! -f $vcf_all.tbi ]; then
        tabix -p vcf ${chr}.vcf.gz
    else
        ln -s $vcf_all.tbi ${chr}.vcf.gz.tbi
    fi
    zcat ${chr}.vcf.gz | head -n 8000 | egrep '^#' > ${chr}_raw.vcf
    tabix --regions <(awk -v OFS="\t" 'NR>1 {print $1,$2,$2}' $coding_variants) ${chr}.vcf.gz >> ${chr}_raw.vcf
    bcftools sort ${chr}_raw.vcf | bgzip -c - > ${chr}_sorted.vcf.gz
    tabix -p vcf ${chr}_sorted.vcf.gz
    bcftools view --types snps ${chr}_sorted.vcf.gz| bgzip -c - > ${chr}_snps.vcf.gz
    bcftools annotate -x '^INFO/AC,^INFO/AF,^INFO/AN' ${chr}_snps.vcf.gz | bgzip -c - > ${chr}_AF.vcf.gz
    tabix -p vcf ${chr}_AF.vcf.gz
    rm ${chr}_raw.vcf ${chr}_snps.vcf.gz ${chr}_sorted.vcf.gz ${chr}_sorted.vcf.gz.tbi ${chr}.vcf.gz ${chr}.vcf.gz.tbi
done

bcftools concat -o ALL.vcf.gz \
    chr1_AF.vcf.gz \
    chr10_AF.vcf.gz \
    chr11_AF.vcf.gz \
    chr12_AF.vcf.gz \
    chr13_AF.vcf.gz \
    chr14_AF.vcf.gz \
    chr15_AF.vcf.gz \
    chr16_AF.vcf.gz \
    chr17_AF.vcf.gz \
    chr18_AF.vcf.gz \
    chr19_AF.vcf.gz \
    chr2_AF.vcf.gz \
    chr20_AF.vcf.gz \
    chr21_AF.vcf.gz \
    chr22_AF.vcf.gz \
    chr3_AF.vcf.gz \
    chr4_AF.vcf.gz \
    chr5_AF.vcf.gz \
    chr6_AF.vcf.gz \
    chr7_AF.vcf.gz \
    chr8_AF.vcf.gz \
    chr9_AF.vcf.gz \
    chrX_AF.vcf.gz

header='##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##ALT=<ID=NON_REF,Description="Represents any possible alternative allele at this location">
##FILTER=<ID=LowQual,Description="Low quality">
##FILTER=<ID=VQSRTrancheINDEL99.00to100.00+,Description="Truth sensitivity tranche level for INDEL model at VQS Lod < -85077.1808">
##FILTER=<ID=VQSRTrancheINDEL99.00to100.00,Description="Truth sensitivity tranche level for INDEL model at VQS Lod: -85077.1808 <= x < -2.3474">
##FILTER=<ID=VQSRTrancheSNP99.80to100.00+,Description="Truth sensitivity tranche level for SNP model at VQS Lod < -343144.6144">
##FILTER=<ID=VQSRTrancheSNP99.80to100.00,Description="Truth sensitivity tranche level for SNP model at VQS Lod: -343144.6144 <= x < -13.4687">
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency, for each ALT allele, in the same order as listed">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
##contig=<ID=chr1,length=248956422>
##contig=<ID=chr10,length=133797422>
##contig=<ID=chr11,length=135086622>
##contig=<ID=chr12,length=133275309>
##contig=<ID=chr13,length=114364328>
##contig=<ID=chr14,length=107043718>
##contig=<ID=chr15,length=101991189>
##contig=<ID=chr16,length=90338345>
##contig=<ID=chr17,length=83257441>
##contig=<ID=chr18,length=80373285>
##contig=<ID=chr19,length=58617616>
##contig=<ID=chr2,length=242193529>
##contig=<ID=chr20,length=64444167>
##contig=<ID=chr21,length=46709983>
##contig=<ID=chr22,length=50818468>
##contig=<ID=chr3,length=198295559>
##contig=<ID=chr4,length=190214555>
##contig=<ID=chr5,length=181538259>
##contig=<ID=chr6,length=170805979>
##contig=<ID=chr7,length=159345973>
##contig=<ID=chr8,length=145138636>
##contig=<ID=chr9,length=138394717>
##contig=<ID=chrX,length=156040895>
##reference=file:///gpfs/commons/datasets/old-nygc-resources/GRCh38_1000genomes/GRCh38_full_analysis_set_plus_decoy_hla.fa
##bcftools_normVersion=1.3.1-98-ga6a7829+htslib-1.3.1-64-g74bcfd7
##bcftools_normCommand=norm -f /gpfs/commons/datasets/old-nygc-resources/GRCh38_1000genomes/GRCh38_full_analysis_set_plus_decoy_hla.fa -m-both /gpfs/commons/groups/compbio/projects/1KGP_3202_SNV_INDEL_phasing/flagging/HWE_ME_final/flagged_vcf/CCDG_14151_B01_GRM_WGS_2020-08-05_chr1.recalibrated_variants.flagged.vcf.gz; Date=Mon Sep 28 11:34:23 2020
##bcftools_viewVersion=1.3.1-98-ga6a7829+htslib-1.3.1-64-g74bcfd7
##bcftools_viewCommand=view -G -O z -o /gpfs/commons/home/usevani/compbio/CCDG/Project_CCDG_14151_B01_GRM_WGS/final_annotated_vcfs/tmp_dir_annt/CCDG_14151_B01_GRM_WGS_2020-08-05_chr1.recalibrated_variants.annotated.normalize.vcf.gz -; Date=Mon Sep 28 11:34:23 2020
##SnpEffVersion="4.3i (build 2016-12-15 22:33), by Pablo Cingolani"
##SnpEffCmd="SnpEff  -noStats -lof GRCh38.86 /gpfs/commons/home/usevani/compbio/CCDG/Project_CCDG_14151_B01_GRM_WGS/final_annotated_vcfs/tmp_dir_annt/CCDG_14151_B01_GRM_WGS_2020-08-05_chr1.recalibrated_variants.annotated.normalize.vcf.gz "
##bcftools_viewCommand=view /gpfs/commons/home/usevani/compbio/CCDG/Project_CCDG_14151_B01_GRM_WGS/final_annotated_vcfs/tmp_dir_annt/CCDG_14151_B01_GRM_WGS_2020-08-05_chr1.recalibrated_variants.annotated.repeats.vcf.gz; Date=Mon Sep 28 21:16:47 2020
##bcftools_viewVersion=1.23+htslib-1.23
##bcftools_viewCommand=view --types snps chr1_sorted.vcf.gz; Date=Tue Jul 28 10:00:51 2026
##bcftools_annotateVersion=1.23+htslib-1.23
##bcftools_annotateCommand=annotate -x ^INFO/AC,^INFO/AF,^INFO/AN chr1_snps.vcf.gz; Date=Tue Jul 28 10:00:53 2026
##bcftools_concatVersion=1.23+htslib-1.23
##bcftools_concatCommand=concat -o ALL.vcf.gz chr1_AF.vcf.gz chr10_AF.vcf.gz chr11_AF.vcf.gz chr12_AF.vcf.gz chr13_AF.vcf.gz chr14_AF.vcf.gz chr15_AF.vcf.gz chr16_AF.vcf.gz chr17_AF.vcf.gz chr18_AF.vcf.gz chr19_AF.vcf.gz chr2_AF.vcf.gz chr20_AF.vcf.gz chr21_AF.vcf.gz chr22_AF.vcf.gz chr3_AF.vcf.gz chr4_AF.vcf.gz chr5_AF.vcf.gz chr6_AF.vcf.gz chr7_AF.vcf.gz chr8_AF.vcf.gz chr9_AF.vcf.gz chrX_AF.vcf.gz; Date=Tue Jul 28 10:36:56 2026
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO'

echo -e "$header" > reorder.vcf
zgrep -e '^[^#]'  ALL.vcf.gz >> reorder.vcf
bgzip reorder.vcf

mv reorder.vcf.gz dsc-pileup-ref.vcf.gz
tabix -p vcf dsc-pileup-ref.vcf.gz

mv dsc-pileup-ref.vcf.gz $OUTDIR
mv dsc-pileup-ref.vcf.gz.tbi $OUTDIR
