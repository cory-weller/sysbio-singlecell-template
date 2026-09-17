#!/usr/bin/env bash

export BAM=${1}
export OUTNAME=${2}
export REF_FASTA=${3}
export THREADS=${4}

cd $TMPDIR

mkdir -p $(dirname ${OUTNAME})


pileup() {
    contig=$1
    bcftools mpileup -r ${contig} -d 8000 -f ${REF_FASTA} ${BAM} | bcftools call --ploidy GRCh38 -v -m -Oz -o ${contig}.vcf.gz
}

export -f pileup

if [ ! -f ${BAM}.bai ]; then
    samtools index --threads ${THREADS} ${BAM}
fi

# Build array of chr to iterate over with parallel
chrs=( chr{1..22} chrX )

# Execute as parallel
parallel -j ${THREADS} pileup ::: ${chrs[@]}

if [ ! -f chr1.vcf.gz ]; then
    echo "ERROR: VCF file does not exist from bcftools mpileup"
    exit 1
fi


# Build array of files to concatenate
declare -a files
chr_ordered=(chr1 chr{10..19} chr2 chr{20..22} chr{3..9} chrX)
for i in ${chr_ordered[@]}; do
    files+=("${i}.vcf.gz")
done

bcftools concat ${files[@]} -Oz -o ALL.vcf.gz && \
bcftools +fill-tags ALL.vcf.gz -Oz -o ALL_AF.vcf.gz -- -t AF && \
rm ALL.vcf.gz && \
rm chr*.vcf.gz

cp ALL_AF.vcf.gz ${OUTNAME}