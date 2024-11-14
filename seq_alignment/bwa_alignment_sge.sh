#!/bin/bash

#$ -N bwa              # Job name
# -o bwa_logs/output  # Path for standard output
# -e bwa_logs/error   # Path for standard error
#$ -l h_rt=04:00:00    # Requested runtime (hh:mm:ss)
#$ -l h_vmem=16G       # Memory per slot
#$ -cwd                # Run job from the current working directory
#$ -V                  # Export environment variables to the job

PATH=$HOME/nobackup/samtools/bin/:$PATH
PATH=$HOME/nobackup/bcftools/bin/:$PATH
PATH=$HOME/nobackup/software/bwa/:$PATH

BWA_EXEC="bwa mem"

FASTA="/home/home01/scrb/nobackup/genomes/reference/submission.contigs.fasta"
FASTAGZ=${FASTA}.gz
FASTQ_DIR="/home/home01/scrb/nobackup/genomes/sequences/${FASTQ_DIR}"

# FASTQ="/home/home01/scrb/nobackup/genomes/sequences/"
# "/home/rodrigo/Postdoc/GWAS/data/genomes/sequences/01/1_L4_1.fq.gz"

OUTPUT_FOLDER="/home/home01/scrb/nobackup/genomes/alignments_paired_end"

# Check if FASTQ_DIR is set as an environment variable
if [ -z "$FASTQ_DIR" ]; then
    echo "Error: FASTQ_DIR environment variable is not set."
    exit 1
fi

for FASTQ_PREFIX in $(for FILE in `ls ${FASTQ_DIR}/*fq.gz`; do echo ${FILE%_*fq.gz}; done | uniq); do
   
    echo "$FASTQ_PREFIX"
    PREFIX=$(basename $FASTQ_PREFIX)
    PREFIX=${OUTPUT_FOLDER}/${PREFIX}
    SAMFILE=${PREFIX}.sam.gz
    
    echo "Prefix is ${PREFIX}"
    FASTQ_1=${FASTQ_PREFIX}_1.fq.gz
    FASTQ_2=${FASTQ_PREFIX}_2.fq.gz

    echo $SAMFILE

    # Check if the SAM file already exists
    if [ -f "$SAMFILE" ]; then
        echo "File $SAMFILE already exists. Skipping job for $FASTQ_1 and $FASTQ_2..."
    else
        echo "Processing $FASTQ_1..."
        BWA_CMD='bwa mem "$FASTAGZ" "${FASTQ_1}" "${FASTQ_2}" | gzip -3 > "$SAMFILE"'
        echo $BWA_CMD
        bwa mem "$FASTAGZ" "${FASTQ_1}" "${FASTQ_2}" | gzip -3 > "$SAMFILE"
        # $BWA_CMD
    fi

    if [ -f "${PREFIX}-sorted.bam" ]; then
        echo "BAM ${PREFIX}.bam already sorted and indexed. Skipping..."
    else
        samtools view -S -b ${PREFIX}.sam.gz > ${PREFIX}.bam
        samtools sort ${PREFIX}.bam -o ${PREFIX}-sorted.bam
        samtools index ${PREFIX}-sorted.bam
    fi

    bcftools mpileup -f $FASTA -d 20000 ${PREFIX}-sorted.bam | bcftools call --ploidy 1 -mv -Oz -o ${PREFIX}.vcf.gz

done
