HISAT2_INDEXES=/home/MOKA/antonio/RNA-SeqTools/indexes/grch38/
FASTQ_FILES_R1=$(ls *_1_val_1.fq.gz |cut -d "_" -f1)
FASTQ_FILES_R2=$(ls *_2_val_2.fq.gz |cut -d "_" -f1)

for i in $FASTQ_FILES_R1; do hisat2 -q -p 20 -x $HISAT2_INDEXES/genome -1 ${i}_1_val_1.fq.gz -2 ${i}_2_val_2.fq.gz | samtools view -Sbh > ${i}.bam; done
