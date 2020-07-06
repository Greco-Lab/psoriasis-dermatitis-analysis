HISAT2_INDEXES=/home/MOKA/antonio/RNA-SeqTools/indexes/grch38/
FASTQ_FILES=$(ls *.fastq.gz |cut -d "." -f1)

for i in $FASTQ_FILES; do hisat2 -q -p 20 -x $HISAT2_INDEXES/genome -U ${i}.fastq.gz | samtools view -Sbh > ${i}.bam; done
