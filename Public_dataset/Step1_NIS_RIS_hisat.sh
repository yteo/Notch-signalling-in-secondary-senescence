module load hisat2 
module load samtools
# fastq file obtained from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1861895 (NIS and RIS)
# align single end reads NIS,RIS to hg19 using hisat2
for file in *_1.fastq.gz
do
	hisat2 -q -x /users/nneretti/data/annotation/hg19/hg19  \
	-U $file1 -S ${file%_1.fastq.gz}.sam -p 20 --dta-cufflinks

	samtools view -bS ${file%_1.fastq.gz}.sam |samtools sort -n - -o \
	${file%_1.fastq.gz}_sorted
	rm ${file%_1.fastq.gz}.sam


	# htseq count for bam files


	python htseq-count --stranded=no -f bam ${file%_1.fastq.gz}_sorted \
	-r name /data/gencode.v16.annotation.gtf \
	> ${file%_1.fastq.gz}_sorted_htseq.txt 

done
