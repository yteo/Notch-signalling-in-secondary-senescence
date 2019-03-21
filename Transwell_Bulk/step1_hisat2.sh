module load samtools

for file in *_R1_001.fastq.gz
	con=$(echo $file|cut -d "_" -f1)
	file1=${file}
	file2=${file%_R1_001.fastq.gz}_R2_001.fastq.gz
	/gpfs/scratch/yteo/hisat2/bin/hisat2 -q -x /data/hg19 \
	-1 $file1 -2 $file2 -S ${con}.sam -p 20
	samtools view -bS ${con}.sam |samtools sort -n - -o ${con}_sorted.bam
	rm ${con}.sam

	python htseq-count -f bam --stranded=no \
	${con}_sorted.bam -r name \
	/data/hg19.annotation.gtf > ${con}_htseq.txt

done
