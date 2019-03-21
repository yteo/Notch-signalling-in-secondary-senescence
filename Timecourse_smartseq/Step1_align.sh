module load samtools
module load hisat2/2.0.1-beta
module load trimgalore/0.4.0
module load cutadapt/1.3


for sample in *_1.fastq.gz
do
	trim_galore $sample ${sample%_1.fastq.gz}_2.fastq.gz --paired --path_to_cutadapt /gpfs/runtime/opt/cutadapt/1.3/bin/cutadapt
done


for sample in *_1_val_1.fq.gz
do
	file=$(echo ${sample%_1_val_1.fq.gz})
	python Step1.py -infile1 ${sample} -infile2 ${sample%_1_val_1.fq.gz}_2_val_2.fq.gz -outputprefix ${file} -o ${file} -t 20
done


# cat all mappingstat for each cell together to create mappedreadsstat.txt
