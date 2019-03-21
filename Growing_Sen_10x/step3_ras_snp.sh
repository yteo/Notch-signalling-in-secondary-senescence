module load samtools/1.2
module load freebayes/2014Dec15
# freebayes
for file in *_merged_2_0_1/
	do 
	cd /data/$file/outs
	samtools view -h possorted_genome_bam.bam "chr11:534287-534289" > ${file%/}_ras.sam
	mkdir Ras
	mv ${file%/}_ras.sam ./Ras/
	cd Ras
	samtools view -bS ${file%/}_ras.sam > ${file%/}_ras.bam
	samtools index ${file%/}_ras.bam
	bamtools split -tag CB -in ${file%/}_ras.bam

	for file2 in *-1.bam
	do 
		samtools index $file2
		freebayes -f /data/hg19_neo_GFP_puro.fa -C 1 ${file2}|grep "chr11"|grep -w "534288" > ${file2%.bam}.vcf; done; for file3 in *.vcf; do a=$(awk '$4 == "C" && $5 == "A" {print $0}' ${file3}|cut -f 10|cut -d ":" -f3,5|cut -d ":" -f1); b=$(awk '$4 == "C" && $5 == "A" {print $0}' ${file3}|cut -f 10|cut -d ":" -f3,5|cut -d ":" -f2)
		echo -e ${file3%.vcf}"\t"$a"\t"$b >> /data/${file%/}SNP_Ras.txt
	done
done
