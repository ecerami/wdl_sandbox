echo "Aligning to e-coli genome"
bowtie2 -x indexes/e_coli -U simulated_reads/sim_reads_1.fq -S alignments/sim_reads_aligned_1.sam

echo "Converting SAM to BAM"
samtools view -b -S -o alignments/sim_reads_aligned_1.bam alignments/sim_reads_aligned_1.sam

echo "Sorting BAM"
samtools sort -o alignments/sim_reads_aligned_1.sorted.bam alignments/sim_reads_aligned_1.bam

echo "Indexing BAM"
samtools index alignments/sim_reads_aligned_1.sorted.bam

echo "Calling variants"
samtools mpileup -g -f genomes/NC_008253.fna alignments/sim_reads_aligned_1.sorted.bam > alignments/sim_variants_1.bcf

bcftools call -c -v alignments/sim_variants_1.bcf > alignments/sim_variants_1.vcf

echo "Checking files created"
ls -l alignments