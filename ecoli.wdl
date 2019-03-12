task align2Genome {
	File simulated_read_file
	String genomeDir
	
	command {
		bowtie2 -x ${genomeDir}/indexes/e_coli -U ${simulated_read_file} -S sim_reads_aligned.sam
	}

	output {
		File sam_file = "sim_reads_aligned.sam"
	}

	runtime {
		docker: 'biocontainers/bowtie2:v2.2.9_cv2' 
	}
}

task convertSam2Bam {
	File sam_file

	command {
		samtools view -b -S -o sim_reads_aligned.bam ${sam_file}
	}

	output {
		File aligned_bam_file = "sim_reads_aligned.bam"
	}

	runtime {
		docker: 'drjimbo/cidc-bwa-samtools:latest' 
	}
}

task sortBam {
	File bam_file

	command {
		samtools sort -o sim_reads_aligned.sorted.bam ${bam_file}
	}

	output {
		File aligned_sorted_bam_file = "sim_reads_aligned.sorted.bam"
	}

	runtime {
		docker: "drjimbo/cidc-bwa-samtools:latest" 
	}
}

task indexBam {
	File aligned_sorted_bam_file

	command {
		samtools index ${aligned_sorted_bam_file} sim_reads_aligned.sorted.bam.bai
	}
	
	output {
		File aligned_sorted_bam_index_file = "sim_reads_aligned.sorted.bam.bai"
	}	

	runtime {
		docker: "drjimbo/cidc-bwa-samtools:latest" 
	}
}

task callVariants1 {
	File aligned_sorted_bam_file
	String genomeDir
	
	command {
		samtools mpileup -g -f ${genomeDir}/genomes/NC_008253.fna ${aligned_sorted_bam_file} > variants.bcf
	}

	output {
		File bcf_variants_file = "variants.bcf"
	}

	runtime {
		docker: "drjimbo/cidc-bwa-samtools:latest" 
	}
}

task callVariants2 {
	File bcf_variants_file
	
	command {
		bcftools call -c -v ${bcf_variants_file} > variants.vcf
	}

	output {
		File sam = "variants.vcf"
	}

	runtime {
		docker: "biocontainers/bowtie2:v2.2.9_cv2" 
	}
}

workflow ecoliWorkflow {
	Array[File] simulated_read_files
	String genomeDir

	scatter(simulated_read_file in simulated_read_files) {
		call align2Genome {
			input: simulated_read_file=simulated_read_file, genomeDir=genomeDir
		}
		call convertSam2Bam { 
			input: sam_file=align2Genome.sam_file
		}
		call sortBam {
			input: bam_file=convertSam2Bam.aligned_bam_file
		}
		call indexBam {
			input: aligned_sorted_bam_file=sortBam.aligned_sorted_bam_file
		}
		call callVariants1 {
			input: aligned_sorted_bam_file=sortBam.aligned_sorted_bam_file, genomeDir=genomeDir
		}
		call callVariants2 {
			input: bcf_variants_file=callVariants1.bcf_variants_file
		}
	}
}
