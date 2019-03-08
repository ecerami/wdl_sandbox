task convertSam2Bam {
	File sam_input

	command {
		samtools view -b -S -o sim_reads_aligned.bam ${sam_input}
	}

	output {
		File aligned_bam = "sim_reads_aligned.bam"
	}

	runtime {
		docker: "mgibio/samtools:1.9" 
	}
}

task sortBam {
	File bam_input

	command {
		samtools sort -o sim_reads_aligned.sorted.bam ${bam_input}
	}

	output {
		File aligned_sorted_bam = "sim_reads_aligned.sorted.bam"
	}

	runtime {
		docker: "mgibio/samtools:1.9" 
	}
}

task indexBam {
	File aligned_sorted_bam_input

	command {
		samtools index ${aligned_sorted_bam_input} sim_reads_aligned.sorted.bam.bai
	}
	
	output {
		File aligned_sorted_bam_index = "sim_reads_aligned.sorted.bam.bai"
	}	

	runtime {
		docker: "mgibio/samtools:1.9" 
	}
}

workflow ecoliWorkflow {
	Array[File] sam_files

	scatter(sam_file in sam_files) {
		call convertSam2Bam { 
			input: sam_input=sam_file
		}
		call sortBam {
			input: bam_input=convertSam2Bam.aligned_bam
		}
		call indexBam {
			input: aligned_sorted_bam_input=sortBam.aligned_sorted_bam
		}
	}
}
