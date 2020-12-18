#Daniel Ence
#oct. 8, 2020

rule samtools_index:
	input:
		"results/merged_lane_bams/{sample}.merged.bam"
	output:
		"results/merged_lane_bams/{sample}.merged.bam.bai"
	log:
		"logs/samtools_index_sorted.{sample}.log"
	shell:
		"unset TMPDIR; module load samtools; samtools index {input}"
