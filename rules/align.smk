def get_fq(wildcards):
	if config["trimming"]["skip"]:
		# no trimming, use raw reads
		return units.loc[(wildcards.sample, wildcards.unit), ["fq1", "fq2"]].dropna()
	else:
		# yes trimming, use trimmed data
		if not is_single_end(**wildcards):
			# paired-end sample
			return expand("results/trimmed/{sample}-{unit}.{group}.fastq.gz",
				group=[1, 2], **wildcards)
        # single end sample
		return "results/trimmed/{sample}-{unit}.fastq.gz".format(**wildcards)

rule star_align:
	input:
		sample=get_fq
	output:
        # see STAR manual for additional output files
		"results/star/{sample}-{unit}Aligned.out.bam"
	log:
		"logs/star/{sample}-{unit}.log"
	params:
		# path to STAR reference genome index
		star_dir="results/star/{sample}-{unit}",
		index=config["ref"]["index"],
        # optional parameters
		extra=config["params"]["star"]
	threads: 24
    #wrapper:
    #    "0.19.4/bio/star/align"
	
	shell:
		"unset TMPDIR; module load star/2.7.5c; set -euo pipefail;  " +
		"STAR {params.extra} --runThreadN {threads} --genomeDir {params.index} " +
		"--readFilesIn {input.sample} --outSAMtype BAM Unsorted " +
		"--outFileNamePrefix {params.star_dir} --outStd Log  > {log} 2>&1"
