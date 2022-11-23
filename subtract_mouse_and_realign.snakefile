#subtract_mouse_and_realign.snakefile
#Anna-Lisa Doebley
#Template made 2021-09-24
#Ha Lab
#Fred Hutchinson Cancer Research Center

"""
#before running snakemake, do in tmux terminal:
ml snakemake/5.19.2-foss-2019b-Python-3.7.4
ml BWA/0.7.17-foss-2018b
ml SAMtools/1.10-GCCcore-8.3.0
ml picard/2.18.29-Java
ml Java/11.0.2
ml GATK/4.1.4.1-GCCcore-8.3.0-Java-11
ml R/3.6.2-foss-2019b-fh1


#command to run snakemake (remove -np at end when done validating):
snakemake -s subtract_mouse_and_realign.snakefile --latency-wait 60 --keep-going --cluster-config config/cluster_slurm.yaml --cluster "sbatch -p {cluster.partition} --mem={cluster.mem} -t {cluster.time} -c {cluster.ncpus} -n {cluster.ntasks} -o {cluster.output} -J {cluster.JobName}" -j 40 -np

#output file marked as temp is deleted after all rules that use it as an input are completed
"""

configfile: "config/samples.yaml"
configfile: "config/config.yaml"

rule all:
	input:
		#concatRef metrics
		expand(config['results_path']+"/{samples}/{samples}_ConcatRef_ASM.txt", samples=config["samples"]),
		expand(config['results_path']+"/{samples}/{samples}_ConcatRef_wgs_metrics.txt", samples=config["samples"]),
		
		#results from human realignment
		expand(config['results_path']+"/{samples}/{samples}_marked_dup_metrics.txt", samples=config["samples"]),
		expand(config['results_path']+"/{samples}/{samples}_recalibration_data.table", samples=config["samples"]),
		expand(config['results_path']+"/{samples}/{samples}_recalibrated.bam", samples=config["samples"]),
		expand(config['results_path']+"/{samples}/{samples}_recalibrated.bam.bai", samples=config["samples"]),
		expand(config['results_path']+"/{samples}/{samples}_ASM.txt", samples=config["samples"]),
		expand(config['results_path']+"/{samples}/{samples}_wgs_metrics.txt", samples=config["samples"]),
		expand(config['results_path']+"/{samples}/{samples}_isize.txt", samples=config["samples"]),
		expand(config['results_path']+"/{samples}/{samples}_isize.pdf", samples=config["samples"]),
		
		#results from mouse realignment


rule unmap:
	input:
		bam_file = lambda wildcards: config["samples"][wildcards.samples]
	output:
		fastq1 = temp(config['results_path']+"/concatRef_fastqs/{samples}_fastq1.fq.gz"),
		fastq2 = temp(config['results_path']+"/concatRef_fastqs/{samples}_fastq2.fq.gz"),
		unpaired_fastq = temp(config['results_path']+"/concatRef_fastqs/{samples}_unpaired.fq.gz")
	params:
		java=config["java"],
		picard_jar = config["picard_jar"],
		input_reference_genome = config["input_reference_genome"] #required for cram input
	shell:
		"{params.java} -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx85G \
		-jar {params.picard_jar} SamToFastq I={input.bam_file} \
		TMP_DIR="+config['results_path']+"/concatRef_tmps \
		REFERENCE_SEQUENCE={params.input_reference_genome} \
		FASTQ={output.fastq1} \
		SECOND_END_FASTQ={output.fastq2} \
		UNPAIRED_FASTQ={output.unpaired_fastq}"


rule map_to_ConcatRef:
	input:
		fastq1=config['results_path']+"/concatRef_fastqs/{samples}_fastq1.fq.gz",
		fastq2=config['results_path']+"/concatRef_fastqs/{samples}_fastq2.fq.gz"
	output:
		temp(config['results_path']+"/{samples}/{samples}_ConcatRef_unsorted.bam")
	params:
		reference_genome=config["ConcatRef_genome"],
		bwa=config["bwa"],
		samtools=config["samtools"],
		bwa_threads=config["bwa_threads"]
	shell:
		"({params.bwa} mem -t {params.bwa_threads} -M \
		-R '@RG\\tID:no_id\\tLB:no_library\\tPL:no_platform\\tPU:no_unit\\tSM:{wildcards.samples}' \
		{params.reference_genome} \
		{input.fastq1} {input.fastq2} | {params.samtools} view -b - > {output})"


rule sort_concatRef_by_coord:
	input:
		config['results_path']+"/{samples}/{samples}_ConcatRef_unsorted.bam"
	output:
		temp(config['results_path']+"/{samples}/{samples}_ConcatRef_sorted.bam")
	params:
		samtools=config["samtools"]
	shell:
		"({params.samtools} sort -o {output} {input})"


rule index_ConcatRef_bam:
	input:
		sorted_bam = config['results_path']+"/{samples}/{samples}_ConcatRef_sorted.bam"
	output:
		sorted_bam_index = temp(config['results_path']+"/{samples}/{samples}_ConcatRef_sorted.bam.bai")
	params:
		samtools=config["samtools"]
	shell:
		"{params.samtools} index {input.sorted_bam} {output.sorted_bam_index}"


rule get_ConcatRef_wgs_metrics:
	input:
		sorted_bam = config['results_path']+"/{samples}/{samples}_ConcatRef_sorted.bam",
		bai_index = config['results_path']+"/{samples}/{samples}_ConcatRef_sorted.bam.bai"
	output:
		config['results_path']+"/{samples}/{samples}_ConcatRef_wgs_metrics.txt"
	params:
		reference_genome=config["ConcatRef_genome"],
		java=config["java"],
		picard_jar=config["picard_jar"]
	run: 
		shell("{params.java} -jar {params.picard_jar} CollectWgsMetrics \
		I={input.sorted_bam} \
		O={output} \
		R={params.reference_genome} \
		TMP_DIR="+config['results_path']+"/concatRef_tmps \
		COUNT_UNPAIRED=true \
		USE_FAST_ALGORITHM=true \
		INCLUDE_BQ_HISTOGRAM=true")
		

rule get_ConcatRef_alignment_metrics:
	input:
		sorted_bam = config['results_path']+"/{samples}/{samples}_ConcatRef_sorted.bam",
		bai_index = config['results_path']+"/{samples}/{samples}_ConcatRef_sorted.bam.bai"
	output:
		config['results_path']+"/{samples}/{samples}_ConcatRef_ASM.txt"
	params:
		reference_genome=config["ConcatRef_genome"],
		java=config["java"],
		picard_jar=config["picard_jar"]
	shell:
		"{params.java} -jar {params.picard_jar} CollectAlignmentSummaryMetrics \
		R={params.reference_genome} \
		I={input.sorted_bam} \
		O={output} \
		TMP_DIR="+config['results_path']+"/concatRef_tmps"

rule remove_mouse:
	input:
		sorted_bam = config['results_path']+"/{samples}/{samples}_ConcatRef_sorted.bam",
		bai_index = config['results_path']+"/{samples}/{samples}_ConcatRef_sorted.bam.bai"
	output:
		sorted_bam_step_1 = temp(config['results_path']+"/{samples}/{samples}_ConcatRef_sorted.human.bam"),
		sorted_bam_step_1_index = temp(config['results_path']+"/{samples}/{samples}_ConcatRef_sorted.human.bam.bai"),
		sorted_clean_bam = temp(config['results_path']+"/{samples}/{samples}_ConcatRef_sorted.human.cleaned.bam"),
		sorted_clean_bam_index = temp(config['results_path']+"/{samples}/{samples}_ConcatRef_sorted.human.cleaned.bam.bai")
	params:
		path = config['results_path']+"/{samples}/",
		filename = "{samples}_ConcatRef_sorted.bam",
		tag = config["tag"]
	shell:
		"sh ./scripts/mod_pipe_ConcatRef.sh {params.path} {params.filename} {params.tag}"

#new rule
rule remove_human:
	input:
		sorted_bam = config['results_path']+"/{samples}/{samples}_ConcatRef_sorted.bam",
		bai_index = config['results_path']+"/{samples}/{samples}_ConcatRef_sorted.bam.bai"
	output:
		sorted_bam_step_1 = temp(config['results_path']+"/{samples}/{samples}_ConcatRef_sorted.mouse.bam"),
		sorted_bam_step_1_index = temp(config['results_path']+"/{samples}/{samples}_ConcatRef_sorted.mouse.bam.bai"),
		sorted_clean_bam = temp(config['results_path']+"/{samples}/{samples}_ConcatRef_sorted.mouse.cleaned.bam"),
		sorted_clean_bam_index = temp(config['results_path']+"/{samples}/{samples}_ConcatRef_sorted.mouse.cleaned.bam.bai")
	params:
		path = config['results_path']+"/{samples}/",
		filename = "{samples}_ConcatRef_sorted.bam",
		human_tag = config["human_tag"]
	shell:
		"sh ./scripts/mod_pipe_ConcatRef_mouse.sh {params.path} {params.filename} {params.human_tag}"


#the following rules realign post mouse subtraction
rule unmap_human_cleaned:
	input:
		bam_file = config['results_path']+"/{samples}/{samples}_ConcatRef_sorted.human.cleaned.bam"
	output:
		fastq1 = temp(config['results_path']+"/human_fastqs/{samples}_fastq1.fq.gz"),
		fastq2 = temp(config['results_path']+"/human_fastqs/{samples}_fastq2.fq.gz"),
		fastqUnpaired = temp(config['results_path']+"/human_fastqs/{samples}_fastq_unpaired.fq.gz")
	params:
		java=config["java"],
		picard_jar = config["picard_jar"]
	shell:
		"{params.java} -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx70G -jar \
		{params.picard_jar} SamToFastq \
		I={input.bam_file} \
		TMP_DIR="+config['results_path']+"/human_tmps \
		FASTQ={output.fastq1} \
		SECOND_END_FASTQ={output.fastq2} \
		UNPAIRED_FASTQ={output.fastqUnpaired}"

rule map_to_human_ref:
	input:
		fastq1=config['results_path']+"/human_fastqs/{samples}_fastq1.fq.gz",
		fastq2=config['results_path']+"/human_fastqs/{samples}_fastq2.fq.gz"
	output:
		temp(config['results_path']+"/{samples}/{samples}_human.unsorted.bam")
	params:
		reference_genome=config["human_reference_genome"],
		bwa=config["bwa"],
		samtools=config["samtools"],
		bwa_threads=config["bwa_threads"]
	shell:
		"({params.bwa} mem -t {params.bwa_threads} -M \
		-R '@RG\\tID:no_id\\tLB:no_library\\tPL:no_platform\\tPU:no_unit\\tSM:{wildcards.samples}' \
		{params.reference_genome} \
		{input.fastq1} {input.fastq2} | {params.samtools} view -b - > {output})"

#the following rules realign post HUMAN subtraction
rule unmap_mouse_cleaned:
	input:
		bam_file = config['results_path']+"/{samples}/{samples}_ConcatRef_sorted.mouse.cleaned.bam"
	output:
		fastq1 = temp(config['results_path']+"/mouse_fastqs/{samples}_fastq1.fq.gz"),
		fastq2 = temp(config['results_path']+"/mouse_fastqs/{samples}_fastq2.fq.gz"),
		fastqUnpaired = temp(config['results_path']+"/mouse_fastqs/{samples}_fastq_unpaired.fq.gz")
	params:
		java=config["java"],
		picard_jar = config["picard_jar"]
	shell:
		"{params.java} -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx70G -jar \
		{params.picard_jar} SamToFastq \
		I={input.bam_file} \
		TMP_DIR="+config['results_path']+"/human_tmps \
		FASTQ={output.fastq1} \
		SECOND_END_FASTQ={output.fastq2} \
		UNPAIRED_FASTQ={output.fastqUnpaired}"

rule map_to_mouse_ref:
	input:
		fastq1=config['results_path']+"/mouse_fastqs/{samples}_fastq1.fq.gz",
		fastq2=config['results_path']+"/mouse_fastqs/{samples}_fastq2.fq.gz"
	output:
		temp(config['results_path']+"/{samples}/{samples}_mouse.unsorted.bam")
	params:
		reference_genome=config["mouse_reference_genome"],
		bwa=config["bwa"],
		samtools=config["samtools"],
		bwa_threads=config["bwa_threads"]
	shell:
		"({params.bwa} mem -t {params.bwa_threads} -M \
		-R '@RG\\tID:no_id\\tLB:no_library\\tPL:no_platform\\tPU:no_unit\\tSM:{wildcards.samples}' \
		{params.reference_genome} \
		{input.fastq1} {input.fastq2} | {params.samtools} view -b - > {output})"

rule mark_dups_cleaned_human:
	input:
		config['results_path']+"/{samples}/{samples}_human.unsorted.bam"
	output:
		bam=temp(config['results_path']+"/{samples}/{samples}_human.dedup.bam"),
		metrics=protected(config['results_path']+"/{samples}/{samples}_human.marked_dup_metrics.txt")
	params:
		samtools=config["samtools"],
		bwa_threads=config["bwa_threads"]
	shell:
		"({params.samtools} samtools sort -n -m 2G -@ {params.bwa_threads} {input} \
    	| samtools fixmate -m -@ {params.bwa_threads} - - \
    	| samtools sort -m 2G -@ {params.bwa_threads} - \
    	| samtools markdup -r -s -f {output.metrics} --barcode-name -@ {params.bwa_threads} \
        - {output.bam})"

rule mark_dups_cleaned_mouse:
	input:
		config['results_path']+"/{samples}/{samples}_mouse.unsorted.bam"
	output:
		bam=temp(config['results_path']+"/{samples}/{samples}_mouse.dedup.bam"),
		metrics=protected(config['results_path']+"/{samples}/{samples}_mouse.marked_dup_metrics.txt")
	params:
		samtools=config["samtools"],
		bwa_threads=config["bwa_threads"]
	shell:
		"({params.samtools} samtools sort -n -m 2G -@ {params.bwa_threads} {input} \
    	| samtools fixmate -m -@ {params.bwa_threads} - - \
    	| samtools sort -m 2G -@ {params.bwa_threads} - \
    	| samtools markdup -r -s -f {output.metrics} --barcode-name -@ {params.bwa_threads} \
        - {output.bam})"

rule mark_dups_cleaned_human:
	input:
		config['results_path']+"/{samples}/{samples}_human.dedup.bam"
	output:
		config['results_path']+"/{samples}/{samples}_human.bed"
	params:
		samtools=config["samtools"]
	shell:
		"({params.samtools} samtools view -B -h {input} \
| bedtools bamtobed -bedpe -i stdin \
| cut -f1,2,6,7 | sort -k1,1 -k2n,2n -k3n,3n \
| awk -v OFS='\t' '{len = $3 - $2 ; print $0, len }' > {output})"

rule mark_dups_cleaned_mouse:
	input:
		config['results_path']+"/{samples}/{samples}_mouse.dedup.bam"
	output:
		config['results_path']+"/{samples}/{samples}_mouse.bed"
	params:
		samtools=config["samtools"]
	shell:
		"({params.samtools} samtools view -B -h {input} \
| bedtools bamtobed -bedpe -i stdin \
| cut -f1,2,6,7 | sort -k1,1 -k2n,2n -k3n,3n \
| awk -v OFS='\t' '{len = $3 - $2 ; print $0, len }' > {output})"