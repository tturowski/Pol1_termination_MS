import re, os, subprocess

# command used to create processing env
# mamba create -n processing -y -c conda-forge -c bioconda deeptools pyCRAC salmon star subread fastx_toolkit samtools

#paths
path = "bifx_processing/"
name_elem = '.fasta.gz'  #write file ending here

#references
STAR_INDEX = "seq_references/EF4.74_STAR_index/"
GTF = "seq_references/Saccharomyces_cerevisiae.EF4.74.shortChNames_with_PolIII_transcripts_extended_slop_intergenic_sort.gtf"

#parsing file names and preparatory jobs
SAMPLES = [n.replace(name_elem,"") for n in os.listdir(path) if n.endswith(name_elem)]
print(SAMPLES)

#SnakeMake pipeline

########## OUTPUTS ##########

rule all:
	input:
		STAR_INDEX+"index.done",
		expand("02_alignment/{sample}.bam",sample=SAMPLES),
		expand("02_alignment/{sample}.bam.bai",sample=SAMPLES),
		"cleaninig.done",
		"03_FetaureCounts/featureCounts_multimappers.list",
		"03_FetaureCounts/featureCounts_uniq.list",
		expand("04_BigWig/{sample}_raw_plus.bw",sample=SAMPLES),
		expand("04_BigWig/{sample}_raw_minus.bw",sample=SAMPLES),
		expand("04_BigWig/{sample}_CPM_plus.bw",sample=SAMPLES),
		expand("04_BigWig/{sample}_CPM_minus.bw",sample=SAMPLES),
		expand("04_BigWig/{sample}.sam",sample=SAMPLES)
		# expand("04_BigWig/{sample}_PROFILE_read_fwd.bw",sample=SAMPLES),
		# expand("04_BigWig/{sample}_PROFILE_read_rev.bw",sample=SAMPLES),
		# expand("04_BigWig/{sample}_PROFILE_5end_fwd.bw",sample=SAMPLES),
		# expand("04_BigWig/{sample}_PROFILE_5end_rev.bw",sample=SAMPLES),
		# expand("04_BigWig/{sample}_PROFILE_3end_fwd.bw",sample=SAMPLES),
		# expand("04_BigWig/{sample}_PROFILE_3end_rev.bw",sample=SAMPLES),
		# expand("04_BigWig/{sample}_PROFILE_3end_polyA_fwd.bw",sample=SAMPLES),
		# expand("04_BigWig/{sample}_PROFILE_3end_polyA_rev.bw",sample=SAMPLES)	

########## ALIGNMENT ##########

rule genome_generate:
	input:
		fasta_file = "seq_references/Saccharomyces_cerevisiae.EF4.74.dna.toplevel.shortChrNames.fasta",
		gtf_file = GTF
	output:
		touch(STAR_INDEX+"index.done"),
		index_check = STAR_INDEX+"SAindex",
		outdir = directory(STAR_INDEX)
	conda:
		"envs/processing.yml"
	shell:
		"STAR --runThreadN 10 --genomeSAindexNbases 10 --runMode genomeGenerate --genomeDir {output.outdir} --genomeFastaFiles {input.fasta_file} --sjdbGTFfile {input.gtf_file}"

rule load_genome:
	input:
		index_check = STAR_INDEX+"SAindex",
		check_file = STAR_INDEX+"index.done"
	params:
		index_dir = STAR_INDEX
	output:
		touch('loading.done')
	conda:
		"envs/processing.yml"
	shell:
		"STAR --genomeLoad LoadAndExit --genomeDir {params.index_dir}"

rule align:
	input:
		reads = "bifx_processing/{sample}.fasta.gz",
		idx = 'loading.done'
	params:
		index_dir = STAR_INDEX,
		prefix = "02_alignment/{sample}_STAR_"
	output:
		bam = "02_alignment/{sample}_STAR_Aligned.out.bam"
	conda:
		"envs/processing.yml"
	shell:
		"STAR --outFileNamePrefix {params.prefix} --readFilesCommand zcat --genomeDir {params.index_dir} --genomeLoad LoadAndKeep --outSAMtype BAM Unsorted --readFilesIn {input.reads}"

rule unload_genome:
	# Delete the loading.done flag file otherwise subsequent runs of the pipeline 
	# will fail to load the genome again if STAR alignment is needed.
	input:
		bam = expand("02_alignment/{sample}_STAR_Aligned.out.bam",sample=SAMPLES),
		idx = 'loading.done'
	params:
		index_dir = STAR_INDEX
	output:
		"logs/STARunload_Log.out",
	conda:
		"envs/processing.yml"
	shell:
		"""
		STAR --genomeLoad Remove --genomeDir {params.index_dir} --outFileNamePrefix logs/STARunload_
		rm {input.idx}
		"""
ruleorder: align > sort
########## POSTPROCESSING ##########

rule sort:
	input:
		log = "logs/STARunload_Log.out", #wait to unload_genome
		bam = "02_alignment/{sample}_STAR_Aligned.out.bam"
	output:
		bam = "02_alignment/{sample}.bam",
		bai = "02_alignment/{sample}.bam.bai"
	conda:
		"envs/processing.yml"
	shell:
		"""
		samtools sort {input.bam} > {output.bam}
		samtools index {output.bam}
		"""

rule clean:
	input:
		expand("02_alignment/{sample}.bam.bai",sample=SAMPLES)
	output:
		touch("cleaninig.done")
	conda:
		"envs/processing.yml"
	shell:
		"""
		rm -r logs/
		rm Aligned.out.sam Log.final.out Log.out Log.progress.out SJ.out.tab
		"""

rule featureCounts:
	input:
		bam = expand("02_alignment/{sample}.bam",sample=SAMPLES) #use list of files
	output:
		multi = "03_FetaureCounts/featureCounts_multimappers.list",
		uniq = "03_FetaureCounts/featureCounts_uniq.list"
	params:
		gtf=GTF
	conda:
		"envs/processing.yml"
	shell:
		"""
		featureCounts -M -s 1 -a {params.gtf} -o {output.multi} {input.bam}
		featureCounts -s 1 -a {params.gtf} -o {output.uniq} {input.bam}
		"""

rule bamcoverage_CPM:
	input:
		bam = "02_alignment/{sample}.bam",
		bai = "02_alignment/{sample}.bam.bai"
	output:
		bwP = "04_BigWig/{sample}_CPM_plus.bw",
		bwM = "04_BigWig/{sample}_CPM_minus.bw"
	conda:
		"envs/processing.yml"
	shell:
		"""
		bamCoverage --bam {input.bam} -of bigwig -o {output.bwP} --filterRNAstrand reverse --normalizeUsing CPM --binSize 1
		bamCoverage --bam {input.bam} -of bigwig -o {output.bwM} --filterRNAstrand forward --normalizeUsing CPM --binSize 1
		"""

rule bamcoverage_raw:
	input:
		bam = "02_alignment/{sample}.bam",
		bai = "02_alignment/{sample}.bam.bai"
	output:
		bwP = "04_BigWig/{sample}_raw_plus.bw",
		bwM = "04_BigWig/{sample}_raw_minus.bw"
	conda:
		"envs/processing.yml"
	shell:
		"""
		bamCoverage --bam {input.bam} -of bigwig -o {output.bwP} --filterRNAstrand reverse --binSize 1
		bamCoverage --bam {input.bam} -of bigwig -o {output.bwM} --filterRNAstrand forward --binSize 1
		"""

rule bam2sam:
	input:
		bam = "02_alignment/{sample}.bam",
		bai = "02_alignment/{sample}.bam.bai"
	output:
		sam = "04_BigWig/{sample}.sam"
	conda:
		"envs/processing.yml"
	shell:
		"""
		samtools view -h {input.bam} > {output.sam}
		"""

# rule trxtools_5end:
# 	input:
# 		"04_BigWig/{sample}.sam"
# 	output:
# 		"04_BigWig/{sample}_PROFILE_5end_fwd.bw",
# 		"04_BigWig/{sample}_PROFILE_5end_rev.bw"
# 	conda:
# 		"envs/processing.yml"
# 	shell:
# 		"""
# 		SAM2profilesGenomic.py -u 5end -f {input}
# 		"""

# rule trxtools_3end:
# 	input:
# 		"04_BigWig/{sample}.sam"
# 	output:
# 		"04_BigWig/{sample}_PROFILE_3end_fwd.bw",
# 		"04_BigWig/{sample}_PROFILE_3end_rev.bw",
# 		"04_BigWig/{sample}_PROFILE_3end_polyA_fwd.bw",
# 		"04_BigWig/{sample}_PROFILE_3end_polyA_rev.bw"
# 	conda:
# 		"envs/processing.yml"
# 	shell:
# 		"""
# 		SAM2profilesGenomic.py -u 3end -n -f {input}
# 		"""

# rule trxtools_reads:
# 	input:
# 		"04_BigWig/{sample}.sam"
# 	output:
# 		"04_BigWig/{sample}_PROFILE_read_fwd.bw",
# 		"04_BigWig/{sample}_PROFILE_read_rev.bw"
# 	conda:
# 		"envs/processing.yml"
# 	shell:
# 		"""
# 		SAM2profilesGenomic.py -u read -f {input}
# 		"""