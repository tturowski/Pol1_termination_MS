import re, os, subprocess
import pandas as pd

# command used to create processing env
# mamba create -n processing -y -c conda-forge -c bioconda deeptools pyCRAC salmon star subread fastx_toolkit samtools

#paths
path = "00_raw/"
name_elem = '.fastq.gz'  #write file ending here

#references
STAR_INDEX = "/home/tomasz.turowski/seq_references/mm10/mm10_STAR_index/"
GTF = "/home/tomasz.turowski/seq_references/mm10/ENCFF871VGR.gtf"

#parsing file names and preparatory jobs
longName = [n.strip(name_elem) for n in os.listdir(path) if n.endswith(name_elem)]
adaptors = [n.split("_")[0] for n in longName]
barcodes = [n.split("_")[1] for n in longName]
SAMPLES = ["_".join(n.split("_")[2:]) for n in longName]

def parseBarcode(bc):
    toGrep = "".join(re.findall(r'([ATCG])', bc))
    [firstN,postN]=[len(i) for i in bc.split(toGrep)]
    return firstN, toGrep, postN

def readLengths(n):
	f = path+n+name_elem
	command = "zcat "+f+" | head -40 | fastqReadsLength.awk | cut -f1"
	return int(subprocess.run([command],shell=True, stdout=subprocess.PIPE).stdout.decode('utf-8'))

df_names = pd.DataFrame({
	'adaptor' : adaptors,
	'barcode' : barcodes,
	'bcLen'	  : [len(bc)+1 for bc in barcodes], #+1 for fastx_trimmer
	'longName' : longName,
	'readLength': [readLengths(n) for n in longName],
	'name'	  : SAMPLES}).set_index("name")

print(df_names)
df_names.to_csv('names.tab', sep="\t")

d1_name = df_names['longName'].to_dict()
d2_bcLen = df_names['bcLen'].to_dict()
d3_as = df_names['adaptor'].to_dict()
d4_readLen = df_names['readLength'].to_dict()

#SnakeMake pipeline

########## OUTPUTS ##########

rule all:
	input:
		expand("01_preprocessing/01e_{sample}_3end.fasta.gz",sample=SAMPLES),
		expand("02_alignment/{sample}.bam",sample=SAMPLES),
		expand("02_alignment/{sample}.bam.bai",sample=SAMPLES),
		"cleaninig.done",
		"03_FetaureCounts/featureCounts_multimappers.list",
		"03_FetaureCounts/featureCounts_uniq.list",
		expand("04_BigWig/{sample}_raw_plus.bw",sample=SAMPLES),
		expand("04_BigWig/{sample}_raw_minus.bw",sample=SAMPLES),
		expand("04_BigWig/{sample}_CPM_plus.bw",sample=SAMPLES),
		expand("04_BigWig/{sample}_CPM_minus.bw",sample=SAMPLES)

########## PREPROCESSING ##########

def bcFile(wildcards):
	sample_name = wildcards.sample
	return path+d1_name[sample_name]+name_elem

rule QC:
	input:
		bcFile
	params:
		"01_preprocessing/01a_{sample}_QC"
	output:
		"01_preprocessing/01a_{sample}_QC.fastq.gz"
	conda:
		"envs/flexbar.yml"
	shell:
		"flexbar -r {input} -t {params} -q TAIL -qf i1.8 -qt 20 -z GZ"

rule collapsing:
	input:
		"01_preprocessing/01a_{sample}_QC.fastq.gz"
	output:
		"01_preprocessing/01b_{sample}_comp.fasta.gz"
	conda:
		"envs/processing.yml"
	shell:
		"gunzip -c {input} | fastx_collapser | gzip > {output}"

#functions for debarcoding
def bcLen(wildcards):
	sample_name = wildcards.sample
	return d2_bcLen[sample_name]

rule debarcoding:
	input:
		"01_preprocessing/01b_{sample}_comp.fasta.gz"
	params:
		bcLen = bcLen
	output:
		"01_preprocessing/01c_{sample}_deBC.fasta.gz"
	conda:
		"envs/processing.yml"
	shell:
		"gunzip -c {input} | fastx_trimmer -f {params.bcLen} | gzip > {output}"

#functions for adaptor sequence
def adSeq(wildcards):
	sample_name = wildcards.sample
	return d3_as[sample_name]

rule flexbar_3end_trimming:
	input:
		"01_preprocessing/01c_{sample}_deBC.fasta.gz"
	params:
		adSeq = adSeq
	output:
		"01_preprocessing/01d_{sample}_flexbar.fasta.gz"
	conda:
		"envs/flexbar.yml" 
	shell:
		"flexbar -r {input} -t {params} -as params.adSeq -ao 4 -u 3 -m 7 -n 4 -bt RIGHT -z GZ"

def maxLen(wildcards):
	sample_name = wildcards.sample
	readLen = d4_readLen[sample_name]
	bcLen = d2_bcLen[sample_name]
	return readLen-bcLen-4

rule length_filtering:
	input:
		"01_preprocessing/01d_{sample}_flexbar.fasta.gz"
	params:
		maxLen = maxLen
	output:
		"01_preprocessing/01e_{sample}_3end.fasta.gz"
	shell:
		"zcat {input} | ./scripts/lenFilterFastaMax.awk -v var={params.maxLen} | gzip > {output}"

########## ALIGNMENT ##########

rule load_genome:
	input:
		index_check = STAR_INDEX+"SAindex",
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
		reads = "01_preprocessing/01e_{sample}_3end.fasta.gz",
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

rule BigWigs_CPM:
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

rule BigWigs_raw:
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
