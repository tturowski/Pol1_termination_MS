import re, os, subprocess
import pandas as pd

# command used to create processing env
# mamba create -n processing -y -c conda-forge -c bioconda deeptools pyCRAC salmon star subread fastx_toolkit samtools
# run
# snakemake -c12 -j12 --use-conda -slurm -s SM_CRACprocessing3end_all.smk


#paths
path = "00_raw/"
name_elem = '.fastq.gz'  #write file ending here

#references
STAR_INDEX = "seq_references/EF4.74_STAR_index/"
GTF = "seq_references/Saccharomyces_cerevisiae.EF4.74.shortChNames_with_PolIII_transcripts_extended_slop_intergenic_sort.gtf"

#parsing file names and preparatory jobs
# longName = [n.strip(name_elem) for n in os.listdir(path) if n.endswith(name_elem) and 'Rpa190HTP_wt_none_6' in n]
longName = [n.replace(name_elem,"") for n in os.listdir(path) if n.endswith(name_elem)]
barcodes = [n.split("_")[0] for n in longName]
SAMPLES = ["_".join(n.split("_")[1:]) for n in longName]

def parseBarcode(bc):
    toGrep = "".join(re.findall(r'([ATCG])', bc))
    [firstN,postN]=[len(i) for i in bc.split(toGrep)]
    return firstN, toGrep, postN

def readLengths(n):
	f = path+n+name_elem
	command = "zcat "+f+" | head -40 | scripts/fastqReadsLength.awk | cut -f1"
	return int(subprocess.run([command],shell=True, stdout=subprocess.PIPE).stdout.decode('utf-8'))

df_names = pd.DataFrame({
	'barcode' : barcodes,
	'bcLen'	  : [len(bc)+1 for bc in barcodes], #+1 for fastx_trimmer
	'longName' : longName,
	'readLength': [readLengths(n) for n in longName],
	'name'	  : SAMPLES}).set_index("name")

print(df_names)
df_names.to_csv('names.tab', sep="\t")

d1_name = df_names['longName'].to_dict()
d2_bcLen = df_names['bcLen'].to_dict()
d4_readLen = df_names['readLength'].to_dict()

#SnakeMake pipeline

########## OUTPUTS ##########

rule all:
	input:
		expand("01_preprocessing/01e_{sample}_3end.fasta.gz",sample=SAMPLES),
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

########## PREPROCESSING ##########

def bcFile(wildcards):
	sample_name = wildcards.sample
	return path+d1_name[sample_name]+name_elem

# droped QC to avoid keeping PCR duplicates artificially
# rule QC:
# 	input:
# 		bcFile
# 	params:
# 		"01_preprocessing/01a_{sample}_QC"
# 	output:
# 		"01_preprocessing/01a_{sample}_QC.fastq.gz"
# 	conda:
# 		"envs/flexbar.yml"
# 	shell:
# 		"flexbar -r {input} -t {params} -q TAIL -qf i1.8 -qt 20 -z GZ"

rule collapsing:
	input:
		bcFile
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

rule flexbar_3end_trimming:
	input:
		"01_preprocessing/01c_{sample}_deBC.fasta.gz"
	params:
		adSeq = "TGGAATTCTCGGGTGCCAAGGC",
		out_prefix = "01_preprocessing/01d_{sample}_flexbar"
	output:
		"01_preprocessing/01d_{sample}_flexbar.fasta.gz"
	conda:
		"envs/flexbar.yml" 
	shell:
		"flexbar -r {input} -t {params.out_prefix} -as {params.adSeq} -ao 4 -u 3 -m 7 -n 4 -bt RIGHT -z GZ"

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

rule align:
	input:
		reads = "01_preprocessing/01e_{sample}_3end.fasta.gz",
	params:
		index_dir = STAR_INDEX,
		prefix = "02_alignment/{sample}_STAR_"
	output:
		bam = "02_alignment/{sample}_STAR_Aligned.out.bam"
	conda:
		"envs/processing.yml"
	shell:
		"STAR --outFileNamePrefix {params.prefix} --readFilesCommand zcat --genomeDir {params.index_dir} --genomeLoad LoadAndKeep --outSAMtype BAM Unsorted --readFilesIn {input.reads} --limitOutSJcollapsed 2000000"


ruleorder: align > sort
########## POSTPROCESSING ##########

rule sort:
	input:
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

# rule clean:
# 	input:
# 		expand("02_alignment/{sample}.bam.bai",sample=SAMPLES)
# 	output:
# 		touch("cleaninig.done")
# 	conda:
# 		"envs/processing.yml"
# 	shell:
# 		"""
# 		rm -r logs/
# 		rm Aligned.out.sam Log.final.out Log.out Log.progress.out SJ.out.tab
# 		"""

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