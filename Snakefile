configfile:
	"config_files/config.yaml"

sample_ids, = glob_wildcards(config['reads']+"/{sample}.R1.fastq.gz")
REF, = glob_wildcards(config['ref']+"/{reference}.fa")
outputdir = config['output']

print(sample_ids,)
print(REF,)


rule all:
	input:
		expand("{output}/snippyout/{reference}/{sample}.out", sample=sample_ids, reference=REF, output = config['output']),		
		expand("{output}/fasttree/{reference}.clean.fullcore.tree", reference=REF, output = config['output']),
		expand("{output}/snp_dists/{reference}.pairwise_snps.csv", reference=REF, output = config['output'])		

rule snippy_run:
	input:
		ref = config['ref']+"/{reference}.fa",
		r1 = config['reads']+"/{sample}.R1.fastq.gz",
		r2 = config['reads']+"/{sample}.R2.fastq.gz"
	output: 
		directory("{output}/snippyout/{reference}/{sample}.out")
	conda: 
		"config_files/snippy.yaml"
	priority:
		10
	shell: 
		"snippy --outdir {output} --ref {input.ref} --R1 {input.r1} --R2 {input.r2}"


rule snippy_core:
	input:
		ref = config['ref']+"/{reference}.fa",	
	output:	
		"{output}/core/{reference}.full.aln"
	conda:
		"config_files/snippy.yaml"
	priority:
		9
	shell:
		"""
		snippy-core --prefix {wildcards.output}/core/{wildcards.reference} --ref {input.ref} {wildcards.output}/snippyout/{wildcards.reference}/*.out
		"""
	
rule snippy_clean:
	input:
		"{output}/core/{reference}.full.aln"
	output:
		"{output}/core/{reference}.clean.full.aln"
	conda: 
		"config_files/snippy.yaml"	
	shell:
		"snippy-clean_full_aln {input} > {output}"

rule gubbins:
	input:
		"{output}/core/{reference}.clean.full.aln"
	output:
		"{output}/gubbins/{reference}.filtered_polymorphic_sites.fasta"
	params:
		prefix = "{output}/gubbins/{reference}",
		filt = config["gubbins"]["params"]
	conda:
		"config_files/gubbins.yaml"
	resources:
		mem_mb=config["gubbins"]["memory"]
	shell:
		"""
		run_gubbins.py -u -v {params.filt} -p {params.prefix} {input}
		rm {wildcards.reference}.clean.full.aln.seq.joint.txt
		"""
rule snp_sites:
	input:
		"{output}/gubbins/{reference}.filtered_polymorphic_sites.fasta"
	output:
		"{output}/snp_sites/{reference}.clean.fullcore.aln"
	conda: 
		"config_files/snippy.yaml"
	shell:
		"snp-sites -c {input} > {output}"

rule snp_dists:
	input:
		"{output}/snp_sites/{reference}.clean.fullcore.aln"
	output:
		"{output}/snp_dists/{reference}.pairwise_snps.csv"
	conda: 
		"config_files/snp_dists.yaml"
	shell: 
		"snp-dists -c {input} > {output}"

rule fasttree:
	input:
		"{output}/snp_sites/{reference}.clean.fullcore.aln"
	output:
		"{output}/fasttree/{reference}.clean.fullcore.tree"
	conda:
		"config_files/fasttree.yaml"
	shell:
		"fasttree -gtr -nt {input} > {output}"






