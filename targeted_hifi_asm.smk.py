import os 
import sys
from Bio import SeqIO
import itertools
import re
shell.prefix("set -eo pipefail; ")


configfile: "targeted_hifi_asm.yaml"
PLOT_SCRIPT="/net/eichler/vol26/home/mvollger/projects/nucfreq/NucPlot.py"
TRF="/net/eichler/vol26/projects/sda_assemblies/nobackups/software/SDA/externalRepos/trf-v4.09/bin/trf"
MIN_ASM_LEN=10000


regions = { "_".join(line.strip().split()) : "{}:{}-{}".format(*line.strip() .split()) for line in open(config["regions"]) }
RGNS = list(regions.keys())
print(regions)

SMS=list(config["samples"].keys())
print(SMS)

MIN_OVLS = [ "{:05d}".format(x) for x in range(500,5250,250) ]

wildcard_constraints:
	SM="|".join(SMS),
	RGN="|".join(RGNS),
	MIN_OVL="\d+",
	SM1="|".join(SMS),
	SM2="|".join(SMS),
	SM3="|".join(SMS),

rule all:
	input:
		"Targeted_HiFi_Asm/assembly.stats.txt",


#
# get fastq input by region and sample 
#
def get_bam(wildcards):
	return( config["samples"][str(wildcards.SM)] )

def get_region_size(wildcards):
	RGN = str(wildcards.RGN)	
	tokens = RGN.split("_")
	#print(tokens)
	length = int(tokens[-1]) - int(tokens[-2])
	return max(int(length/3), 1000)

def get_rgn(wildcards):
	return( regions[str(wildcards.RGN)] )

rule fastq:
	input:
		get_bam,
	output:
		bam=temp("Targeted_HiFi_Asm/asm/{RGN}/{SM}.bam"),
		fastq=temp("Targeted_HiFi_Asm/asm/{RGN}/{SM}.fastq"),
	params:
		rgn = get_rgn,
	shell:"""
samtools merge -R {params.rgn} {output.bam} {input}
samtools fastq {output.bam} > {output.fastq}
"""

def get_min_ovl(wildcards):
	return(str(wildcards.MIN_OVL).lstrip("0") )
rule canu_asm:
	input:
		rules.fastq.output.fastq,
	output:
		canu_dir = temp(directory("Targeted_HiFi_Asm/asm/{RGN}/canu_{SM}_{MIN_OVL}/")),
		fasta = temp("Targeted_HiFi_Asm/asm/{RGN}/temp_{SM}_{MIN_OVL}.contigs.fasta"),
	threads: 4
	params:
		gsize = get_region_size,
		min_ovl = get_min_ovl, 
	shell:"""
module load canu/1.9
canu \
	-d {output.canu_dir} \
    -p {wildcards.SM} \
    useGrid=false \
	maxThreads={threads} \
	genomeSize={params.gsize} \
	batOptions="-eg 0.0 -eM 0.0 -mo {params.min_ovl} -dg 3 -db 3 -dr 1 -ca 50 -cp 5"  \
    -pacbio-hifi {input} && \
	mv {output.canu_dir}/{wildcards.SM}.contigs.fasta {output.fasta}
"""
#contigFilter="2 {MIN_ASM_LEN} 1.0 0.5 4" \

rule rename_asm:
	input:
		fasta = rules.canu_asm.output.fasta,
		canu_dir = rules.canu_asm.output.canu_dir, 
	output:
		fasta = "Targeted_HiFi_Asm/asm/{RGN}/{SM}.{MIN_OVL}.contigs.fasta",
		fai =   "Targeted_HiFi_Asm/asm/{RGN}/{SM}.{MIN_OVL}.contigs.fasta.fai",
	run:
		recs = list(SeqIO.parse(input["fasta"], "fasta"))
		# order by size 
		recs.sort(key = (lambda rec: len(rec.seq) ), reverse=True)

		out = []
		for rec in recs:
			if( len(recs) == 1 or MIN_ASM_LEN <= len(rec.seq)):
				rec.id = wildcards["SM"] + "__" + wildcards["RGN"] + "__" + rec.id
				out.append(rec)
		SeqIO.write(out, output["fasta"], "fasta")
		shell("samtools faidx {output.fasta}")


rule plots:
	input:
		fastq =	rules.fastq.output.fastq,
		fasta = rules.rename_asm.output.fasta,
		fai   = rules.rename_asm.output.fai,
	output:
		png = "Targeted_HiFi_Asm/asm/{RGN}/{SM}.{MIN_OVL}.plots.png",
		bam = temp("Targeted_HiFi_Asm/asm/{RGN}/{SM}.{MIN_OVL}.bam"),
		pai = temp("Targeted_HiFi_Asm/asm/{RGN}/{SM}.{MIN_OVL}.bam.bai"),
	threads: 4
	shell:"""
pbmm2 align --log-level DEBUG --preset SUBREAD --min-length 5000 -j {threads} \
	{input.fasta} {input.fastq} | \
	samtools view -F 2308 -u - | samtools sort - > {output.bam} && \
	samtools index {output.bam} && \
	{PLOT_SCRIPT} {output.bam} {output.png}
"""

rule stats:	
	input:
		fai = expand(rules.rename_asm.output.fai, SM=SMS, RGN=RGNS, MIN_OVL=MIN_OVLS),
		plots = expand(rules.plots.output.png, SM=SMS, RGN=RGNS, MIN_OVL=MIN_OVLS),
	output:
		"Targeted_HiFi_Asm/assembly.stats.txt",
	shell:"""
head -n 100 {input.fai} > {output}
"""



def get_asms(wildcards):
	SM1 = str(wildcards.SM1)
	SM2 = str(wildcards.SM2)
	SM3 = str(wildcards.SM3)
	RGN = str(wildcards.RGN)
	return(sorted(expand(rules.rename_asm.output.fasta, RGN=RGN, SM = [SM1, SM2, SM3]) ))

rule dot_plot:
	input:
		get_asms,
	output:
		pdf = "Targeted_HiFi_Asm/dot_plot/{RGN}/{SM1}__{SM2}__{SM3}.pdf",
		fasta = temp("Targeted_HiFi_Asm/dot_plot/{RGN}/{SM1}__{SM2}__{SM3}.fasta"),
	shell:"""
cat {input}  > {output.fasta}
python /net/eichler/vol27/projects/ruiyang_projects/nobackups/vntr_project/dotplots/dotplots/DotplotMain.py \
   	{output.fasta} -o {output.pdf} -w 100
"""


rule dot_plots:
	input:
		pdfs = expand( "Targeted_HiFi_Asm/dot_plot/{RGN}/chm1__chm13__{SM}.pdf", RGN=RGNS, SM=SMS) 
		#pdfs = [ "Targeted_HiFi_Asm/dot_plot/{SM1}__{SM2}__{SM3}.pdf".format(SM1=a,SM2=b,SM3=c) for a,b,c in itertools.combinations(sorted(SMS), 3) ]
	output:
		"Targeted_HiFi_Asm/dot_plot/done.txt"
	shell:"""
touch {output}
"""


rule trf:
	input:
		fasta = rules.rename_asm.output.fasta,
	output:
		trf = "Targeted_HiFi_Asm/asm/{RGN}/{SM}.trf.txt",
		tmptrf = temp("Targeted_HiFi_Asm/asm/{RGN}/tmp.{SM}.trf.txt"),
	run:
		shell("{TRF} {input.fasta} 2 7 7 80 10 50 2000 -ngs > {output.tmptrf}")
		out = "\t".join('contig start end PeriodSize CopyNumber ConsensusSize PercentMatches PercentIndels Score A C G T Entropy Motif Sequence L_Flank R_Flank'.split() ) + "\n"
		contig = ""
		for line in open(output["tmptrf"]).readlines():
			tokens = line.strip().split()
			if( tokens[0][0] == "@"):
				contig = tokens[0][1:]	
			else:
				out += "\t".join( [contig] + tokens ) + "\n"
		open(output["trf"], "w+").write(out)	



rule merge_trf:
	input:
		expand(rules.trf.output.trf, SM=SMS, RGN=RGNS)




rule run_msa:
	input:
	output:
	shell:"""
/net/eichler/vol26/projects/sda_assemblies/nobackups/software/prank/bin/prank
"""






