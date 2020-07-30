import os 
import sys
from Bio import SeqIO
import itertools
import re
import tempfile
import glob
import pandas as pd

SCRIPTS = os.path.dirname(workflow.snakefile) + "/scripts"
ENV = os.path.dirname(workflow.snakefile) + "/env.cfg"
shell.prefix("source {ENV}; set -eo pipefail; ")


YAK="/net/eichler/vol26/projects/koren_hifi_asm/nobackups/heng_hifiasm/yak/yak"
HIFIASM="/net/eichler/vol26/projects/koren_hifi_asm/nobackups/heng_hifiasm/hifiasm/hifiasm"
TRF="/net/eichler/vol26/projects/sda_assemblies/nobackups/software/SDA/externalRepos/trf-v4.09/bin/trf"


# set the tmp dir
SSD_TMP_DIR = "/data/scratch/ssd"
if "TMPDIR" in os.environ:
    TMPDIR = os.environ['TMPDIR']
elif "TMPDIR" in config:
    TMPDIR = config['TMPDIR']
elif os.path.exists(SSD_TMP_DIR):
	TMPDIR = SSD_TMP_DIR
else:
    TMPDIR = tempfile.gettempdir()



configfile: "targeted_trio.yaml"

# define the samples 
SMS=list(config.keys())

#
# define the regions 
#
regions = {}
for line in open(config["regions"]):
	if(line[0] == "#"):
		continue 
	t = line.strip().split()
	if(len(t) != 4):
		key = "_".join(t)
	else:
		key = t[3]	
	regions.setdefault(key, [])
	regions[key].append( (t[0], int(t[1]), int(t[2])) )
RGNS=list(regions.keys())
SMS.remove("regions")

#
# define parameters 
#

# define ref
REF=config["ref"]; SMS.remove("ref")

# set up max threads
THREADS=16 
if("threads" in SMS): THREADS=config["threads"]; SMS.remove("threads")

# set up min asm length 
MIN_ASM_LEN=25000
if("minasmlen" in SMS): MIN_ASM_LEN=config["minasmlen"]; SMS.remove("minasmlen")

#
# pull out the trios
#
TRIOS = {}
for SM in config["trios"]:
	TRIOS[SM] = {"pat": config["trios"][SM]["pat"], "mat": config["trios"][SM]["mat"] }
SMS.remove("trios")
TRIO_SMS = list(TRIOS.keys())

#
# Get the reads
#
READ_SMS = {} # hifi inputs
ILL_SMS = {}  # illumina inputs
for SM in SMS:
	if("fofn" in config[SM] ):
		READ_SMS[SM] = [line.strip() for line in open(config[SM]["fofn"]).readlines()]
	if("illumina" in config[SM]):
		ILL_SMS[SM] = config[SM]["illumina"]

#
# function that keeps temp files if in DEBUG mode
#
DEBUG=True
def tempd(f):
	if(DEBUG): return(f)
	return(temp(f))

PARENTS = ["pat", "mat"]

wildcard_constraints:
	SM="|".join(SMS),
	RGN="|".join(RGNS),
	PAR = "|".join(PARENTS),


#
# setup working dir
#
OUTDIR="Targeted_trio"
if("outdir" in SMS): OUTDIR=config["outdir"]; SMS.remove("outdir")
workdir: OUTDIR




rule all:
	input:
		hap1s = expand("asm/{RGN}/{SM}.asm.hap1.p_ctg.gfa.fasta", SM=TRIO_SMS, RGN=RGNS),
		hap2s = expand("asm/{RGN}/{SM}.asm.hap2.p_ctg.gfa.fasta", SM=TRIO_SMS, RGN=RGNS),
		pngs  = expand("asm/{RGN}/{SM}.cov.png", SM=TRIO_SMS, RGN=RGNS),
		beds  = expand("asm/{RGN}/{SM}.to_ref.bed", SM=TRIO_SMS, RGN=RGNS),
	shell:"""
seq_stats.py -r {input.hap1s} {input.hap2s} 2> /dev/null | column -t 
"""

#
# Align the reads
# 
def get_read(wildcards):
        SM = str(wildcards.SM)
        ID = int(str(wildcards.ID))
        return(READ_SMS[SM][ID])

rule align:
        input:
                ref=REF,
                reads=get_read,
        output:
                bam=temp("temp/{SM}_{ID}.bam"),
        benchmark:
                "logs/align_{SM}_{ID}.b"
        resources:
                mem = 8,
                smem = 4
        threads: THREADS
        shell:"""
SAM_TMP="{TMPDIR}/temp_{wildcards.SM}_{wildcards.ID}"
rm -f $SAM_TMP*
pbmm2 align --unmapped --log-level DEBUG --preset SUBREAD -j {threads} {input.ref} {input.reads} | samtools view -F 256 -u - | samtools sort -T $SAM_TMP -m {resources.smem}G -@ {threads} - > {output.bam} 
"""


#
# Merge the alignments 
#
def get_bams(wildcards):
        SM = str(wildcards.SM)
        IDS = list(range(len(READ_SMS[SM])))
        bams = expand("temp/" + SM + "_{ID}.bam", ID=IDS)
        return(bams)

rule merge:
	input:
		bams = get_bams,
	output:
		bam=protected("alignments/{SM}.bam"),
		bai=protected("alignments/{SM}.bam.bai"),
	benchmark:
		"logs/merge_{SM}.b"
	resources:
		mem = 2,
	threads: min(16, THREADS)
	shell:"""
samtools merge -@ {threads} {output.bam} {input.bams} && samtools index {output.bam}
"""

rule get_coverage:
	input:
		bam = rules.merge.output.bam,
		bai = rules.merge.output.bai,
		fai = REF + ".fai",
	output:
		cov = protected("alignments/{SM}.coverage.tbl"),
	benchmark:
		"logs/coverage_{SM}.b"
	resources:
		mem = 64,
	threads: 1
	shell:"""
bedtools coverage -bed -mean -sorted \
	-g {input.fai} \
	-a <(awk '{{if($2 > 1000000){{print $1"\t"0"\t"$2}}}}' {input.fai}) \
	-b {input.bam}  \
	> {output.cov}
"""

rule make_alns:
	input:
		bams=expand(rules.merge.output.bam,SM=SMS),
		bais=expand(rules.merge.output.bai,SM=SMS),
		covs=expand(rules.get_coverage.output.cov,SM=SMS),
	output:
		done="alignments/done.txt",
	resources:
		mem = 2,
	threads: 1
	shell:"""
touch {output.done}
"""


#
# Get fastq input by region and sample 
#
def get_rgn(wildcards):
	rgns = regions[str(wildcards.RGN)] 
	rtn = ""
	for rgn in rgns: rtn += "'{}:{}-{}' ".format(*rgn)
	return(rtn) 
		

rule fastq:
	input:
		bam=ancient(rules.merge.output.bam),
		bai=ancient(rules.merge.output.bai),
	output:
		bam=tempd("asm/{RGN}/temp/{SM}.bam"),
		fastq=tempd("asm/{RGN}/temp/{SM}.fastq"),
		fai=tempd("asm/{RGN}/temp/{SM}.fastq.fai"),
	params:
		rgn = get_rgn,
	shell:"""
#samtools merge -R {params.rgn} {output.bam} {input.bam}
samtools view -b {input.bam} {params.rgn} > {output.bam}
samtools fastq {output.bam} > {output.fastq}
samtools faidx {output.fastq}
"""



def get_parent(wc):
	if(wc.SM in ILL_SMS):
		return( ILL_SMS[wc.SM] )
	else:
		return( (rules.merge.output.bam).format(SM=wc.SM) )


rule yak_reads:
	input:
		reads = get_parent, 
	output:
		reads = temp("yak/{SM}.{PAR}.reads"),
	threads: 4
	run:
		if(input.reads[-5:] == ".cram"):
			shell("samtools fasta -@ {threads} {input.reads} > {output.reads}")
		elif(input.reads[-4:] == ".bam"):
			shell("samtools fasta -@ {threads} {input.reads} > {output.reads}")
		elif(input.reads[-9:] == ".fastq.gz"):
			shell("ln -s {input.reads} {output.reads}")
		else:
			print("Bad input for yak")
	
rule yak:
	input:
		reads = rules.yak_reads.output.reads, 
	output:
		yak = protected("yak/{SM}.{PAR}.yak"),
	threads: 32
	shell:"""
{YAK} count -K 1000000000 -k31 -b37 -t{threads} -o {output.yak} {input.reads}
"""


def get_mat(wc):
	PAR="mat"
	SM = TRIOS[wc.SM][PAR] 
	yak = (rules.yak.output.yak).format(SM=SM, PAR=PAR, RGN=wc.RGN)
	return(yak)

def get_pat(wc):
	PAR="pat"
	SM = TRIOS[wc.SM][PAR] 
	yak = (rules.yak.output.yak).format(SM=SM, PAR=PAR, RGN=wc.RGN)
	return(yak)

rule hifiasm:
	input:
		hifi = rules.fastq.output.fastq,
		mat_yak = get_mat,
		pat_yak = get_pat,
	output:
		dipy = "asm/{RGN}/{SM}.asm.dip.r_utg.gfa",
		dipn = "asm/{RGN}/{SM}.asm.dip.r_utg.noseq.gfa",
		hap1 = "asm/{RGN}/{SM}.asm.hap1.p_ctg.gfa",
		hap1n = "asm/{RGN}/{SM}.asm.hap1.p_ctg.noseq.gfa",
		hap2 = "asm/{RGN}/{SM}.asm.hap2.p_ctg.gfa",
		hap2n = "asm/{RGN}/{SM}.asm.hap2.p_ctg.noseq.gfa",
		ec = "asm/{RGN}/{SM}.asm.ec.bin",
		reverse = "asm/{RGN}/{SM}.asm.ovlp.reverse.bin",
		source = "asm/{RGN}/{SM}.asm.ovlp.source.bin",
	threads: THREADS
	shell:"""
{HIFIASM} -o asm/{wildcards.RGN}/{wildcards.SM}.asm -t {threads} \
	-1 {input.mat_yak} -2 {input.pat_yak} {input.hifi}
"""


rule make_fasta:	
	input:
		hap1 = rules.hifiasm.output.hap1,
		hap2 = rules.hifiasm.output.hap2,
	output:
		hap1 = rules.hifiasm.output.hap1 + ".fasta",
		hap2 = rules.hifiasm.output.hap2 + ".fasta",
		hap1fai = rules.hifiasm.output.hap1 + ".fasta.fai",
		hap2fai = rules.hifiasm.output.hap2 + ".fasta.fai",
	threads:1 
	shell:"""
awk '/^S/{{print ">"$2"_{wildcards.SM}_mat\\n"$3}}' {input.hap1} > {output.hap1} && samtools faidx {output.hap1}
awk '/^S/{{print ">"$2"_{wildcards.SM}_pat\\n"$3}}' {input.hap2} > {output.hap2} && samtools faidx {output.hap2}
"""


rule nucplot:
	input:
		hap1 = rules.make_fasta.output.hap1,
		hap2 = rules.make_fasta.output.hap2,
		fastq= rules.fastq.output.fastq,
	output:
		fasta = temp("asm/{RGN}/{SM}.tmp.fasta"),
		bam = "asm/{RGN}/{SM}.bam",
		bai = "asm/{RGN}/{SM}.bam.bai",
		png = "asm/{RGN}/{SM}.cov.png",
	threads: THREADS
	shell:"""
cat {input.hap1} {input.hap2} > {output.fasta}
pbmm2 align --log-level DEBUG --preset SUBREAD --min-length 5000 -j {threads} \
	{output.fasta} {input.fastq} | samtools view -F 2308 -u - | samtools sort - > {output.bam}
samtools index {output.bam}
{SCRIPTS}/NucPlot.py --width 32 --height 5 --soft -c 1000 {output.bam} {output.png}
"""


rule map_to_ref:
	input:
		ref = REF,
		fasta = rules.nucplot.output.fasta,
	output:
		bam = "asm/{RGN}/{SM}.to_ref.bam",
		bed = "asm/{RGN}/{SM}.to_ref.bed",
	threads:THREADS
	shell:"""
minimap2 --eqx -Y -ax asm20 -t {threads} -r 50000 \
	{input.ref} {input.fasta} | samtools view -F 260 -u - | samtools sort - > {output.bam}

samtools index {output.bam}
~/projects/hifi_asm/scripts/samIdentity.py --bed {output.bam} > {output.bed}
"""



