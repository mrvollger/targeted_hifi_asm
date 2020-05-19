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


config = {
	"HG00733" :{
		"hifi":"/net/eichler/vol27/projects/hgsvc/nobackups/data/sequencing/HiFi/HG00733/ccs/fastq.fofn",
		"pat":"/net/eichler/vol26/projects/koren_hifi_asm/nobackups/heng_hifiasm/pat_ill.fasta",
		"mat":"/net/eichler/vol26/projects/koren_hifi_asm/nobackups/heng_hifiasm/mat_ill.fasta"
		}
	}



SMS=list(config.keys())
PARENTS = ["pat", "mat"]
wildcard_constraints:
	SM= "|".join(SMS),
	PAR = "|".join(PARENTS),

workdir: "hifiasm_out"


rule all:
	input:
		hap1s = expand("{SM}.asm.hap1.p_ctg.gfa.fasta", SM=SMS),
		hap2s = expand("{SM}.asm.hap2.p_ctg.gfa.fasta", SM=SMS),


def get_illumina(wc):
	return(config[wc.SM][wc.PAR] )

rule yac:
	input:
		illumina = get_illumina, 
	output:
		yak = "temp/{SM}.{PAR}.yac"
	threads: 32
	shell:"""
{YAK} count -k31 -b37 -t{threads} -o {output.yak} {input.illumina}
"""


def get_hifi(wc):
	return(config[wc.SM]["hifi"] )

rule hifiasm:
	input:
		hifi=get_hifi,
		pat_yak = "temp/{SM}.pat.yac",
		mat_yak = "temp/{SM}.mat.yac",
	output:
		dipy = "{SM}.asm.dip.r_utg.gfa",
		dipn = "{SM}.asm.dip.r_utg.noseq.gfa",
		hap1 = "{SM}.asm.hap1.p_ctg.gfa",
		hap1n = "{SM}.asm.hap1.p_ctg.noseq.gfa",
		hap2 = "{SM}.asm.hap2.p_ctg.gfa",
		hap2n = "{SM}.asm.hap2.p_ctg.noseq.gfa",
		ec = "{SM}.asm.ec.bin",
		reverse = "{SM}.asm.ovlp.reverse.bin",
		source = "{SM}.asm.ovlp.source.bin",
	threads:60
	shell:"""
{HIFIASM} -o {wildcards.SM}.asm -t {threads} \
	-1 {input.pat_yak} -2 {input.mat_yak} $(cat {input.hifi})
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
awk '/^S/{{print ">"$2"\n"$3}}' {input.hap1} > {output.hap1} && samtools faidx {output.hap1}
awk '/^S/{{print ">"$2"\n"$3}}' {input.hap2} > {output.hap2} && samtools faidx {output.hap2}
"""




