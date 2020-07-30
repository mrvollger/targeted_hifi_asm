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



configfile: "targeted_hifi_asm.yaml"
TRF="/net/eichler/vol26/projects/sda_assemblies/nobackups/software/SDA/externalRepos/trf-v4.09/bin/trf"


# define the samples 
SMS=list(config.keys())

# define ref
REF=config["ref"]; SMS.remove("ref")

# define the regions 
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
	



RGNS = list(regions.keys())
SMS.remove("regions")

#
# setup working dir
#
OUTDIR="Targeted_HiFi_Asm"
if("outdir" in SMS): OUTDIR=config["outdir"]; SMS.remove("outdir")
workdir: OUTDIR


# set up max threads
THREADS=16 
if("threads" in SMS): THREADS=config["threads"]; SMS.remove("threads")

# set up min asm length 
MIN_ASM_LEN=25000
if("minasmlen" in SMS): MIN_ASM_LEN=config["minasmlen"]; SMS.remove("minasmlen")

# set up minimum overlaps 
MIN_OVLS=[2000]
if("minovls" in SMS): MIN_OVLS = config["minovls"]; SMS.remove("minovls")

# set up sample to read dictionaries 
REF_SMS = []
FOFN_SMS = []
READ_SMS = {}
CMD_SMS = []
for SM in SMS:
	if("ref" in config[SM]):
		REF_SMS.append(SM)
	if("fofn" in config[SM] ):
		FOFN_SMS.append(SM)
		READ_SMS[SM] = [line.strip() for line in open(config[SM]["fofn"]).readlines()]

# function that keeps temp files if in DEBUG mode
DEBUG=True
def tempd(f):
	if(DEBUG): return(f)
	return(temp(f))

wildcard_constraints:
	SM="|".join(SMS),
	RGN="|".join(RGNS),
	MIN_OVL="\d+",
	SM1="|".join(SMS),
	SM2="|".join(SMS),
	SM3="|".join(SMS),

rule all:
	input:
		"assembly.stats.txt",

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


#
# Run the local assembly 
# 
rule region_size:
	input:
		fastq = rules.fastq.output.fastq,
		fai = rules.fastq.output.fai,
		cov = rules.get_coverage.output.cov,
	output:
		size = temp( "asm/{RGN}/temp_{SM}_region_size.txt" ),
	threads: 4
	run:
		fastq = pd.read_csv(input["fai"], sep="\t", names=["read", "length", "x", "y", "z", "w"])
		cov =	pd.read_csv(input["cov"], sep="\t", names=["contig", "start", "end", "coverage"] )
		total = fastq.length.sum()
		median = cov.coverage.median()
		size = int(1.5*total/median)
		print(total/10**6, median, size)
		out = open(output["size"],"w+")
		out.write( "{}\n".format(size) )
		out.close()


def get_min_ovl(wildcards):
	return(str(wildcards.MIN_OVL).lstrip("0") )

rule canu_asm:
	input:
		fastq = rules.fastq.output.fastq,
		fai = rules.fastq.output.fai,
		size = rules.region_size.output.size,
	output:
		canu_dir = temp(directory("asm/{RGN}/canu_{SM}_{MIN_OVL}/")),
		fasta = temp("asm/{RGN}/temp_{SM}_{MIN_OVL}.contigs.fasta"),
		readNames = "asm/{RGN}/read_info/{SM}_{MIN_OVL}.readNames.txt",
		readToTig = "asm/{RGN}/read_info/{SM}_{MIN_OVL}.layout.readToTig",
	threads: 4
	params:
		min_ovl = get_min_ovl, 
	shell:'''
canu \
	-d {output.canu_dir} \
    -p {wildcards.SM} \
    useGrid=false \
	maxThreads={threads} \
	genomeSize=$(cat {input.size}) \
	batOptions="-eg 0.0 -eM 0.0 -mo {params.min_ovl} -dg 3 -db 3 -dr 1 -ca 50 -cp 5" \
    -pacbio-hifi {input.fastq} && \
	mv {output.canu_dir}/{wildcards.SM}.contigs.fasta {output.fasta} && \
    mv {output.canu_dir}/{wildcards.SM}.contigs.layout.readToTig {output.readToTig} && \
    mv {output.canu_dir}/{wildcards.SM}.seqStore/readNames.txt {output.readNames}

#batOptions="-eg 0.0 -eM 0.0 -mo {params.min_ovl} -dg 3 -db 3 -dr 1 -ca 50 -cp 5"
#batOptions="-eg 0.0003 -sb 0.01 -dg 0 -db 3 -dr 0 -ca 50 -cp 5 -mo {params.min_ovl}" 
'''

rule rename_asm:
	input:
		fastq =	rules.fastq.output.fastq,
		fasta = rules.canu_asm.output.fasta,
		canu_dir = rules.canu_asm.output.canu_dir, 
		readNames = rules.canu_asm.output.readNames,
		readToTig = rules.canu_asm.output.readToTig,
	output:
		fasta = "asm/{RGN}/temp/{SM}.{MIN_OVL}.contigs.fasta",
		fai =   "asm/{RGN}/temp/{SM}.{MIN_OVL}.contigs.fasta.fai",
		outdir = directory("asm/{RGN}/temp/{SM}.{MIN_OVL}/"),
	run:
		# get reads
		recs = SeqIO.to_dict(SeqIO.parse(input["fastq"], "fastq"))
		# get read names by idx
		readNames = { line.strip().split()[0] : line.strip().split()[1] for line in open(input["readNames"] ) }
		print()
		# get reads for each tig
		tigsToReads = {}
		for line in open(input["readToTig"]):
			if(line[0] == "#"): continue 
			readid, tigid = line.strip().split()[0:2]
			readname = readNames[readid]
			tigname = "tig{:08d}".format(int(tigid))
			if(tigname not in tigsToReads):
				tigsToReads[tigname] = []
			# add records to the tig
			tigsToReads[tigname].append(recs[readname])	

		recs = list(SeqIO.parse(input["fasta"], "fasta"))
		# order by size 
		recs.sort(key = (lambda rec: len(rec.seq) ), reverse=True)

		out = []
		for rec in recs:
			if( len(recs) == 1 or MIN_ASM_LEN <= len(rec.seq)):
				tmpid = str(rec.id)
				# rename contig 
				rec.id = wildcards["SM"] + "__" + wildcards["RGN"] + "__" + tmpid
				# write contig assembly to sub direcotry 
				SeqIO.write( [rec], "{}/{}.fasta".format(output["outdir"], rec.id), "fasta")
				# save reads to sub direcotry 
				SeqIO.write( tigsToReads[tmpid], "{}/{}.reads.fastq".format(output["outdir"], rec.id), "fastq")
				# add contig to final fasta for the region sample 
				out.append(rec)

		SeqIO.write(out, output["fasta"], "fasta")
		shell("samtools faidx {output.fasta}")


rule plots:
	input:
		fastq =	rules.fastq.output.fastq,
		fasta = rules.rename_asm.output.fasta,
		fai   = rules.rename_asm.output.fai,
		outdir = rules.rename_asm.output.outdir, 
	output:
		png = "asm/{RGN}/temp/{SM}.{MIN_OVL}.plots.png",
		bam = temp("asm/{RGN}/temp/{SM}.{MIN_OVL}.plots.bam"),
		bai = temp("asm/{RGN}/temp/{SM}.{MIN_OVL}.plots.bam.bai"),
	threads: int(THREADS/4)
	run:
		reads = glob.glob(input["outdir"] + "/*.reads.fastq")
		contigs = glob.glob(input["outdir"] + "/*.fasta")
		for contig, read in zip( sorted(contigs), sorted(reads) ):
			bam = contig + ".bam"
			png = contig + ".png"
			print(contig, read, bam)

			shell(f"""pbmm2 align --log-level DEBUG --preset SUBREAD --min-length 5000 -j {threads} \
				{contig} {read} | samtools view -F 2308 -u - | samtools sort - > {bam} """)
			shell(f"samtools index {bam} && {SCRIPTS}/NucPlot.py --height 5 --soft -c 100 {bam} {png}")
	
		#shell("samtools merge {output.bam} {input.outdir}/*.bam && samtools index {output.bam}")
		shell("""pbmm2 align --log-level DEBUG --preset SUBREAD --min-length 5000 -j {threads} \
			{input.fasta} {input.fastq} | samtools view -F 2308 -u - | samtools sort - > {output.bam}""")
		shell("samtools index {output.bam} && {SCRIPTS}/NucPlot.py --height 5 --soft -c 100 {output.bam} {output.png}")


rule pick_best_asm:
	input:
		fasta = expand("asm/{{RGN}}/temp/{{SM}}.{MIN_OVL}.contigs.fasta", MIN_OVL=MIN_OVLS),
		fai = expand("asm/{{RGN}}/temp/{{SM}}.{MIN_OVL}.contigs.fasta.fai", MIN_OVL=MIN_OVLS),
		png = expand("asm/{{RGN}}/temp/{{SM}}.{MIN_OVL}.plots.png", MIN_OVL=MIN_OVLS),
	output:
		fasta = "asm/{RGN}/{SM}.best.contigs.fasta",
		fai = "asm/{RGN}/{SM}.best.contigs.fasta.fai",
		png = "asm/{RGN}/{SM}.best.plots.png",
	threads: 1
	run:
		best=None
		bestscore = 0
		for fasta, fai, png in zip(input["fasta"], input["fai"], input["png"]):
			df = pd.read_csv(fai, header=None, sep="\t")
			score = sum(df[1])/len(df[1])		
			if(score > bestscore):
				best = (fasta, fai, png)
				bestscore=score
			print(score, fasta)
		print(best, fasta)
		shell(f"ln {best[0]} {{output.fasta}}")
		shell(f"ln {best[1]} {{output.fai}}")
		shell(f"ln {best[2]} {{output.png}}")

rule best_dot_plot:
	input:
		fasta = rules.pick_best_asm.output.fasta,
	output:
		pdf = "asm/{RGN}/{SM}.pdf",
	shell:"""
python /net/eichler/vol27/projects/ruiyang_projects/nobackups/vntr_project/dotplots/dotplots/DotplotMain.py \
   	{input.fasta} -o {output.pdf} -w 100
"""



rule stats:	
	input:
		fasta = expand(rules.pick_best_asm.output.fasta, SM=SMS, RGN=RGNS),
		fai = expand(rules.pick_best_asm.output.fai, SM=SMS, RGN=RGNS),
		png = expand(rules.pick_best_asm.output.png, SM=SMS, RGN=RGNS),
		pdf = expand(rules.best_dot_plot.output.pdf, SM=SMS, RGN=RGNS)
	output:
		"assembly.stats.txt",
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
		pdf = "dot_plot/{RGN}/{SM1}__{SM2}__{SM3}.pdf",
		fasta = temp("dot_plot/{RGN}/{SM1}__{SM2}__{SM3}.fasta"),
	shell:"""
cat {input}  > {output.fasta}
python /net/eichler/vol27/projects/ruiyang_projects/nobackups/vntr_project/dotplots/dotplots/DotplotMain.py \
   	{output.fasta} -o {output.pdf} -w 100
"""


rule dot_plots:
	input:
		pdfs = expand( "dot_plot/{RGN}/chm1__chm13__{SM}.pdf", RGN=RGNS, SM=SMS) 
		#pdfs = [ "dot_plot/{SM1}__{SM2}__{SM3}.pdf".format(SM1=a,SM2=b,SM3=c) for a,b,c in itertools.combinations(sorted(SMS), 3) ]
	output:
		"dot_plot/done.txt"
	shell:"""
touch {output}
"""


rule trf:
	input:
		fasta = rules.rename_asm.output.fasta,
	output:
		trf = "asm/{RGN}/{SM}.trf.txt",
		tmptrf = temp("asm/{RGN}/tmp.{SM}.trf.txt"),
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






