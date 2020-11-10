import os 
import sys
from Bio import SeqIO
import itertools
import re
import tempfile
import glob
import pandas as pd

SCRIPTS = os.path.dirname(workflow.snakefile) + "/scripts"
BIN = os.path.dirname(workflow.snakefile) + "/bin"
ENV = os.path.dirname(workflow.snakefile) + "/env.cfg"
shell.prefix("source {ENV}; set -eo pipefail; ")
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



configfile: "TAT.yaml"

# define ref
REF=config["ref"]

if("gff" in config): GFF  = config["gff"]

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

# set up max threads
THREADS=16 
if("threads" in config): THREADS=config["threads"]

# set up min asm length 
MIN_ASM_LEN=25000
if("minasmlen" in config): MIN_ASM_LEN=config["minasmlen"]

# set up minimum overlaps 
MIN_OVLS=[2000]
if("minovls" in config): MIN_OVLS = config["minovls"]


#
# Get the reads
#
df = pd.read_csv(config["samples"], sep="\t")
READ_SMS = {}
TRIOS = {}
ILL_SMS = {}  # illumina inputs
SMS = []
for idx,row in df.iterrows():
    SM = row["sample"]
    fofn = row["fofn"]
    READ_SMS[SM] = [line.strip() for line in open(fofn).readlines()]
    if(row.isnull().values.any()):
        SMS.append(SM)
    else:
        TRIOS[SM] = {"pat": row["pat"], "mat": row["mat"] }
TRIO_SMS = list(TRIOS.keys())
PARENTS = ["pat", "mat"]

print(TRIO_SMS,len(TRIO_SMS))
print(SMS, len(SMS))

#
# setup working dir
#
OUTDIR="TAT_Asm"
if("outdir" in SMS): OUTDIR=config["outdir"]; SMS.remove("outdir")
workdir: OUTDIR


# function that keeps temp files if in DEBUG mode
DEBUG=True
def tempd(f):
	if(DEBUG): return(f)
	return(temp(f))

wildcard_constraints:
	SM="|".join(SMS+TRIO_SMS),
	RGN="|".join(RGNS),
	MIN_OVL="\d+",
	SM1="|".join(SMS),
	SM2="|".join(SMS),
	SM3="|".join(SMS),
	PAR = "|".join(PARENTS),

localrules: get_coverage, yak, yak_reads, 

rule all:
    input:
        stat = "assembly.stats.txt",
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
        ref=REF,
        bams = get_bams,
    output:
        cram=protected("alignments/{SM}.cram"),
        crai=protected("alignments/{SM}.cram.crai"),
    benchmark:
        "logs/merge_{SM}.b"
    resources:
        mem = 2,
    threads: min(16, THREADS)
    shell:"""
samtools merge -@ {threads} --write-index --reference {input.ref} -O CRAM -f {output.cram} {input.bams}
"""

rule get_coverage:
	input:
		cram = rules.merge.output.cram,
		crai = rules.merge.output.crai,
		fai = REF + ".fai",
		ref = REF,
	output:
		cov = "alignments/{SM}.mosdepth.summary.txt",
		dist = temp("alignments/{SM}.mosdepth.global.dist.txt"),
	benchmark:
		"logs/coverage_{SM}.b"
	resources:
		mem = 8,
	threads: 4
	shell:"""
{BIN}/mosdepth -f {input.ref} -t {threads} -x -n alignments/{wildcards.SM} {input.cram}
"""

rule make_alns:
	input:
		crams=expand(rules.merge.output.cram,SM=SMS+TRIO_SMS),
		crais=expand(rules.merge.output.crai,SM=SMS+TRIO_SMS),
		covs=expand(rules.get_coverage.output.cov,SM=SMS+TRIO_SMS),
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
		cram=ancient(rules.merge.output.cram),
		crai=ancient(rules.merge.output.crai),
	output:
		bam=tempd("asm/{RGN}/temp/{SM}.bam"),
		fastq=tempd("asm/{RGN}/temp/{SM}.fastq"),
		fai=tempd("asm/{RGN}/temp/{SM}.fastq.fai"),
	params:
		rgn = get_rgn,
	shell:"""
samtools view -b {input.cram} {params.rgn} > {output.bam}
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
		gsize = temp( "asm/{RGN}/temp_{SM}_region_size.txt" ),
	threads: 4
	run:
		fastq = pd.read_csv(input["fai"], sep="\t", names=["read", "length", "x", "y", "z", "w"])
		cov =	pd.read_csv(input["cov"], sep="\t")
		total = fastq.length.sum()
		median = cov["mean"].median()
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
		gsize = rules.region_size.output.gsize,
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
	genomeSize=$(cat {input.gsize}) \
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


rule stats:	
    input:
        fasta = expand(rules.pick_best_asm.output.fasta, SM=SMS, RGN=RGNS),
        fai = expand(rules.pick_best_asm.output.fai, SM=SMS, RGN=RGNS),
        png = expand(rules.pick_best_asm.output.png, SM=SMS, RGN=RGNS),
    output:
        "assembly.stats.txt",
    shell:"""
head -n 100 {input.fai} > {output}
"""
#pdf = expand(rules.best_dot_plot.output.pdf, SM=SMS, RGN=RGNS)



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




#########################################
#########################################
######### TRIO ASSEMBLY #################
#########################################
#########################################

def get_parent(wc):
    return(TRIOS[wc.SM][wc.PAR])

rule yak_reads:
    input:
        reads = get_parent, 
    output:
        reads = temp("yak/{SM}.{PAR}.reads.gz"),
    threads: 16
    run:
        if(input.reads[-5:] == ".cram"):
            shell("samtools fasta -@ {threads} {input.reads} | pigz -p {threads} > {output.reads}")
        elif(input.reads[-4:] == ".bam"):
            shell("samtools fasta -@ {threads} {input.reads} | pigz -p {threads}  > {output.reads}")
        elif(input.reads[-9:] == ".fastq.gz"):
            shell("ln -s {input.reads} {output.reads}")
        else:
            print("Bad input for yak")

rule yak:
	input:
		reads = rules.yak_reads.output.reads, 
	output:
		yak = protected("yak/{SM}.{PAR}.yak"),
	threads: 16
	shell:"""
yak count -K 1000000000 -k31 -b37 -t{threads} -o {output.yak} {input.reads}
"""

rule make_yaks:
    input:
        yaks = expand("yak/{SM}.{PAR}.yak", SM=TRIO_SMS, PAR=PARENTS),


def get_mat(wc):
	PAR="mat"
	yak = (rules.yak.output.yak).format(SM=wc.SM, PAR=PAR)
	return(yak)

def get_pat(wc):
	PAR="pat"
	yak = (rules.yak.output.yak).format(SM=wc.SM, PAR=PAR)
	return(yak)

rule hifiasm:
	input:
		hifi = rules.fastq.output.fastq,
		mat_yak = get_mat,
		pat_yak = get_pat,
	output:
		dipy = temp("asm/{RGN}/{SM}.asm.dip.r_utg.gfa"),
		dipn = temp("asm/{RGN}/{SM}.asm.dip.r_utg.noseq.gfa"),
		hap1 = temp("asm/{RGN}/{SM}.asm.hap1.p_ctg.gfa"),
		hap1n = temp("asm/{RGN}/{SM}.asm.hap1.p_ctg.noseq.gfa"),
		hap2 = temp("asm/{RGN}/{SM}.asm.hap2.p_ctg.gfa"),
		hap2n = temp("asm/{RGN}/{SM}.asm.hap2.p_ctg.noseq.gfa"),
		ec = temp("asm/{RGN}/{SM}.asm.ec.bin"),
		hreverse = temp("asm/{RGN}/{SM}.asm.ovlp.reverse.bin"),
		source = temp("asm/{RGN}/{SM}.asm.ovlp.source.bin"),
	threads: THREADS
	shell:"""
hifiasm -o asm/{wildcards.RGN}/{wildcards.SM}.asm -t {threads} \
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


def get_asms(wc):
    if(wc.SM in TRIO_SMS):
        return()
    else:
        return()


rule map_to_ref:
	input:
		ref = REF,
		fasta = rules.nucplot.output.fasta,
	output:
		bam = "asm/{RGN}/{SM}.to_ref.bam",
		bed = "asm/{RGN}/{SM}.to_ref.bed",
	threads:THREADS
	shell:"""
minimap2 -R '@RG\tID:{wildcards.SM}_ID\tSM:{wildcards.SM}' --eqx -Y -ax asm20 -t {threads} -r 50000 \
	{input.ref} {input.fasta} | samtools view -F 260 -u - | samtools sort - > {output.bam}

samtools index {output.bam}
~mvollger/projects/hifi_asm/scripts/samIdentity.py --bed {output.bam} > {output.bed}
"""



#################################################################3
# LIFTOFF
#################################################################3


rule liftoff:
    input:
        r = REF,
        gff = GFF,
        t = rules.nucplot.output.fasta,
    output:
        gff = "genes/{SM}.gff3",
        unmapped = "genes/{SM}.unmapped.gff3",
        temp = directory('genes/temp.{SM}')
    params:
        rgn=get_rgn
    threads: 32
    resources:
        mem=8,
    shell:"""
~mvollger/.local/bin/liftoff -dir {output.temp} -sc 0.95 -copies -p {threads} -r {input.r} -t {input.t} -g {input.gff} -o {output.gff} -u {output.unmapped}
"""


rule orf_gff:
    input:
        gff = rules.liftoff.output.gff,
        fasta =  rules.nucplot.output.fasta,
    output:
        gff="genes/{SM}.orf_only.gff3",
    threads: 1
    resources:
        mem=8,
    shell:"""
gffread --adj-stop -C -F -g {input.fasta} {input.gff} \
        | gffread --keep-comments -F -J -g {input.fasta} /dev/stdin \
        | gffread --keep-comments -F -M -K /dev/stdin \
        > {output.gff}
"""

rule orf_bb:
    input:
        gff = rules.orf_gff.output.gff,
        fai = rules.nucplot.output.fasta+".fai",
    output:
        bb="genes/{SM}.orf_only.bb",
        bed="genes/{SM}.orf_only.bed",
    threads: 1
    resources:
        mem=8,
    shell:"""
{SDIR}/scripts/AddUniqueGeneIDs.py {input.gff} | \
        gff3ToGenePred -geneNameAttr=gene_name -useName /dev/stdin /dev/stdout | \
        genePredToBigGenePred /dev/stdin /dev/stdout | \
        awk -F $'\t' '{{ t = $4; $4 = $13; $13 = t; print; }}' OFS=$'\t' | \
        bedtools sort -i - > {output.bed} 

bedToBigBed -extraIndex=name,name2 -type=bed12+7 -tab -as={SDIR}/templates/bigGenePred.as {output.bed} {input.fai} {output.bb}
"""

rule gff_index:
	input:
		rules.liftoff.output.gff
	output:
		gff = "genes/{SM}.gff3.gz",
		tbi = "genes/{SM}.gff3.gz.tbi",
	threads: 1
	resources:
		mem=8,
	shell:'''
bedtools sort -i {input} | bgzip > {output.gff}
tabix {output.gff}
'''

rule genes:
	input:
		gffs = expand(rules.gff_index.output.gff, SM=[SM]),
		orf = expand(rules.orf_gff.output, SM=[SM]),
		bbs = expand(rules.orf_bb.output, SM=[SM]),




