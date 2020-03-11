import os 
import sys
from Bio import SeqIO
import itertools
import re
import tempfile
import glob

shell.prefix("set -eo pipefail; ")

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
regions = { "_".join(line.strip().split()) : "{}:{}-{}".format(*line.strip() .split()) for line in open(config["regions"]) }
RGNS = list(regions.keys())
SMS.remove("regions")

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
                bam=temp("Targeted_HiFi_Asm/temp/{SM}_{ID}.bam"),
        benchmark:
                "Targeted_HiFi_Asm/logs/align_{SM}_{ID}.b"
        resources:
                mem = 8,
                smem = 4
        threads: THREADS
        shell:"""
SAM_TMP="{TMPDIR}/temp_{wildcards.SM}_{wildcards.ID}"
rm -f $SAM_TMP*
pbmm2 align --log-level DEBUG --preset SUBREAD -j {threads} {input.ref} {input.reads} | samtools view -F 260 -u - | samtools sort -T $SAM_TMP -m {resources.smem}G -@ {threads} - > {output.bam} 
"""


#
# Merge the alignments 
#
def get_bams(wildcards):
        SM = str(wildcards.SM)
        IDS = list(range(len(READ_SMS[SM])))
        bams = expand("Targeted_HiFi_Asm/temp/" + SM + "_{ID}.bam", ID=IDS)
        return(bams)

rule merge:
	input:
		bams = get_bams,
	output:
		bam="Targeted_HiFi_Asm/alignments/{SM}.bam",
		bai="Targeted_HiFi_Asm/alignments/{SM}.bam.bai",
	benchmark:
		"Targeted_HiFi_Asm/logs/merge_{SM}.b"
	threads: min(16, THREADS)
	shell:"""
samtools merge -@ {threads} {output.bam} {input.bams} && samtools index {output.bam}
"""


#
# Get fastq input by region and sample 
#
def get_rgn(wildcards):
	return( regions[str(wildcards.RGN)] )

rule fastq:
	input:
		bam=rules.merge.output.bam,
		bai=rules.merge.output.bai,
	output:
		bam=temp("Targeted_HiFi_Asm/asm/{RGN}/{SM}.bam"),
		fastq=temp("Targeted_HiFi_Asm/asm/{RGN}/{SM}.fastq"),
	params:
		rgn = get_rgn,
	shell:"""
samtools merge -R '{params.rgn}' {output.bam} {input.bam}
samtools fastq {output.bam} > {output.fastq}
"""


#
# Run the local assembly 
# 
def get_region_size(wildcards):
	RGN = str(wildcards.RGN)	
	tokens = RGN.split("_")
	length = int(tokens[-1]) - int(tokens[-2])
	return(int(length))

def get_min_ovl(wildcards):
	return(str(wildcards.MIN_OVL).lstrip("0") )

rule canu_asm:
	input:
		rules.fastq.output.fastq,
	output:
		canu_dir = temp(directory("Targeted_HiFi_Asm/asm/{RGN}/canu_{SM}_{MIN_OVL}/")),
		fasta = temp("Targeted_HiFi_Asm/asm/{RGN}/temp_{SM}_{MIN_OVL}.contigs.fasta"),
		readNames = "Targeted_HiFi_Asm/asm/{RGN}/read_info/{SM}_{MIN_OVL}.readNames.txt",
		readToTig = "Targeted_HiFi_Asm/asm/{RGN}/read_info/{SM}_{MIN_OVL}.layout.readToTig",
	threads: 4
	params:
		gsize = get_region_size,
		min_ovl = get_min_ovl, 
	shell:'''
canu \
	-d {output.canu_dir} \
    -p {wildcards.SM} \
    useGrid=false \
	maxThreads={threads} \
	genomeSize={params.gsize} \
	batOptions="-eg 0.0 -eM 0.0 -mo {params.min_ovl} -dg 3 -db 3 -dr 1 -ca 50 -cp 5"  \
    -pacbio-hifi {input} && \
	mv {output.canu_dir}/{wildcards.SM}.contigs.fasta {output.fasta} && \
    mv {output.canu_dir}/{wildcards.SM}.contigs.layout.readToTig {output.readToTig} && \
    mv {output.canu_dir}/{wildcards.SM}.seqStore/readNames.txt {output.readNames}
'''

rule rename_asm:
	input:
		fastq =	rules.fastq.output.fastq,
		fasta = rules.canu_asm.output.fasta,
		canu_dir = rules.canu_asm.output.canu_dir, 
		readNames = rules.canu_asm.output.readNames,
		readToTig = rules.canu_asm.output.readToTig,
	output:
		fasta = "Targeted_HiFi_Asm/asm/{RGN}/{SM}.{MIN_OVL}.contigs.fasta",
		fai =   "Targeted_HiFi_Asm/asm/{RGN}/{SM}.{MIN_OVL}.contigs.fasta.fai",
		outdir = directory("Targeted_HiFi_Asm/asm/{RGN}/{SM}.{MIN_OVL}/"),
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
		fasta = rules.rename_asm.output.fasta,
		fai   = rules.rename_asm.output.fai,
		outdir = rules.rename_asm.output.outdir, 
	output:
		png = "Targeted_HiFi_Asm/asm/{RGN}/{SM}.{MIN_OVL}.plots.png",
		bam = temp("Targeted_HiFi_Asm/asm/{RGN}/{SM}.{MIN_OVL}.plots.bam"),
		bai = temp("Targeted_HiFi_Asm/asm/{RGN}/{SM}.{MIN_OVL}.plots.bam.bai"),
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
			shell(f"samtools index {bam} && NucPlot.py {bam} {png}")
	
		shell("samtools merge {output.bam} {input.outdir}/*.bam && samtools index {output.bam}")
		shell("NucPlot.py --soft -c 100 {output.bam} {output.png}")	

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






