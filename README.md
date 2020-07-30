# targeted_hifi_asm


## Example config file
Must be named `targeted_hifi_asm.yaml`
```
ref: /net/eichler/vol26/projects/sda_assemblies/nobackups/assemblies/hg38/ucsc.hg38.no_alts.fasta
regions: regions.bed
threads: 48
minovls:
    - 500
    - 1000
    - 2000
chm13:
    fofn: /net/eichler/vol27/projects/sequence_data/nobackups/human/CHM13/PacBioHiFi/20_kbp_insert_hifi_beta/fastq.fofn
chm1:
    fofn: /net/eichler/vol27/projects/sequence_data/nobackups/human/CHM1/PacBioHiFi/chm1_10kbp_hifi.fofn
```

## Example region file
Must have contig, start, and end. If you include a 4th row it uses it as a region name. If multiple rows have the same region name reads from all the regions are grouped together and assembled. 

```
chr6    160494344       160686570       LPA
chr7    100945872       100972010
chr1    146081967       146328717       NOTCH2
chr1    148532846       148779515       NOTCH2
chr1    120705588       120836006       NOTCH2
chr1    149366551       149480964       NOTCH2
chr1    119989248       120163868       NOTCH2
```
ref: /net/eichler/vol26/projects/sda_assemblies/nobackups/assemblies/hg38/ucsc.hg38.no_alts.fasta
regions: regions.50kbp.slop.bed
threads: 8


## Example config for targeted trio assembly ##
Runs basically the same way but with the `hifiasm_trio.smk` snakemake.

Must have a file called `targeted_trio.yaml`. 

```
ref: /net/eichler/vol26/projects/sda_assemblies/nobackups/assemblies/hg38/ucsc.hg38.no_alts.fasta
regions: regions.bed
threads: 8


#
# HGSVC
#
HG00512:
    illumina: /net/eichler/vol26/projects/1kg_high_coverage/nobackups/data/Illumina/cram/HG00512.final.cram
HG00513:
    illumina: /net/eichler/vol26/projects/1kg_high_coverage/nobackups/data/Illumina/cram/HG00513.final.cram
HG00514:
    fofn: /net/eichler/vol27/projects/hgsvc/nobackups/data/sequencing/HiFi/HG00514/ccs/fastq.fofn

HG00731:
    illumina: /net/eichler/vol26/projects/1kg_high_coverage/nobackups/data/Illumina/cram/HG00731.final.cram
HG00732:
    illumina: /net/eichler/vol26/projects/1kg_high_coverage/nobackups/data/Illumina/cram/HG00732.final.cram
HG00733:
    fofn: /net/eichler/vol27/projects/hgsvc/nobackups/data/sequencing/HiFi/HG00733/ccs/fastq.fofn

NA19238:
    illumina: /net/eichler/vol26/projects/1kg_high_coverage/nobackups/data/Illumina/cram/NA19238.final.cram
NA19239:
    illumina: /net/eichler/vol26/projects/1kg_high_coverage/nobackups/data/Illumina/cram/NA19239.final.cram
NA19240:
    fofn: /net/eichler/vol27/projects/hgsvc/nobackups/data/sequencing/HiFi/NA19240/ccs/fastq.fofn


#
# HPRC
#
HG003:
    illumina: /net/eichler/vol27/projects/sequence_data/nobackups/human/HG003/HG003.fastq.gz
HG004:
    illumina: /net/eichler/vol27/projects/sequence_data/nobackups/human/HG004/HG004.fastq.gz
HG002:
    fofn: /net/eichler/vol27/projects/sequence_data/nobackups/human/HG002/PacBioHiFi/HPRC_HiFi/downsample/fastq.fofn

HG01121:
    illumina: /net/eichler/vol26/projects/1kg_high_coverage/nobackups/data/Illumina/cram/HG01121.final.cram
HG01122:
    illumina: /net/eichler/vol26/projects/1kg_high_coverage/nobackups/data/Illumina/cram/HG01122.final.cram
HG01123:
    fofn: /net/eichler/vol27/projects/hprc/backups/data/sequencing/HiFi/HG01123/ccs/fastq.fofn

HG01256:
    illumina: /net/eichler/vol26/projects/1kg_high_coverage/nobackups/data/Illumina/cram/HG01256.final.cram
HG01257:
    illumina: /net/eichler/vol26/projects/1kg_high_coverage/nobackups/data/Illumina/cram/HG01257.final.cram
HG01258:
    fofn: /net/eichler/vol27/projects/hprc/backups/data/sequencing/HiFi/HG01258/ccs/fastq.fofn

NA12892:
    illumina: /net/eichler/vol26/projects/1kg_high_coverage/nobackups/data/Illumina/cram/NA12892.final.cram
NA12891:
    illumina: /net/eichler/vol26/projects/1kg_high_coverage/nobackups/data/Illumina/cram/NA12891.final.cram
NA12878: 
    fofn: /net/eichler/vol27/projects/hprc/nobackups/data/sequencing/HiFi/NA12878/ccs/fastq.fofn

HG01357:
    illumina: /net/eichler/vol26/projects/1kg_high_coverage/nobackups/data/Illumina/cram/HG01357.final.cram
HG01356:
    illumina: /net/eichler/vol26/projects/1kg_high_coverage/nobackups/data/Illumina/cram/HG01356.final.cram
HG01358:
    fofn: /net/eichler/vol27/projects/hprc/backups/data/sequencing/HiFi/HG01358/ccs/fastq.fofn

HG01360:
    illumina: /net/eichler/vol26/projects/1kg_high_coverage/nobackups/data/Illumina/cram/HG01360.final.cram
HG01359:
    illumina: /net/eichler/vol26/projects/1kg_high_coverage/nobackups/data/Illumina/cram/HG01359.final.cram
HG01361:
    fofn: /net/eichler/vol27/projects/hprc/backups/data/sequencing/HiFi/HG01361/ccs/fastq.fofn

HG01889:
    illumina: /net/eichler/vol26/projects/1kg_high_coverage/nobackups/data/Illumina/cram/HG01889.final.cram
HG01890:
    illumina: /net/eichler/vol26/projects/1kg_high_coverage/nobackups/data/Illumina/cram/HG01890.final.cram
HG01891:
    fofn: /net/eichler/vol27/projects/hprc/backups/data/sequencing/HiFi/HG01891/ccs/fastq.fofn

HG02256:
    illumina: /net/eichler/vol26/projects/1kg_high_coverage/nobackups/data/Illumina/cram/HG02256.final.cram
HG02255:
    illumina: /net/eichler/vol26/projects/1kg_high_coverage/nobackups/data/Illumina/cram/HG02255.final.cram
HG02257:
    fofn: /net/eichler/vol27/projects/hprc/backups/data/sequencing/HiFi/HG02257/ccs/fastq.fofn

HG02485:
    illumina: /net/eichler/vol26/projects/1kg_high_coverage/nobackups/data/Illumina/cram/HG02485.final.cram
HG02484:
    illumina: /net/eichler/vol26/projects/1kg_high_coverage/nobackups/data/Illumina/cram/HG02484.final.cram
HG02486:
    fofn: /net/eichler/vol27/projects/hprc/backups/data/sequencing/HiFi/HG02486/ccs/fastq.fofn

HG02558:
    illumina: /net/eichler/vol26/projects/1kg_high_coverage/nobackups/data/Illumina/cram/HG02558.final.cram
HG02557:
    illumina: /net/eichler/vol26/projects/1kg_high_coverage/nobackups/data/Illumina/cram/HG02557.final.cram
HG02559:
    fofn: /net/eichler/vol27/projects/hprc/backups/data/sequencing/HiFi/HG02559/ccs/fastq.fofn

HG02571:
    illumina: /net/eichler/vol26/projects/1kg_high_coverage/nobackups/data/Illumina/cram/HG02571.final.cram
HG02570:
    illumina: /net/eichler/vol26/projects/1kg_high_coverage/nobackups/data/Illumina/cram/HG02570.final.cram
HG02572:
    fofn: /net/eichler/vol27/projects/hprc/backups/data/sequencing/HiFi/HG02572/ccs/fastq.fofn


# establish the trio structure for each sample
trios:
    NA12878:
        mat: NA12892
        pat: NA12891
    HG00733: 
        mat: HG00732
        pat: HG00731
    HG00514:
        mat: HG00513
        pat: HG00512
    NA19240:
        mat: NA19238
        pat: NA19239
    HG01123:
        mat: HG01122
        pat: HG01121
    HG002:
        mat: HG004
        pat: HG003
    HG01258:
        mat: HG01257
        pat: HG01256
    HG01358:
        mat: HG01357
        pat: HG01356
    HG01361:
        mat: HG01360
        pat: HG01359
    HG01891:
        mat: HG01889
        pat: HG01890
    HG02257:
        mat: HG02256
        pat: HG02255
    HG02486:
        mat: HG02485
        pat: HG02484
    HG02559:
        mat: HG02558
        pat: HG02557
    HG02572:
        mat: HG02571
        pat: HG02570
```

