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

```
