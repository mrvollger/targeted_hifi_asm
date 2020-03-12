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
