# Data aquisition and pre-processing

## Files

The files contained within this directory include:
- LF_embryo.fa; small RNA passengers length filtered to 15-30nt
- 26G_embryo.fa; 26G sRNA sequences. Length filtered to 26 starting with a G
- 22G_embryo.fa; 22G sRNA sequences. Length filtered to 22 starting with a G

## Aquisition and Processing

1) Wild-Type Embryo Data Table (Sample ID: [GSM801363](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM801363)) was downloaded from GEO Accession Number [GSE32336](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE32366) as GSM801363-2802.txt

The table has been processed according to [Fischer et al., 2011](https://pubmed.ncbi.nlm.nih.gov/22102828/) and contains unique sequences and counts obtained from the FASTQ seqeunceing file (filtered for Q>20)

2) Header details were manually removed and sequences converted to FASTA format

```
tail -n +4 GSM801363-2802.txt > GSM801363_WTembryo.txt 
awk '{print ">" NR; print $0}' GSM801363_WTembryo.txt > GSM801363_rawseqs_WTembryo.txt 
```

3) Length filtering was performed with a [NGS TOOLBOX](https://www.smallrnagroup.uni-mainz.de/software/TBr2.zip) script: TBr2_length-filter.pl

```
perl TBr2_length-filter.pl -i GSM801363_WTembryo.txt -o LF_embryo.fa -min 15 -max 30
perl TBr2_length-filter.pl -i GSM801363_WTembryo.txt -o 26G_embryo.fa -min 26 -max 26
perl TBr2_length-filter.pl -i GSM801363_WTembryo.txt -o 22G_embryo.fa -min 22 -max 22
```

4) 22G and 26G seqeunces were then selected using GREP
```
grep '^G' -B 1 fileXXX.fa | sed '/--/d' > new_fileXXX.fa
```
