# Data aquisition and pre-processing

## Contents:
- [Experimental Data](#Experimental-Data)
- [Simulated Data](#Simulated-Data)

## Experimental Data

### Files

The files contained within ```experimental_data/``` include:
- *LF_embryo.fa*; small RNA passengers length filtered to 15-30nt
- *26G_embryo.fa*; 26G sRNA sequences. Length filtered to 26 starting with a G
- *22G_embryo.fa*; 22G sRNA sequences. Length filtered to 22 starting with a G

### Aquisition and Processing

1) Wild-Type Embryo Data Table (Sample ID: [GSM801363](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM801363)); GEO Accession Number: [GSE32336](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE32366) was downloaded as *GSM801363-2802.tab*

The table has been processed according to [Fischer et al. (2011)](https://pubmed.ncbi.nlm.nih.gov/22102828/) and contains unique sequences and counts obtained from the FASTQ seqeunceing file (filtered for Q>20)

2) Header details were removed and sequences converted to FASTA format

```
tail -n +4 GSM801363-2802.tab > GSM801363_WTembryo.tab 

awk '{print ">" NR; print $0}' GSM801363_WTembryo.tab > GSM801363_rawseqs_WTembryo.fa 
```

3) Length filtering was performed with the [NGS TOOLBOX](https://www.smallrnagroup.uni-mainz.de/software/TBr2.zip) script, TBr2_length-filter.pl, from the [small RNA group, Mainz University](https://www.smallrnagroup.uni-mainz.de/)

```
perl TBr2_length-filter.pl -i GSM801363_rawseqs_WTembryo.fa -o LF_embryo.fa -min 15 -max 30
perl TBr2_length-filter.pl -i GSM801363_rawseqs_WTembryo.fa -o 26_embryo.fa -min 26 -max 26
perl TBr2_length-filter.pl -i GSM801363_rawseqs_WTembryo.fa -o 22_embryo.fa -min 22 -max 22
```

4) 22G and 26G seqeunces were then selected using GREP
```
grep '^G' -B 1 26_embryo.fa | sed '/--/d' > 26G_embryo.fa
grep '^G' -B 1 22_embryo.fa | sed '/--/d' > 22G_embryo.fa
```

## Simulated Data

### Files

The files contained within ```simulated_data/``` include:
- *ref.fasta*; referece sequence (row 1 in *alignment.csv*)
- *reads.fasta*; query sequences (rows 2-29 in *alignment.csv*)
- [*alignment.pdf*](./simulated_data/alignment.pdf); a visual alignment of how the reads align to the genome

These can be run to test that stepRNA is working on your system and provide a small dataset to help understand the output.
