# 1 Download the required files

# Reads
# miRNA sequences
wget https://www.mirbase.org/ftp/CURRENT/mature.fa.gz -P feat_files
gunzip feat_files/mature.fa.gz
#Extract thaliana fasta records only, any way you want
#feat_files/thaliana_miRNA.fasta

# 2 Filter for 15 - 30 nt length sRNAs
#Already pre-filtered to 8443883 reads

# 3 Filter out any miRNA sequences
python3 ~/stepRNA/stepRNA/miRNAfilter.py \
	--nonmatches \
	reads/GSM1845210_trim.fasta \
	feat_files/thaliana_miRNA.fasta \
	reads/GSM1845210_trim_miRNAfilter.fasta 
#238269 miRNA reads removed; 8205614 remaining

# Collapse (using NGStoolbox 

perl NGStoolbox/TBr2_collapse.pl \
	-i reads/GSM1845210_trim_miRNAfilter.fasta \
	-o reads/GSM1845210_trim_miRNAfilter_collapsed.fasta
# 2908404 sequences left after collapsing

# 4 Run stepRNA
stepRNA --reference reads/GSM1845210_trim_miRNAfilter_collapsed.fasta \
	--reads reads/GSM1845210_trim_miRNAfilter_collapsed.fasta \
	--directory stepRNAoutput/miRNAfilter \
	--name WT \
	--make_unique 

# 1186873 aligned
# 376924 with a 24nt passenger
# 3' overhang is between 0 and 2 (173325, 142424, 107354)
# 5' overhang is either side of 0q

### Repeat for DCL mutant

# 7252698 reads input
python3 ~/stepRNA/stepRNA/miRNAfilter.py \
	--nonmatches \
	reads/GSM1845222_trim.fasta \
	feat_files/thaliana_miRNA.fasta \
	reads/GSM1845222_trim_miRNAfilter.fasta 
#744377 miRNA reads removed; 6508321 remaining

# Collapse (using NGStoolbox) 

perl ../NGStoolbox/TBr2_collapse.pl \
	-i reads/GSM1845222_trim_miRNAfilter.fasta \
	-o reads/GSM1845222_trim_miRNAfilter_collapsed.fasta
# 2733702 sequences left after collapsing

# 4 Run stepRNA
stepRNA --reference reads/GSM1845222_trim_miRNAfilter_collapsed.fasta \
	--reads reads/GSM1845222_trim_miRNAfilter_collapsed.fasta \
	--directory stepRNAoutput/miRNAfilter \
	--name DCL \
	--make_unique 
# 1113644 aligned


#Repeat on 24nt sRNAs only, observed to be of interest...

#WT
perl NGStoolbox/TBr2_length-filter.pl \
	-i GSM1845210_trim_miRNAfilter_collapsed.fasta \
	-o GSM1845210_trim_miRNAfilter_24nt_collapsed.fasta \
	-min 24 -max 24

stepRNA --reference reads/GSM1845210_trim_miRNAfilter_24nt_collapsed.fasta \
	--reads reads/GSM1845210_trim_miRNAfilter_collapsed.fasta \
	--min_score 14 \
	--directory stepRNAoutput/miRNAfilter \
	--name WT_24nt 

#DCL
perl NGStoolbox/TBr2_length-filter.pl \
	-i GSM1845222_trim_miRNAfilter_collapsed.fasta \
	-o GSM1845222_trim_miRNAfilter_24nt_collapsed.fasta \
	-min 24 -max 24

stepRNA --reference reads/GSM1845222_trim_miRNAfilter_24nt_collapsed.fasta \
	--reads reads/GSM1845222_trim_miRNAfilter_collapsed.fasta \
	--min_score 14 \
	--directory stepRNAoutput/miRNAfilter \
	--name DCL_24nt
