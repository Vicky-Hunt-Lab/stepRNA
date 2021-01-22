
# Code to test out parameters for bowtie2

# Using the mock data I have made:
# ~50% should align and pass the test
# May be a little more due to random sequences

# Length is 20 nucleotides with a 3 base overhang

# Therefore the 'best' score should be set to 17 for local alignment

bowtie2-build dsRNA_ref_seqs.fasta reference

bowtie2 -x reference \
-U dsRNA_lib_seqs.fasta \
-f -N 0 -L 10 \
--local \
--ma 2 --mp 25,25 --rdg 25,1 --rfg 25,1 \
--score-min L,1,0 \
-S dsRNA_alignment.sam

# The above worked with 48 aligning 0 times

bowtie2 -x reference \
-U dsRNA_lib_seqs.fasta \
-f -N 0 -L 10 \
--local \
--ma 2 --mp 25,25 --rdg 25,1 --rfg 25,1 \
--score-min G,17,0 \
-S dsRNA_alignment.sam

bowtie2 -x reference \
-U dsRNA_lib_seqs.fasta \
-f -N 0 -L 17 \
-a -D 100 -R 10 \
--no-1mm-upfront \
--local \
--ignore-quals \
--ma 2 --mp 25,25 --score-min L,17,0 \
--nofw \
-S dsRNA_alignment.sam

bowtie2 -x reference \
-U dsRNA_lib_seqs.fasta \
-f -N 0 -L 10 \
--no-1mm-upfront \
--local \
--ignore-quals \
--ma 2 --mp 25,25 --score-min L,17,0 \
--nofw \
-S dsRNA_alignment.sam



# What seems to be working

bowtie2 -x reference -U dsRNA_lib_seqs.fasta \
-f -N 0 -L 17 \
--nofw --no-1mm-upfront \
--local --ma 3 --score-min L,51,0 \
-S dsRNA_alignment.sam

