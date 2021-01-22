
# Code to test out parameters for bowtie2

# Using the mock data I have made:
# ~50% should align and pass the test
# May be a little more due to random sequences

# Length is 20 nucleotides with a 3 base overhang

# Therefore the 'best' score should be set to 17 for local alignment

bowtie2-build dsRNA_ref_seqs.fasta reference

# What seems to be working

bowtie2 -x reference -U dsRNA_lib_seqs.fasta \
-f -N 0 -L 17 \
--no-1mm-upfront \
--local --ma 3 --score-min L,51,0 \
-S dsRNA_alignment.sam

