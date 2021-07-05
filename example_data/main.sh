#!/usr/bin/env bash
#Data is collected and named sensibly in data/
#If trying to replicate cd into the example_data/ directory and run the following commands


#### Ref Uncollapsed; Reads Collapsed; No Spike ####

stepRNA -r biological_data/26G_uncol.fasta -q biological_data/Pass_col.fasta -d Uncol_NoSpike -n 26G
stepRNA -r biological_data/22G_uncol.fasta -q biological_data/Pass_col.fasta -d Uncol_NoSpike -n 22G

#### Ref Collapsed; Reads Collapsed; No Spike ####

stepRNA -r biological_data/26G_col.fasta -q biological_data/Pass_col.fasta -d Col_NoSpike -n 26G
stepRNA -r biological_data/22G_col.fasta -q biological_data/Pass_col.fasta -d Col_NoSpike -n 22G

#Randomly select a reference sequence from each file
get_seeded_random()
{
  seed="$1"
  openssl enc -aes-256-ctr -pass pass:"$seed" -nosalt \
    </dev/zero 2>/dev/null
}

samtools view Uncol_NoSpike/26G_AlignmentFiles/26G_passed.bam | shuf -n 1 --random-source=<(get_seeded_random 100)
samtools view Uncol_NoSpike/22G_AlignmentFiles/22G_passed.bam | shuf -n 1 --random-source=<(get_seeded_random 100)

grep -A 1 'Ref_151559' 26G_uncol.fasta
grep -A 1 'Ref_54200' 22G_uncol.fasta
#26G = Ref_151559 ; GTACAAGCGTGCTGGCGTCGAACAGA 
#22G = Ref_54200 ; GGAAAAGATAACCAGTGTCATT
#These were added to make_ref_main.py

### Break ###
# Now run make_ref_main.py to generate the spike in datasets
mkdir spikein_data
python3 -m make_ref_main 


#Generate queries and reference SPIKE-IN files
cat biological_data/26G_uncol.fasta spikein_data/26G_uncol_Refspike.fasta > spikein_data/26G_uncol_spike.fasta
cat biological_data/22G_uncol.fasta spikein_data/22G_uncol_Refspike.fasta > spikein_data/22G_uncol_spike.fasta
cat biological_data/Pass_col.fasta spikein_data/26G_uncol_Passspike.fasta > spikein_data/26G_Pass_combined.fasta
cat biological_data/Pass_col.fasta spikein_data/22G_uncol_Passspike.fasta > spikein_data/22G_Pass_combined.fasta

cat biological_data/26G_col.fasta spikein_data/24_col_Refspike.fasta > spikein_data/26G_col_spike.fasta
cat biological_data/22G_col.fasta spikein_data/24_col_Refspike.fasta > spikein_data/22G_col_spike.fasta
cat biological_data/Pass_col.fasta spikein_data/24_Passspike.fasta > spikein_data/24_Pass_combined.fasta


### Run stepRNA ### 

stepRNA -r spikein_data/26G_uncol_spike.fasta -q spikein_data/26G_Pass_combined.fasta -d Uncol_Spike -n 26G
stepRNA -r spikein_data/22G_uncol_spike.fasta -q spikein_data/22G_Pass_combined.fasta -d Uncol_Spike -n 22G


#### Ref Collapsed; Reads Collapsed; Spike ####

stepRNA -r spikein_data/26G_col_spike.fasta -q spikein_data/24_Pass_combined.fasta -d Col_Spike -n 26G
stepRNA -r spikein_data/22G_col_spike.fasta -q spikein_data/24_Pass_combined.fasta -d Col_Spike -n 22G

