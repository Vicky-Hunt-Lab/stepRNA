#!/usr/bin/env python3
import make_references
import random

def make_overhanginfo(overhanginfo, overhanglist):
    for n, end in overhanglist:
        overhanginfo.append([n, ((end[1], '3prime'), (end[0], '5prime'))])


if __name__ == '__main__':
    overhanglist = [
    [50,(0,-3)],
    [50,(-1,0)],
    [400,(2,-2)],
    [50,(-2,0)],
    [25,(1,-2)],
    
    [25,(-2,-4)],
    [25,(3,-2)],
    [50,(-1,-1)],
    [25,(2,-3)],
    [25,(-3,-2)],

    [25,(2,-1)],
    [25,(-3,3)],
    [15,(4,-4)],
    [75,(3,1)],
    [100,(2,0)],

    [200,(1,-1)],
    [50,(0,1)],
    [25,(-4,2)],
    [25,(-1,-2)],
    [50,(0,-2)],

    [50,(1,-3)]
    ]
    overhanginfo = []
    make_overhanginfo(overhanginfo, overhanglist)
    extended_sim = make_references.main(20000, 1000, 21, overhanginfo, 'spikein_data/refgen_extended.fasta', 'spikein_data/sRNA_extended.fasta', 'spikein_data/refRNA_extended.fasta', seed=100)

    #Spike in 24 sequence
    spike24 = [
    [5000, (-1, -2)]
    ]
    spike24info = []
    make_overhanginfo(spike24info, spike24)
    spike24seqs = make_references.main(20000, 5000, 24, spike24info, 'spikein_data/24_spike_genome.fasta','spikein_data/24_Passspike.fasta', 'spikein_data/24_col_Refspike.fasta', seed=100, baseexclude = ['G', 'T'])


    #Spike in 26G sequence (from main.sh)
    spike26 = [
    [5000, (0,-1)]
    ]
    spike26info = []
    make_overhanginfo(spike26info, spike26)
    random.seed(100)
    spike26seqs = make_references.SequenceGen(0,0)
    spike26seqs.custom_refpool(['GTACAAGCGTGCTGGCGTCGAACAGA'])
    spike26seqs.libsize = 5000
    spike26seqs.makeoverhangquery(spike26info[0][1], spike26info[0][0])
    spike26seqs.seqs.append(make_references.Seq('GTACAAGCGTGCTGGCGTCGAACAGA'))
    spike26seqs.sendtofile('spikein_data/26G_uncol_spike_genome.fasta', 'spikein_data/26G_uncol_Passspike.fasta', 'spikein_data/26G_uncol_Refspike.fasta')

    #Spike in 22G sequence (from main.sh) 
    spike22 = [
    [5000, (-2,-2)]
    ]
    spike22info = []
    make_overhanginfo(spike22info, spike22)
    random.seed(100)
    spike22seqs = make_references.SequenceGen(0,0)
    spike22seqs.custom_refpool(['GGAAAAGATAACCAGTGTCATT'])
    spike22seqs.libsize = 5000
    spike22seqs.makeoverhangquery(spike22info[0][1], spike22info[0][0])
    spike22seqs.seqs.append(make_references.Seq('GGAAAAGATAACCAGTGTCATT'))
    spike22seqs.sendtofile('spikein_data/22G_uncol_spike_genome.fasta', 'spikein_data/22G_uncol_Passspike.fasta', 'spikein_data/22G_uncol_Refspike.fasta')

