#!/usr/bin/env python3

from Bio.Seq import Seq
from Bio import SeqIO
import random
import numpy as np

def nuc_seq(length, gc_bias=0.5):
    '''Create a DNA/RNA sequence of length and with a gc bias (keep to 1% points)'''
    at_bias = gc_bias*100
    gc_bias = (1-gc_bias) * 100
    lcm = np.lcm(int(at_bias), int(gc_bias))
    at_bias = int(lcm / at_bias)
    gc_bias = int(lcm / gc_bias)
    return Seq(''.join(random.choice('CG'*gc_bias + 'TA'*at_bias) for _ in range(length)))

def calc_gcbias(sequence):
    '''Calcualte the GC bias of a sequene'''
    at = sequence.count('A') + sequence.count('T')
    gc = sequence.count('G') + sequence.count('C')
    return at, gc


# Make sequenes to target...

def mock_dsRNA(length, hang=3, n=50):
    targets = []
    non_targets = []
    for i in range(n):
        targets.append(nuc_seq(length))
        non_targets.append(nuc_seq(length))
    def make_comp(sequence, overhang=3):
        sequence = sequence[:-overhang]
        sequence = sequence.reverse_complement()
        sequence = sequence + nuc_seq(overhang)
        return sequence
    library_targets = []
    library_non_targets = []
    for seq in targets:
        library_targets.append(make_comp(seq, overhang=overhang))
        library_non_targets.append(nuc_seq(length))
    return targets, non_targets, library_targets, library_non_targets
    
targets, non_targets, library_targets, library_non_targets = mock_dsRNA(20)

ref_targets = []
ref_non_targets = []
name_lib_targets = []
name_lib_non_targets = []
for x in range(len(targets)):
    ref_targets.append('reference_target' + str(x+1))
    ref_non_targets.append('reference_non-target' + str(x+1))
    name_lib_targets.append('library_target' + str(x+1))
    name_lib_non_targets.append('library_non-target' + str(x+1))


with open('dsRNA_ref_seqs.fasta', 'w') as fout:
    for i in range(len(targets)):
        fout.write('>{}\n'.format(ref_targets[i]))
        fout.write('{}\n'.format(targets[i]))
    for i in range(len(non_targets)):
        fout.write('>{}\n'.format(ref_non_targets[i]))
        fout.write('{}\n'.format(non_targets[i]))

with open('dsRNA_lib_seqs.fasta', 'w') as fout:
    for i in range(len(targets)):
        fout.write('>{}\n'.format(name_lib_targets[i]))
        fout.write('{}\n'.format(library_targets[i]))
    for i in range(len(non_targets)):
        fout.write('>{}\n'.format(name_lib_non_targets[i]))
        fout.write('{}\n'.format(library_non_targets[i]))



