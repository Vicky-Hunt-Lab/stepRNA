from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from random import randint

ref_seq = Seq('GAACTCACTAAGATTGACGATTTTGAT')
lst = []
#3'overhang(5-8)
for read in range(5,9,1):
    lst.append(ref_seq + (read*'T'))
    #5'3nt-overhang; 3'6nt-overhang 
    if read == 6:
        lst.append('TTT' + ref_seq + (read * 'T'))


#5'overhang(3-6)
for read in range(3,3+4,1):
    lst.append((read*'T') + ref_seq)
    #5'5nt-overhang; 3'4nt-underhang
    if read == 5:
        lst.append((read*'T') + ref_seq[:-4])


#3'underhang(1-3)
for read in range(1,1+3,1):
    lst.append(ref_seq[:-read])
    

#5'underhang(4-6)
for read in range(4, 4+3, 1):
    lst.append(ref_seq[read:])
    #5'6nt-underhang; 3'3nt-overhang
    if read == 6:
        lst.append(ref_seq[read:] + 'TTT')
    #5'4nt-underhang, 3'2nt-underhang
    if read == 4:
        lst.append(ref_seq[read:-2])

#5'exact; 3'exact
lst.append(ref_seq)

for x, read in enumerate(lst):
    lst[x] = read.reverse_complement()

bases = ['G','T', 'A', 'C']
#Non-matching sequences
for read in bases:
    lst.append(Seq(read*19))


for x in range(2):
    seq_length = randint(16,30)
    seq = []
    for base in range(seq_length):
        index = randint(0,3)
        seq.append(bases[index])
    lst.append(Seq(''.join(seq)))


with open('reads_4.fasta', 'w') as fout:
    for x, read in enumerate(lst):
        record = SeqRecord(read)
        record.id = 'read_' + str(x+1)
        record.description = ''
        SeqIO.write(record, fout, 'fasta')
