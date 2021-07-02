from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import random
import numpy as np
from collections import defaultdict

class SequenceGen():
    def __init__(self, libsize, nreferences):
        self.libsize = libsize
        self.nref = nreferences
        self.seqs = [] 
        self.refsused = set()
        self.querysamples = []
        self.queryseqs = []
        self.fiveprimeinfo = defaultdict(lambda : 0)
        self.threeprimeinfo = defaultdict(lambda : 0)

    def nuc_seq(self, length, gc_bias=0.5):
        '''Create a DNA/RNA sequence of length and with a gc bias (keep to 1% points)'''
        at_bias = gc_bias*100
        gc_bias = (1-gc_bias) * 100
        lcm = np.lcm(int(at_bias), int(gc_bias))
        at_bias = int(lcm / at_bias)
        gc_bias = int(lcm / gc_bias)
        return Seq(''.join(random.choice('CG'*gc_bias + 'TA'*at_bias) for _ in range(length)))

    def calc_gcbias(self, sequence):
        '''Calcualte the GC bias of a sequene'''
        at = sequence.count('A') + sequence.count('T')
        gc = sequence.count('G') + sequence.count('C')
        return at, gc

    def makelibrary(self, length, gc_bias=0.5):
        for seq in range(self.libsize):
             self.seqs.append(self.nuc_seq(length, gc_bias))
        print('Library generated: {} sequences'. format(self.libsize))

    def refpool(self, baseexclude=None):
        exclude_base = []
        if baseexclude != None:
            for seq in self.seqs:
                if seq[0] not in baseexclude:
                    exclude_base.append(seq)
            print('Exclude base list length:')
            print(len(exclude_base))
            try:
                self.refsamples = random.sample(exclude_base, self.nref)
            except ValueError:
                print('Excluding starting bases has failed. Please try again without excluding bases or set a smaller reference sample size.')
        else:
            self.refsamples = random.sample(self.seqs, self.nref)


    def custom_refpool(self, reflist):
        self.refsamples = []
        for seq in reflist:
            self.refsamples.append(Seq(seq))

    def makeoverhangquery(self, cutinfo, nqueries):
        '''Cutinfo contains tuples of (length, end) for 5prime and 3prime ends'''
        samples = []
        try:
            samples = random.sample(self.refsamples, nqueries)
        except ValueError:
            print('Replicating samples')
            while len(samples) < nqueries:
                samples.extend(random.sample(self.refsamples, 1))
        self.refsused.update(set(samples))
        self.fiveprimeinfo[cutinfo[0][0]] += nqueries 
        self.threeprimeinfo[cutinfo[1][0]] += nqueries
        def cutend(sample, overhanglength, overhangend):
            if overhangend != '3prime':
                return sample[:overhanglength]
            else:
                return sample[-overhanglength:]
        def addend(sample, gnomsample, overhanglength, overhangend):
            addbases = self.nuc_seq(overhanglength)
            if overhangend != '3prime':
                return (sample + addbases, addbases.reverse_complement() + gnomsample)
            else:
                return (addbases + sample, gnomsample + addbases.reverse_complement())
        for x, refsample in enumerate(samples):
            gnomsample = refsample
            querysample = refsample.reverse_complement()
            for end in cutinfo:
                overhanglength = end[0]
                overhangend = end[1]
                if overhanglength < 0:
                    querysample = cutend(querysample, overhanglength, overhangend)
                elif overhanglength > 0:
                    querysample, gnomsample = addend(querysample, gnomsample, overhanglength, overhangend)
            self.querysamples.append(querysample)
            self.queryseqs.append(gnomsample)

        self.overhangprop = nqueries / self.libsize
        
    def makerefgenome(self):
        seqs = set(self.seqs)
        seqsnotused = seqs.difference(self.refsused) 
        seqs = list(seqsnotused)
        random.shuffle(seqs)
        hits = self.queryseqs
        random.shuffle(hits)
        def makeseq(hits, seqs):
            genome = Seq('')
            hitslen, seqslen = len(hits), len(seqs)
            total = hitslen + seqslen
            for _ in range(total):
                if len(hits) > 0:
                    if len(seqs) > 0:
                        if random.randint(1, total) < hitslen:
                            genome += hits.pop() + self.nuc_seq(5)
                        else:
                            genome += seqs.pop() + self.nuc_seq(5)
                    else:
                        genome += hits.pop() + self.nuc_seq(5)
                else:
                    genome += seqs.pop() + self.nuc_seq(5)
            return genome 

        return makeseq(hits, seqs)

    def sendtofile(self, genome_filename, sRNA_filename, refRNA_filename, fileformat='fasta'):
        genome = SeqRecord(self.makerefgenome())
        genome.id = 'SimulatedData_ReferenceGenome'
        genome.description = ''
        genome.format(fileformat)

        def makeSeqRecord(srna, fileformat):
            random.shuffle(srna)
            for x, seq in enumerate(srna):
                seq = SeqRecord(seq)
                seq.id = 'Read_' + str(x+1)
                seq.description = ''
                seq.format(fileformat)
                srna[x] = seq
            return srna

        srna = self.seqs + self.querysamples
        print('sRNA count:', len(srna))
        refrna = self.refsamples
        srna = makeSeqRecord(srna, fileformat)
        refrna = makeSeqRecord(refrna, fileformat)

        with open(genome_filename, 'w') as Genome, open(sRNA_filename, 'w') as sRNA, open(refRNA_filename, 'w') as refRNA:
            SeqIO.write(genome, Genome, fileformat)
            SeqIO.write(srna, sRNA, fileformat)
            SeqIO.write(refrna, refRNA, fileformat)
            
def makeoverhangs(seqlib, querysamplesize, oh):
    seqlib.makeoverhangquery(oh, querysamplesize)
    print('Added a {}nt {}; {}nt {} overhang library with {}% matches.'.format(oh[0][0], oh[0][1], oh[1][0], oh[1][1], seqlib.overhangprop * 100))
    return seqlib

def main(libsize, refsamplesize, seqlength, overhang_info, refname, sRNAname, refRNAname, seed=None, baseexclude=None):
    random.seed(seed)
    seqlib = SequenceGen(libsize, refsamplesize)
    seqlib.makelibrary(seqlength)
    seqlib.refpool(baseexclude)
    for nsample, overhangs in overhang_info:
        makeoverhangs(seqlib, nsample, overhangs)
    seqlib.sendtofile(refname, sRNAname, refRNAname)
    return seqlib

if __name__ == '__main__':
    #overhang_info = [[10, ((0, '3prime'), (-2, '5prime'))], [50, ((2, '3prime'), (-2, '5prime'))]]
    overhang_info = [[50, ((0, '3prime'), (0, '5prime'))]]
    _ = main(100, 50, 15, overhang_info, 'refgenome_test.fasta', 'sRNA_test.fasta', 'refRNA_test.fasta') 
