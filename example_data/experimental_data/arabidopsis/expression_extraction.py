from os import path
import glob

from Bio import SeqIO
import pysam
import matplotlib.pyplot as plt
import numpy

def extract_counts(fastafile):
    counts = {}
    for record in SeqIO.parse(fastafile, 'fasta'):
        try:
            count = int(record.description.split()[-1])
            counts[record.id] = int(count)
        except ValueError:
            raise ValueError('Cannot convert description into a integer. Make sure header info is NAME then COUNT only')

    return counts

def extract_expression(bamfile, counts):
    expression = 0
    for line in pysam.AlignmentFile(bamfile):
        #expression += counts[line.query_name]
        expression += counts[line.reference_name]
    
    return expression

def extract_all_expressions(bamfiles, counts, libsize):
    scalingfactor = libsize / 1000000
    expressions = {}
    for bamfile in bamfiles:
        name = '_'.join(path.basename(bamfile).split('_')[-3:-1])
        expressions[name] = extract_expression(bamfile, counts) / scalingfactor

    return expressions


def make_expressions_dict(fastafile, bamfiles, libsize):
    counts = extract_counts(fastafile)
    expressions = extract_all_expressions(bamfiles, counts, libsize)

    return expressions

def plot_overhangs(expressions_dicts, colors, legend_names = [], width = 0.35, N = None):
    fig, ax = plt.subplots()
    if N is None:
        N = len(expressions_dicts[0].keys())
    ind = numpy.arange(N)
    for i, expressions_dict in enumerate(expressions_dicts):
        #i = i + 1
        p = []
        for key, count in expressions_dict.items():
            key = int(key.split('_')[-1])
            p.append((key, count))
        p.sort()
        x = [] 
        y = []
        for length, count in p:
            x.append(length)
            y.append(count)
        ax.bar(ind + (i*width), y, width, color = colors[i])

    ax.legend(legend_names)
    ax.set_xticks(ind + width / (i+1) )
    ax.set_xticklabels(x)

    plt.xlabel('3\' overhang length')
    plt.ylabel('Expression (RPM)')

    return fig

def main(inputfiles):
    expressions_dicts = [] 
    for fastafile, bamfiles, libsize in inputfiles:
        print(fastafile)
        expressions_dicts.append(make_expressions_dict(fastafile, bamfiles, libsize))

    return expressions_dicts



if __name__ == '__main__':
    #argparse shit

    #Only do with either 3' or 5' currently!!!
    fastafile_WT = 'reads/GSM1845210_trim_miRNAfilter_collapsed.fasta'
    bamfiles_WT = glob.glob('stepRNAoutput/miRNAfilter/WT_AlignmentFiles/*3prime*.bam')
    libsize_WT = 8443883
    fastafile_DCL = 'reads/GSM1845222_trim_miRNAfilter_collapsed.fasta'
    bamfiles_DCL = glob.glob('stepRNAoutput/miRNAfilter/DCL_AlignmentFiles/*3prime*.bam')
    libsize_DCL = 7252698
    
    inputfiles = [(fastafile_WT, bamfiles_WT, libsize_WT), (fastafile_DCL, bamfiles_DCL, libsize_DCL)]
    expressions_dicts = main(inputfiles)
    
    overhang_plot = plot_overhangs(expressions_dicts, ['black', 'red'], legend_names = ['WT', 'DCL Mutant'])
    overhang_plot.savefig('WT_DCL_expression.pdf')
    




def plotting_lists(expressions):
    x = []
    y = []
    p = []
    for length, count in expressions.items():
        l = int(length.split('_')[-1])
        p.append((l, count))
        p.sort()

    for l, count in p:
        x.append(l)
        y.append(count)

    return x, y
