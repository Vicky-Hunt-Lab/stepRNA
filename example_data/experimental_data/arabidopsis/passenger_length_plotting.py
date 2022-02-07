import matplotlib.pyplot as plt
import pandas
import numpy
import pysam
from Bio import SeqIO


def extract_passenger_lengths(pass_len_file):
    pass_len_df = pandas.read_csv(pass_len_file)
    pass_len_dict = {}
    for length, count in zip(pass_len_df['passenger_length'], pass_len_df['passenger_count']):
        pass_len_dict[length] = count

    return pass_len_dict

def plot_passenger_lengths(pass_len_dicts, colors, legend_names = [], width = 0.35, N = None, ylabel = None):
    fig, ax = plt.subplots()
    if ylabel is None:
        ylabel = 'Count'
    if N is None:
        N = len(pass_len_dicts[0].keys())
    ind = numpy.arange(N)
    for i, pass_len_dict in enumerate(pass_len_dicts):
        #i = i + 1
        x =  sorted(list(pass_len_dict.keys()))
        y = []
        for length in x:
            y.append(pass_len_dict[length])
        ax.bar(ind + (i*width), y, width, color = colors[i])

    ax.legend(legend_names)
    ax.set_xticks(ind + width / (i+1) )
    ax.set_xticklabels(x)

    plt.xlabel('Passenger length')
    plt.ylabel(ylabel)

    return fig

def extract_counts(fastafile): 
    counts = {} 
    for record in SeqIO.parse(fastafile, 'fasta'): 
        try: 
            count = int(record.description.split()[-1]) 
            counts[record.id] = int(count) 
        except ValueError:
            raise ValueError('Cannot convert description into a integer. Make sure header info is NAME then COUNT only')

    return counts

def passenger_expression(bamfile, counts, libsize):
    scaling_factor = libsize / 1000000
    passenger_expression_dict = {}
    for line in pysam.AlignmentFile(bamfile):
        passenger_name = line.query_name
        passenger_length = line.query_length
        try:
            passenger_expression_dict[passenger_length] += counts[passenger_name]
        except KeyError:
            passenger_expression_dict[passenger_length] = counts[passenger_name]
    for length, expr in passenger_expression_dict.items():
        passenger_expression_dict[length] = expr / scaling_factor
    return passenger_expression_dict


def main(pass_len_files, figfilename = None):
    if figfilename is None:
        figfilename = 'passengerlengths.pdf'
    pass_len_dicts = []
    for pass_len_file in pass_len_files:
        pass_len_dicts.append(extract_passenger_lengths(pass_len_file))

    pass_len_plot = plot_passenger_lengths(pass_len_dicts, ['black', 'red'], legend_names = ['WT', 'DCL Mutant'])
    pass_len_plot.savefig(figfilename)


if __name__ == '__main__':
    pass_len_files = ['stepRNAoutput/miRNAfilter/WT_passenger_length.csv','stepRNAoutput/miRNAfilter/DCL_passenger_length.csv']

    main(pass_len_files, 'WT_DCL_passengerlengths.pdf')


    # passenger expression
    fastafile = 'reads/GSM1845210_trim_miRNAfilter_collapsed.fasta'
    bamfile = 'stepRNAoutput/miRNAfilter/WT_AlignmentFiles/WT_passed.bam'
    libsize = 8443883
    counts = extract_counts(fastafile)
    passenger_expression_dict = passenger_expression(bamfile, counts, libsize)
    fastafile2 = 'reads/GSM1845222_trim_miRNAfilter_collapsed.fasta'
    bamfile2 = 'stepRNAoutput/miRNAfilter/DCL_AlignmentFiles/DCL_passed.bam'
    libsize2 = 7252698
    counts2 = extract_counts(fastafile2)
    passenger_expression_dict2 = passenger_expression(bamfile2, counts2, libsize2)
    passenger_expression_plot = plot_passenger_lengths([passenger_expression_dict, passenger_expression_dict2], colors = ['black', 'red'], legend_names = ['WT', 'DCL Mutant'], ylabel = 'Expression (TPM)')
