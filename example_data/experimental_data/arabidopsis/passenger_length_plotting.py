import matplotlib.pyplot as plt
import pandas
import numpy


def extract_passenger_lengths(pass_len_file):
    pass_len_df = pandas.read_csv(pass_len_file)
    pass_len_dict = {}
    for length, count in zip(pass_len_df['passenger_length'], pass_len_df['passenger_count']):
        pass_len_dict[length] = count

    return pass_len_dict

def plot_passenger_lengths(pass_len_dicts, colors, legend_names = [], width = 0.35, N = None):
    fig, ax = plt.subplots()
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
    plt.ylabel('Count')

    return fig




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







