import pandas
import matplotlib.pyplot as plt
import numpy


def extract_overhangs(overhang_file):
    overhang_df = pandas.read_csv(overhang_file)
    fiveoverhang_dict = {}
    threeoverhang_dict = {}
    for length, fivecount, threecount in zip(overhang_df['Overhang'], overhang_df['5prime'], overhang_df['3prime']):
        fiveoverhang_dict[length] = fivecount
        threeoverhang_dict[length] = threecount

    return fiveoverhang_dict, threeoverhang_dict

def plot_overhangs(overhang_dicts, colors, legend_names = [], width = 0.35, N = None, overhang_type = None):
    fig, ax = plt.subplots()
    if N is None:
        N = len(overhang_dicts[0].keys())
    ind = numpy.arange(N)
    for i, overhang_dict in enumerate(overhang_dicts):
        #i = i + 1
        if overhang_type is not None:
            overhang_dict = filter_overhang(overhang_dict, overhang_type)
            ind = numpy.arange(len(overhang_dict.keys()))
        x =  sorted(list(overhang_dict.keys()))
        y = []
        for length in x:
            y.append(overhang_dict[length])
        ax.bar(ind + (i*width), y, width, color = colors[i])

    ax.legend(legend_names)
    ax.set_xticks(ind + width / (i+1) )
    ax.set_xticklabels(x)

    plt.xlabel('Passenger length')
    plt.ylabel('Count')

    return fig

def filter_overhang(overhang_dict, overhang_type):
    '''Overhang type is one of underhang, overhang'''
    filter_dict = {}
    if overhang_type == 'overhang':
        for key, count in overhang_dict.items():
            if key <= 0:
                filter_dict[key] = count
    elif overhang_type == 'underhang':
        for key, count in overhang_dict.items():
            if key >= 0:
                filter_dict[key] = count
    return filter_dict


def main(overhang_files, overhang_type = None, figfilename = None):
    if figfilename is None:
        figfilename = 'overhang'
    fiveoverhang_dicts = []
    threeoverhang_dicts = []
    for overhang_file in overhang_files:
        overhangs = extract_overhangs(overhang_file)
        fiveoverhang_dicts.append(overhangs[0])
        threeoverhang_dicts.append(overhangs[1])

    fiveoverhang_plot = plot_overhangs(fiveoverhang_dicts, ['black', 'red'], legend_names = ['WT', 'DCL Mutant'], 
            overhang_type = overhang_type)
    fiveoverhang_plot.savefig(figfilename + '_5prime.pdf')

    threeoverhang_plot = plot_overhangs(threeoverhang_dicts, ['black', 'red'], legend_names = ['WT', 'DCL Mutant'],
            overhang_type = overhang_type)
    threeoverhang_plot.savefig(figfilename + '_3prime.pdf')

if __name__ == '__main__':
    overhang_files = ['stepRNAoutput/miRNAfilter/WT_overhang.csv','stepRNAoutput/miRNAfilter/DCL_overhang.csv']

    main(overhang_files, figfilename = 'WT_DCL_both')
    main(overhang_files, overhang_type = 'overhang', figfilename = 'WT_DCL_overhang')

    overhang_files = ['stepRNAoutput/miRNAfilter/WT_unique_overhang.csv','stepRNAoutput/miRNAfilter/DCL_unique_overhang.csv']

    main(overhang_files, figfilename = 'WT_DCL_unique_both')
    main(overhang_files, overhang_type = 'overhang', figfilename = 'WT_DCL_unique_overhang')
