import csv
from os import path

import matplotlib.pyplot as plt

def process_overhang_csv(overhang_csv, end = 'five', library_size = None, add_overhang_csv = None, add_library_size = None):
    valid_ends = {'five', 'three'}
    if end not in valid_ends:
        raise ValueError(f'Invalid "end" option, select one of:\n\t {", ".join(valid_ends)}')


    if end == 'five':
        end_col = 1
    else:
        end_col = 2
    with open(overhang_csv) as csvin:
        csvreader = csv.reader(csvin)
        header_info = next(csvreader)
        xvals, yvals = process_lines(csvreader, end_col, library_size)

    if add_overhang_csv is None:
        return xvals, yvals

    with open(add_overhang_csv) as csvin:
        print('#'*10)
        csvreader = csv.reader(csvin)
        header_info = next(csvreader)
        add_xvals, add_yvals = process_lines(csvreader, end_col, add_library_size)

    xvals = [xvals, add_xvals]
    yvals = [yvals, add_yvals]

    return xvals, yvals


def process_lines(csvreader, end_col, library_size = None):
    xvals = []
    yvals = []
    for line in csvreader:
        x = int(line[0])
        y = int(line[end_col])
        if library_size is not None:
            y = normalise(y, library_size)
        xvals.append(x)
        yvals.append(y)
    return xvals, yvals

def normalise(count, library_size):
    scaling_factor = library_size / 1000000
    print(count / scaling_factor)
    return count / scaling_factor

def make_plot(xvals, yvals, figname, width):
    try:
        plt.bar(xvals, yvals, width)
        plt.xticks(xvals)
        plt.savefig(figname)
    finally:
        plt.close()

def make_double_plot(xvals, yvals, figname, width, legend_names):
    try:
        plt.bar(xvals[0], yvals[0], align = 'edge', width = width)
        plt.bar(xvals[1], yvals[1], align = 'edge', width = -width)
        plt.xticks(xvals[0])
        plt.legend(legend_names)
        plt.savefig(figname)
    finally:
        plt.close()

def main(overhang_csv, end, figname, library_size, width, add_overhang_csv, add_library_size):
    xvals, yvals = process_overhang_csv(overhang_csv, end, library_size, add_overhang_csv, add_library_size)
    if add_overhang_csv is None:
        make_plot(xvals, yvals, figname, width)
    else:
        legend_1 = path.splitext(path.basename(overhang_csv))[0]
        legend_2 = path.splitext(path.basename(add_overhang_csv))[0]
        make_double_plot(xvals, yvals, figname, width, (legend_1, legend_2))

if __name__ == '__main__':
    from argparse import ArgumentParser

    parser = ArgumentParser('Plot a stepRNA overhang output csv')

    parser.add_argument('overhang_csv', help = 'Path to the overhang csv to plot')
    parser.add_argument('output', help = 'Figure output filepath')
    parser.add_argument('-d', '--duplex_end', help = 'End of duplex to plot. Chose from: five or three',
            default = 'five', choices = {'five', 'three'})
    parser.add_argument('-s', '--library_size', help = 'Optional RPM normalisation using library size', 
            default = None, type = int)
    parser.add_argument('-w', '--width', help = 'Bar plot width', default = 0.8, type = float)
    parser.add_argument('-a', '--add_overhang_csv', help = 'Additonal overhang csv to plot alongside', 
            default = None)
    parser.add_argument('-S', '--add_library_size', help = 'Optional RPM normalisation for additional overhang csv using library size', 
            default = None, type = int)

    args = parser.parse_args()

    if (args.library_size is not None) and (args.add_library_size is None):
        raise ValueError(f'Please specify an additional library size using -S')

    main(args.overhang_csv, args.duplex_end, args.output, args.library_size, args.width, args.add_overhang_csv, args.add_library_size)
