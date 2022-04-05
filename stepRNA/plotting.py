import csv

import matplotlib.pyplot as plt

def process_overhang_csv(overhang_csv, end = 'five'):
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
        xvals, yvals = process_lines(csvreader, end_col)

    return xvals, yvals


def process_lines(csvreader, end_col):
    xvals = []
    yvals = []
    for line in csvreader:
        xvals.append(int(line[0]))
        yvals.append(int(line[end_col]))
    return xvals, yvals

def make_plot(xvals, yvals, figname):
    try:
        plt.bar(xvals, yvals)
        plt.xticks(xvals)
        plt.savefig(figname)
    finally:
        plt.close()


def main(overhang_csv, end, figname):
    xvals, yvals = process_overhang_csv(overhang_csv, end)
    make_plot(xvals, yvals, figname)

if __name__ == '__main__':
    from argparse import ArgumentParser

    parser = ArgumentParser('Plot a stepRNA overhang output csv')

    parser.add_argument('overhang_csv', help = 'Path to the overhang csv to plot')
    parser.add_argument('output', help = 'Figure output filepath')
    parser.add_argument('-d', '--duplex_end', help = 'End of duplex to plot. Chose from: five or three',
            default = 'five', choices = {'five', 'three'})

    args = parser.parse_args()

    main(args.overhang_csv, args.duplex_end, args.output)
