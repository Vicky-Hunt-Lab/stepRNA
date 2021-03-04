import csv

from stepRNA.output import make_csv, make_type_csv, write_to_bam, print_hist, refs_counts

def make_hist(csv_in, logger):
    with open(csv_in) as summary:
        left_dens = []
        right_dens = []
        left_tot = 0
        right_tot = 0
        keys = []
        csv_reader = csv.reader(summary, delimiter=',')
        head = next(csv_reader)
        for line in csv_reader:
            keys.append(line[0])
            left_dens.append(int(line[1]))
            left_tot += int(line[1])
            right_dens.append(int(line[2]))
            right_tot += int(line[2])

    for key in range(len(keys)):
        left_dens[key] = 100 * left_dens[key] / left_tot
        right_dens[key] = 100 * right_dens[key] / right_tot

    #Print histogram of overhangs to terminal...
    logger.write('5\' Overhang Histogram')
    print_hist(left_dens, keys, logger)

    logger.write('3\' Overhang Histogram')
    print_hist(right_dens, keys, logger)

def main(right_dic, left_dic, type_dic, read_len_dic, refs_read_dic, right_unique_dic, left_unique_dic, prefix, logger):
    logger.write('\n#### Summary Tables ####')
    logger.write('\n## Overhang counts ##')
    make_csv([right_dic, left_dic], prefix + '_overhang.csv', ['Overhang','5prime','3prime'], logger, show=True)
    logger.write('\n## Unique overhang counts ##')
    make_csv([right_unique_dic, left_unique_dic], prefix + '_unique_overhang.csv', ['Overhang','5prime','3prime'], logger, show=True)
    logger.write('\n## Overhang types ##')
    make_type_csv(type_dic, prefix + '_overhang_type.csv', ['Classification', 'count'], logger, show=True, sort=False)
    logger.write('\n## Read lengths ##')
    make_type_csv(read_len_dic, prefix + '_passenger_length.csv', ['passenger_length', 'passenger_count'], logger,sort=True)
    make_type_csv(refs_read_dic, prefix + '_passenger_number.csv', ['Passenger_length', 'number'], logger, show=False)
    logger.write('\n')

    logger.write('\nAll aligned reads')
    make_hist(prefix + '_overhang.csv', logger)
    logger.write('\nUnique aligned reads')
    make_hist(prefix + '_unique_overhang.csv', logger)

if __name__ == "__main__":
#Setup argparse with input of JSON dicts and create a logger
    main(right_dic, left_dic, type_dic, read_len_dic, refs_read_dic, right_unique_dic, left_unique_dic, prefix, logger)
