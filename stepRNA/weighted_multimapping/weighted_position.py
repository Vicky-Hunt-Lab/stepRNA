import numpy as np

def process_unique_reads(unique_read_dict):
    '''Process a dictionary of sRNAs with the number of unique reads'''
    total = sum(unique_read_dict.values())
    probabilities = []
    _ = [probabilities.append(unique_count / total) for unique_count in unique_read_dict.values()]

    return probabilities

def process_fractional_reads(fractional_read_dict):
    '''Process a dictionary of sRNAs with the number of fractionalmapped reads to that part'''
    scores = []
    for read, read_dict in fractional_read_dict.items():
        score = 0
        for mmap, count in read_dict.items():
            score += count * (1 / mmap)
        scores.append(score)
    total = sum(scores)
    probabilities = []
    _ = [probabilities.append(score / total) for score in scores]
    return probabilities

def make_weighted_array(probabilities):
    '''From a tuple of probabilities make an array to test a np.random.rand value against'''
    weighted_array = [0] 
    _ = [weighted_array.append(weighted_array[-1] + probability) for probability in probabilities]
    return np.array(weighted_array)

def primary_alignment_block(weighted_array):
    '''Make a random number and determine whether which block it fits into'''
    random_number = np.random.rand()
    #print(random_number)
    select_array = random_number < weighted_array
    return len(select_array) - select_array.sum()

def main_unique(unique_read_dict):
    '''Select a random block to put a read into from a dictionary of unique read counts'''
    probabilities = process_unique_reads(unique_read_dict)
    weighted_array = make_weighted_array(probabilities)
    block = primary_alignment_block(weighted_array)
    return block

def main_fractional(fractional_read_dict):
    '''Select a random block to put a read into from a dictionary of fractional read counts'''
    probabilities = process_fractional_reads(fractional_read_dict)
    weighted_array = make_weighted_array(probabilities)
    block = primary_alignment_block(weighted_array)
    return block

if __name__ == '__main__':
    unique_read_dict = {
            'S1' : 6,
            'S2' : 4,
            'S3' : 2
            }

    main(unique_read_dict)
