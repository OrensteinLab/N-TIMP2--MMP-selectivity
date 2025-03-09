import numpy as np
import pandas as pd
import os

def fastq_reader(filename):
    """
    The function receives a fastq file name and reads it into a dictionary.

    :param filename: string - fastq file name
    :return: dictionary of the sequences in the file
    """

    n = 4   # The number of repetitive lines in the file
    lines_keys = ['Name', 'Sequence', 'Optional', 'Quality']
    seq_dict = {}

    with open(filename, 'r') as fastq_file:
        lines = []
        for line in fastq_file:
            lines.append(line.rstrip())
            if len(lines) == n:
                d = dict(zip(lines_keys, lines))
                seq_dict[d['Name'].split(' ')[0]] = d['Sequence'] #the index of the read
                lines = []
    return seq_dict


def sort_seqs_to_mut(seq_dict, info):
    """
    The function receives dictionary of sequences and returns dictionary sorted by the valid sequences into mutations.

    :param  seq_dict: dictionary containing all the sequences
            info: list containing all the desired information:
                start_motif: string- after the motif the desired sequence will begin
                start_seq_pos: string- after this length of nucleotide the sequence starts
                nuc_length: int- the length of the desired sequence
                WTaa: string- The amino acid sequence of the WT
    :return: dictionary of sequences- the key: The name of the seq
                                      the values: aa_seq - amino acid seq, count_diff- number of differences,
                                       diff- string , mut_pos- list of mutation positions
    """

    start_motif = info[0]
    start_seq_pos = info[1]
    nuc_length = info[2]
    WTaa = info[3]
    mutation_dict = dict()
    for name, seq in seq_dict.items():
        start_motif_pos = seq.find(start_motif)
        position = start_motif_pos + start_seq_pos
        if start_motif_pos != -1 and position + nuc_length <= len(seq):
            relevant_seq = seq[position: position + nuc_length]
            aa_seq = translate_nuc_to_aa(relevant_seq)
            stop_codon = aa_seq.find('*')
            unknown_codon = aa_seq.find('?')
            if stop_codon == -1 and unknown_codon == -1:    # If the seq is valid
                count_diff, diff, mut_pos = compare_seq(seq1=WTaa, seq2=aa_seq)
                mutation_dict[name] = aa_seq, count_diff, diff, mut_pos
    return mutation_dict

def translate_nuc_to_aa(seq):
    """
    The function receives a nucleotide sequence and translates it to an amino acid sequence.

    :param seq: string of nucleotide sequence
    :return: string of amino acid sequence
    """
    dictionary = {
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'AGA': 'R', 'AGG': 'R', 'AAC': 'N', 'AAT': 'N',
        'GAC': 'D', 'GAT': 'D', 'TGC': 'C', 'TGT': 'C',
        'CAA': 'Q', 'CAG': 'Q', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'CAC': 'H', 'CAT': 'H', 'ATA': 'I', 'ATC': 'I',
        'ATT': 'I', 'CTA': 'L', 'CTC': 'L', 'CTG': 'L',
        'CTT': 'L', 'TTA': 'L', 'TTG': 'L', 'AAA': 'K',
        'AAG': 'K', 'ATG': 'M', 'TTC': 'F', 'TTT': 'F',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'AGC': 'S', 'AGT': 'S', 'TCA': 'S', 'TCC': 'S',
        'TCG': 'S', 'TCT': 'S', 'ACA': 'T', 'ACC': 'T',
        'ACG': 'T', 'ACT': 'T', 'TGG': 'W', 'TAC': 'Y',
        'TAT': 'Y', 'GTA': 'V', 'GTC': 'V', 'GTG': 'V',
        'GTT': 'V', 'TAA': '*', 'TAG': '*', 'TGA': '*'
    }
    aa_seq = ""
    if len(seq) % 3 == 0:
        for i in range(0, len(seq), 3):
            codon = seq[i:i + 3]
            if codon not in dictionary:
                aa_seq += '?'
            else:
                aa_seq += dictionary[codon]
    return aa_seq


def sort_mut_by_number(mut_dict, num_of_mut, col_name):
    """
    The function receives a dictionary of sequences, a number of mutations and a desired name
    and returns a series of read-counts for the mutations.

    :param mut_dict: dictionary of sequences
    :param num_of_mut: int - number of mutation
    :param col_name: string - the desired name of the column
    :return: a series of read-counts for the mutations with the number of mutations
    """

    df = pd.DataFrame.from_dict(mut_dict, orient='index', columns=[col_name, 'Number of Mutations',
                                                                   'Type of Mutation', 'Position'])
    if num_of_mut <= 7:
        ind_relevant_mut = df['Number of Mutations'] == num_of_mut
    else:
        ind_relevant_mut = df['Number of Mutations'] >= num_of_mut
    relevant_mut = df.loc[ind_relevant_mut, col_name].value_counts()
    return relevant_mut


def oneHot(string):
    """
    The function receives a sequence, and encode it according to one hot.

    :param string: string of the sequence
    :return: one_hot_encoding of the string
    """

    aa_dict = {"A": "0", "R": "1", "N": "2", "D": "3", "C": "4", "Q": "5", "E": "6", "G": "7", "H": "8",
               "I": "9", "L": "10", "K": "11", "M": "12", "F": "13", "P": "14", "S": "15", "T": "16", "W": "17",
               "Y": "18", "V": "19"}
    int_encode = [aa_dict[aa] for aa in string]
    one_hot_encode = []
    for i in int_encode:
        l = [0 for _ in range(20)]
        l[int(i)] = 1
        one_hot_encode.append(l)
    return one_hot_encode


def compare_seq(seq1, seq2):
    """
    The function receives two sequences and compares them.

    :param seq1: string of sequence
    :param seq2: string of sequence
    :return: the function returns 3 param: count_diff - int - the number of differences found between the 2 sequences,
                                           diff - string - the differences,
                                           position- list of int - the positions of the differences in the sequences
    """

    count_diff = 0
    diff = ''
    position = []
    if len(seq1) == len(seq2):
        for i in range(len(seq1)):
            if seq1[i] != seq2[i]:
                count_diff += 1
                diff = diff + seq1[i] + str(i + 1) + seq2[i]
                position.append(i + 1)
        if count_diff == 0:
            diff = 'WT'
    return count_diff, diff, position


def is_interface_mut(seq, WTseq, interface_mut):
    """
    The function receives a sequence, compares it to the WT sequence and returns True if the differences
    in interface positions.

    :param seq: string of the sequence
    :param WTseq: string of the WT sequence
    :param interface_mut: list of int - the interface positions
    :return: True- if the differences are in interface positions
             False - if the differences are not in interface positions
    """

    _, _, diff = compare_seq(seq1=WTseq, seq2=seq)
    if all(position in interface_mut for position in diff):
        return True
    return False


def create_res_df_filter(data_path, relevant_columns):
    """
    The function takes the total count from each csv files, and merge all them into a single dataframe.

    :param data_path: string - fastq file name
    :param relevant_columns - string - fastq file name
    :return: dataframe of the total count of each gate
    """

    WTaa = 'CSCSPVHPQQAFCNADVVIRAKAVSEKEVDSGNDIYGNPIKRIQYEIKQIKMFKGPEKDIEFIYTAPSSAVCGVSLDVGGKKEYLIAGKAEGDGKMHIT'
    file_names = ['Pre', 'MMP1_L', 'MMP1_H', 'MMP3_L', 'MMP3_H']
    final_names = ['pre', 'MMP-1 G1', 'MMP-1 G4', 'MMP-3 G1', 'MMP-3 G4']

    res_df = None
    for index, file_name in enumerate(file_names):
        curr_data_path = os.path.join(data_path, f"All_mutations_{file_name}.csv")
        df = pd.read_csv(curr_data_path, index_col=0)
        if res_df is None:
            res_df = df[['Mutations number', 'Mut pos', 'full_seq', 'Total count']]
        else:
            res_df = pd.merge(res_df, df[['Total count']], how="left", left_index=True, right_index=True)
        res_df = res_df.rename(columns={'Total count': final_names[index]})

    res_df = res_df.fillna(0)
    seq = res_df['full_seq'].str[3:4] + res_df['full_seq'].str[34:35] + res_df['full_seq'].str[37:38] + \
          res_df['full_seq'].str[67:68] + res_df['full_seq'].str[70:71] + res_df['full_seq'].str[96:97] + \
          res_df['full_seq'].str[98:99]

    res_df.index = seq.values

    wanted_position = [0, 4, 35, 38, 68, 71, 97, 99]
    valid_index = res_df['full_seq'].apply(is_interface_mut, WTseq=WTaa, interface_mut=wanted_position)
    filter_df = res_df[valid_index]
    return pd.concat([filter_df.iloc[:, :3], filter_df.loc[:, relevant_columns]], axis=1)


def generate_label(df,gate):
    """
    Computes the enrichment ratio (ER) for a given gate using log2 transformation.
    log2[(freq var (gate)/freq w.t (gate))/(freq var (pre)/freq w.t ([pre]))
    log2[(𝑀𝑀𝑃 𝑣𝑎𝑟(𝑔𝑎𝑡𝑒)/Total gate)/(𝑀𝑀𝑃 w.t(gate)/Total gate)]/[(𝑀𝑀𝑃 𝑣𝑎𝑟(𝑔𝑎𝑡𝑒)/Total gate)/(𝑀𝑀𝑃 w.t(gate)/Total gate)]]
    log2[[𝑀𝑀𝑃 𝑣𝑎𝑟(𝑔𝑎𝑡𝑒)/𝑀𝑀𝑃 w.t(gate)]/[𝑀𝑀𝑃 𝑣𝑎𝑟(pre)/𝑀𝑀𝑃 w.t(pre)]]

    :param df: DataFrame - contains frequency counts of variants
    :param gate: string - name of the gate column
    :return: DataFrame - computed enrichment ratio values for the specified gate
    """

    all_gates_wt = df.iloc[0, :]
    ER_df = df / all_gates_wt
    ER_df = ER_df.div(ER_df['pre'], axis=0)

    return ER_df[gate]

def prepare_data(df, indices, label_data_log2, weight_column):
    """
    Prepares data for model training by extracting relevant features and encoding sequences.

    :param df: DataFrame - contains sequence data
    :param indices: list - list of sequence identifiers
    :param label_data_log2: DataFrame - log-transformed labels for the sequences
    :param weight_column: string - name of the column containing weights
    :return: tuple - (data, weights, one-hot encoded sequences)
    """

    weights = df.loc[indices, weight_column]
    data = label_data_log2.loc[indices]
    encoded_indices = np.reshape(
        np.array([oneHot(idx) for idx in indices], dtype='int8'),
        (len(indices), 140)
    )
    return data, weights, encoded_indices

def filter_threshold(the_gate, a_df,a_threshold,):
    """
    Filters variants based on a minimum threshold for a given gate.

    :param the_gate: string - gate name to apply thresholding
    :param a_df: DataFrame - contains sequence data
    :param a_threshold: float - threshold value for filtering
    :return: tuple - filtered variant indices and their summed counts
    """

    gate_name = f"{the_gate}_sum"
    rel_index = a_df[the_gate] >= a_threshold
    filter_data = a_df.loc[rel_index].copy()
    filter_data[gate_name] = filter_data['pre'] + filter_data[the_gate]
    sorted_df = filter_data.sort_values(by=gate_name, ascending=False)
    return sorted_df.index, sorted_df[gate_name]

def person_cal(x, y):
    """
    Computes the Pearson correlation coefficient between two arrays.

    :param x: array - first data series
    :param y: array - second data series
    :return: float - Pearson correlation coefficient
    """

    correlation_matrix = np.corrcoef(x, y)
    return correlation_matrix[0, 1]


def model_initiation(k_wargs):
    """
    Initializes a neural network model for sequence analysis.

    :param k_wargs: dictionary - contains hyperparameters:
        - 'dense_size_1': int - size of the first dense layer
        - 'drop_1': float - dropout rate after the first dense layer
        - 'dense_size_2': int - size of the second dense layer
        - 'drop_2': float - dropout rate after the second dense layer
        - 'seed': int - random seed for reproducibility
    :return: Sequential model - compiled Keras model
    """

    from tensorflow.keras.models import Sequential
    from tensorflow.keras.layers import Dense, Dropout, InputLayer

    model = Sequential()
    model.add(InputLayer(input_shape=(140,)))
    model.add(Dense(k_wargs['dense_size_1'], activation='relu'))
    model.add(Dropout(k_wargs['drop_1'], seed=k_wargs['seed']))
    model.add(Dense(k_wargs['dense_size_2'], activation='relu'))
    model.add(Dropout(k_wargs['drop_2'], seed=k_wargs['seed']))
    model.add(Dense(1, activation='linear'))
    return model


def set_random_seeds(seed_value):
    """
    Sets random seeds for TensorFlow, NumPy, and Python's random module to ensure reproducibility.

    :param seed_value: int - seed value for reproducibility
    """

    import random
    import tensorflow as tf
    # Set random seeds for reproducibility
    np.random.seed(seed_value)
    tf.keras.utils.set_random_seed(seed_value)
    random.seed(seed_value)

def split_data_test_val_train(dataframe,size):
    """
    Splits a dataframe into test, validation, and training sets based on a fixed size.

    :param dataframe: DataFrame - dataset to split
    :param size: int - number of samples per split (test and validation)
    :return: tuple - (test set, validation set, training set)
    """

    test = dataframe[:size]
    val = dataframe[size: 2 * size]
    train = dataframe[2 * size:]
    return test, val, train

def shuffle_data(index, ER, wight, seed):
    """
    Shuffles dataset indices, expression ratios (ER), and weights.

    :param index: array - indices of samples
    :param ER: array - expression ratios
    :param wight: array - weights associated with samples
    :param seed: int - seed for randomization
    :return: tuple - shuffled indices, expression ratios, and weights
    """

    # Set seeds for reproducibility
    set_random_seeds(seed)

    # shuffle the train parameters: input, value, and weights
    shuffle_idx = np.random.permutation(len(index))
    shuffle_index = index[shuffle_idx]
    shuffle_value = ER[shuffle_idx]
    shuffle_wight = wight[shuffle_idx]

    return shuffle_index, shuffle_value, shuffle_wight







