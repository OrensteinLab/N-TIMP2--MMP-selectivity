
def raw_data_to_mutations_dict(data_path, all_gates, mut_dict, info):
    """
    The function receives a dictionary of the desired gates and mutations number and returns the dictionary
    full of found sequences

    :param data_path: string of data path
    :param all_gates: dictionary containing all the gates
    :param mut_dict: dictionary containing all the desired number of mutations
    :param info: list containing all the desired information:
                start_motif: string- after the motif the desired sequence will begin
                start_seq_pos: string- after this length of nucleotide the sequence starts
                nuc_length: int- the length of the desired sequence
                WTaa: string- The amino acid sequence of the WT
    :return: mut_dict: dictionary containing all the desired sequences
    """

    for name in all_gates.keys():
        seq_dict = fastq_reader(data_path + name + '.fastq')
        all_gates[name] = sort_seqs_to_mut(seq_dict, info)
        # Creating a data frame of the reads count according to the number of mutations
        for num_of_mut in mut_dict.keys():
            mut_dict[num_of_mut] = pd.concat([mut_dict[num_of_mut],
                                              sort_mut_by_number(all_gates[name], num_of_mut, name)], axis=1, sort=True)
            mut_dict[num_of_mut] = mut_dict[num_of_mut].rename(columns= {'count':name})


def mutations_dict_to_df(mut_dict):
    """
        The function receives a dictionary with valid sequences and returns a data frame containing all variants

        :param mut_dict: dictionary containing all the desired sequences
        :param numerator_gate: str of the desired numerator gate
        :param denominatorr_gate: str of the desired denominator gate
        :return: mut_total_data: dictionary containing all variants information
    """

    mut_total_data = pd.DataFrame()
    for num_of_mut in mut_dict.keys():
        mut_dict[num_of_mut]['mutations number'] = num_of_mut
        mut_total = pd.concat([mut_dict[num_of_mut].loc[:, 'Mutations number'],
                               mut_dictionary[num_of_mut]['full_seq'],
                               mut_dict[num_of_mut].loc[:, 'Mutation'],
                               mut_dict[num_of_mut].loc[:, 'Mut pos'],
                               mut_dict[num_of_mut].loc[:, file_name + '_1'].fillna(0),
                               mut_dict[num_of_mut].loc[:, file_name + '_2'].fillna(0),
                              mut_dict[num_of_mut].loc[:, file_name + '_3'].fillna(0),
                              mut_dict[num_of_mut].loc[:, file_name + '_4'].fillna(0)],
                              axis=1, sort=True)

        mut_total_data = pd.concat([mut_total_data, mut_total], sort=False)

    # Count the total number
    file_columns = [file_name + f'_{i}' for i in range(1, 5)]
    mut_total_data['Total count'] = mut_total_data[file_columns].sum(axis=1).fillna(0)

    # Sort the dataframe according to mutation number, and total count
    mut_total_data = mut_total_data.sort_values(by=['Mutations number','Total count'],
                                                ascending=[True,False])

    return mut_total_data

def filter_valid_variants(df,WTAA):
    """
        The function receives a dictionary with valid sequences and returns a data frame containing all variants

        :param mut_dict: dictionary containing all the desired sequences
        :param numerator_gate: str of the desired numerator gate
        :param denominatorr_gate: str of the desired denominator gate
        :return: mut_total_data: dictionary containing all variants information
    """
    seq = df['full_seq'].str[3:4] + df['full_seq'].str[34:35] + df['full_seq'].str[37:38] + \
          df['full_seq'].str[67:68] + df['full_seq'].str[70:71] + df['full_seq'].str[96:97] + \
          df['full_seq'].str[98:99]
    df.index = seq.values
    wanted_position = [0, 4, 35, 38, 68, 71, 97, 99]
    valid_index = df['full_seq'].apply(is_interface_mut, WTseq=WTAA, interface_mut=wanted_position)
    return df[valid_index]


if __name__ == '__main__':
    import os
    import pandas as pd
    from functions import fastq_reader, sort_seqs_to_mut, sort_mut_by_number, compare_seq, is_interface_mut

    # TIMP library NGS- variables
    WTaa = 'CSCSPVHPQQAFCNADVVIRAKAVSEKEVDSGNDIYGNPIKRIQYEIKQIKMFKGPEKDIEFIYTAPSSAVCGVSLDVGGKKEYLIAGKAEGDGKMHIT'
    WTnuc = 'TGCAGCTGCTCCCCGGTGCACCCGCAACAGGCGTTTTGCAATGCAGATGTAGTGATCAGGGCCAAAGCGGTCAGTGAGAAGGAAGTGGACTCTGGAAAC' \
            'GACATCTATGGCAACCCTATCAAGAGGATCCAGTATGAGATCAAGCAGATAAAGATGTTCAAAGGGCCTGAGAAGGATATAGAGTTTATCTACACGGCC' \
            'CCCTCCTCGGCAGTGTGTGGGGTCTCGCTGGACGTTGGAGGAAAGAAGGAATATCTCATTGCAGGAAAGGCCGAGGGGGACGGCAAGATGCACATCACC'
    start_motif = 'AGAGA'
    nuc_length = len(WTnuc)


    data_path = os.path.join("./raw_data_names/")
    file_path = os.path.join("./pre_process_data_new")
    os.makedirs(file_path, exist_ok=True)

    file_names = ['Pre','MMP1_L','MMP1_H','MMP3_L','MMP3_H']

    for file_name in file_names:
        # Create dictionary for all gates
        all_gates_dict = {}
        for num in range(1,5):
            all_gates_dict[file_name + '_' + str(num)] = {}

        mut_dictionary = {0: pd.DataFrame(), 1: pd.DataFrame(), 2: pd.DataFrame(), 3: pd.DataFrame(), 4: pd.DataFrame(),
                          5: pd.DataFrame(), 6: pd.DataFrame(), 7: pd.DataFrame(), 8: pd.DataFrame()}  # 8-> 8 or bigger

        # Read the whole data to the all_gates_dict, and count each mutation in the mut_dictionary.
        # In mut_dictionary each column represent different fasta file.
        raw_data_to_mutations_dict(data_path, all_gates_dict, mut_dictionary, [start_motif, len(start_motif), nuc_length, WTaa])

        # Upload the Mutations number, Mut pos and full_seq of each variant.
        for num_of_mut in mut_dictionary.keys():
            mut_dictionary[num_of_mut]['Mutations number'] = num_of_mut
            mut_dictionary[num_of_mut]['full_seq'] = mut_dictionary[num_of_mut].index
            mut_dictionary[num_of_mut]['Mutation'] = mut_dictionary[num_of_mut]['full_seq'].apply(
                lambda x: compare_seq(seq1=WTaa, seq2=x)[1])

            if num_of_mut == 0:  # Index 2 contain the mutation position.
                mut_dictionary[num_of_mut]['Mut pos'] = [[0]]
            else:
                mut_dictionary[num_of_mut]['Mut pos'] = mut_dictionary[num_of_mut]['full_seq'].apply(
                    lambda x: compare_seq(seq1=WTaa, seq2=x)[2])

        raw_data = mutations_dict_to_df(mut_dictionary)

        filter_df = filter_valid_variants(raw_data,WTaa)


        filter_df.to_csv(f'{file_path}/All_mutations_{file_name}.csv')
        print(f"finish {file_name}")


