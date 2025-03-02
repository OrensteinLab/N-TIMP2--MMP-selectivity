import pandas as pd
import os
from datetime import datetime
from functions import is_interface_mut
from calcaulte_library_above_variants import load_models_from_paths, generate_prediction


def read_fasta_to_df(fasta_file):
    """
        Reads a FASTA file and creates a DataFrame with sequence IDs and amino acid sequences.

        :param fasta_file: str - Path to the FASTA file.
        :return: pd.DataFrame - DataFrame with columns ['ID', 'full_seq'].
    """
    if not os.path.exists(fasta_file):
        raise FileNotFoundError(f"Error: {fasta_file} not found.")

    sequences = []
    with open(fasta_file, "r") as file:
        seq_id, sequence = None, []
        for line in file:
            line = line.strip()
            if line.startswith(">"):
                if seq_id:
                    sequences.append([seq_id, "".join(sequence)])
                seq_id, sequence = line[1:], []
            else:
                sequence.append(line)
        if seq_id:
            sequences.append([seq_id, "".join(sequence)])

    return pd.DataFrame(sequences, columns=["ID", "full_seq"])


def filter_and_extract_variants(df, positions, WTaa):
    """
    Filters sequences that match the reference length and contain mutations only in valid positions.

    :param df: pd.DataFrame - DataFrame containing sequences.
    :param valid_positions: list - List of valid mutation positions (0-based index).
    :param WTaa: str - Reference sequence.
    :return: pd.DataFrame - Filtered DataFrame with extracted mutation sequences.
    """
    # filter variants not in the size of the protein
    df = df[df['full_seq'].str.len() == len(WTaa)]

    # filter variants that have mutaions not in the valid positions
    df = df[df['full_seq'].apply(lambda seq: is_interface_mut(seq, WTseq=WTaa, interface_mut=positions))]
    df['seq'] = df['full_seq'].apply(lambda x: ''.join([x[i - 1] for i in positions]))
    return df


if __name__ == "__main__":
    # Parameters
    WTaa = 'CSCSPVHPQQAFCNADVVIRAKAVSEKEVDSGNDIYGNPIKRIQYEIKQIKMFKGPEKDIEFIYTAPSSAVCGVSLDVGGKKEYLIAGKAEGDGKMHIT'
    valid_positions = [4, 35, 38, 68, 71, 97, 99]
    date = datetime.now().strftime('%Y%m%d')
    proteins = ['MMP1', 'MMP3']

    # Read and filter FASTA file
    fasta_path = "./file.txt"
    df = read_fasta_to_df(fasta_path)
    df = filter_and_extract_variants(df, valid_positions, WTaa)

    print(f"Found {len(df)} valid sequences:")
    print(df[['ID', 'seq']])

    # Get model path from user
    model_name = input("What is the date of the model?\n"
                        "format:YYYYMMDD (in model_YYYYMMDD_100_best_parameters)\n")
    models_path = os.path.join("./Saved_models", f"model_{model_name}_100_best_parameters")

    if not os.path.exists(models_path):
        raise FileNotFoundError(f"Error: Model path {models_path} does not exist.")

    # Prepare results directory
    result_path = "./User_Predictions"
    os.makedirs(result_path, exist_ok=True)

    # Load models and predict
    model_groups = ['MMP-1 G1', 'MMP-3 G1']
    models = load_models_from_paths(models_path, model_groups)
    mut_predictions = generate_prediction(models, df['seq'].values)

    # Add predictions to DataFrame
    for i, protein in enumerate(proteins):
        df[f"Predicted_log2_ER_to_{protein}"] = [round(val, 3) for val in mut_predictions[i]]


    # Save results
    output_file = os.path.join(result_path, f"{date}_user_predictions.csv")
    df[['ID', 'seq'] + [f"Predicted_log2_ER_to_{p}" for p in proteins]].to_csv(output_file, index=False)

    print(f"Predictions saved in {os.path.abspath(output_file)}")
