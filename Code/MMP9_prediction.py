
def multi_mut_np(variants, WTaa):
    """Vectorized function to check if all relevant positions are mutated."""
    return np.all(variants[:, :7] != WTaa[:7], axis=1).astype(np.uint8)  # Faster than apply()

if __name__ == '__main__':
    import numpy as np
    import pandas as pd
    import tensorflow as tf
    import time
    from itertools import product
    from functions import oneHot  # Ensure oneHot is optimized for batch processing

    # Load the model without the optimizer
    load_model = tf.keras.models.load_model('original_model.h5', compile=False)

    # Define amino acids and number of mutating positions
    aaList = np.array(list("ACDEFGHIKLMNPQRSTVWY"))  # Convert to NumPy array for efficiency
    num_positions = 5  # Number of mutation sites
    WTaa = np.array(list("SINSVHT"))  # Convert to NumPy array for fast comparisons


    # Initialize storage for results
    sum_df_5 = None
    sum_df_2 = None

    start_time = time.time()

    # Loop through lastAA and prelastAA combinations
    for lastAA in aaList:
        for prelastAA in aaList:
            cycle_time = time.time()  # Track time per cycle

            # **Step 1: Efficiently generate all possible sequences**
            combinations = np.array([''.join(comb) for comb in product(aaList, repeat=num_positions)])
            all_variants = np.char.add(combinations, prelastAA + lastAA)  # Append last two residues

            # **Step 2: Convert sequences to one-hot encoding (fully vectorized)**
            seq_one_hot = np.array([oneHot(seq) for seq in all_variants],
                                   dtype=np.uint8).reshape(len(all_variants), 140)
            seq_one_hot = np.hstack(
                [seq_one_hot, np.zeros((seq_one_hot.shape[0], 1), dtype=np.uint8)])  # Add extra column

            # **Step 3: Predict using the model**
            pred_test_reg = load_model.predict(seq_one_hot, batch_size=32768)

            # **Step 4: Store predictions efficiently**
            reg_dataset = pd.DataFrame({
                'variant': all_variants,
                'pred_total': pred_test_reg.flatten()
            })

            # **Step 5: Use NumPy boolean indexing instead of apply()**
            variant_matrix = np.array([list(seq) for seq in all_variants])  # Convert to char matrix for fast comparisons
            is_multi_mut = multi_mut_np(variant_matrix, WTaa)  # Vectorized mutation check

            # **Step 6: Select top predictions**
            filtered_data = reg_dataset[is_multi_mut == 1]
            rslt_df5 = filtered_data.nlargest(5, 'pred_total')

            # **Step 7: Append results efficiently**
            sum_df_5 = rslt_df5 if sum_df_5 is None else pd.concat([sum_df_5, rslt_df5], axis=0)

            print(f"prelastAA {prelastAA}, time: {(time.time() - cycle_time):.2f} sec")

            sum_df_5.to_csv('sum_df_5.csv', index=False)

        print(f"End cycle: lastAA {lastAA}, total time: {(time.time() - start_time) // 60} min")

    # Sort sum_df_2 and sum_df_5 by 'pred_total' in descending order before saving
    sum_df_5 = sum_df_5.sort_values(by="pred_total", ascending=False)

    # Save the sorted DataFrames
    sum_df_5.to_csv('MMP9_highest_ER.csv', index=False)

