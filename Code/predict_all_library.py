
def generate_all_combinations(options,pos_num, pre, last):
    """
    Generates all possible sequences of a given length using provided characters,
    appending fixed prefix and suffix, and converts them to one-hot encoding.

    :param options: iterable - characters used for sequence generation
    :param pos_num: int - number of variable positions in the sequence
    :param pre: str - fixed prefix appended to each sequence
    :param last: str - fixed suffix appended to each sequence
    :return: np.ndarray - one-hot encoded sequences of shape (num_combinations, 140)
    """

    combinations = itertools.product(options, repeat=pos_num)
    x = np.array([''.join(comb) + pre + last for comb in combinations])
    all_pos = pd.Series(x)
    return np.reshape(np.array(list(all_pos.apply(oneHot)), dtype='int8'), (len(all_pos), 140))


def load_models_from_paths(models_path, model_groups):
    """
    Loads trained models from specified subdirectories under the given path.

    :param models_path: str - base directory containing model subdirectories
    :param model_groups: list of str - names of subdirectories containing models
    :return: list of list - nested list where each inner list contains loaded models
    """

    models = []
    for group in model_groups:
        group_path = os.path.join(models_path, group)
        curr_models = [load_model(os.path.join(group_path, model)) for model in os.listdir(group_path)]
        models.append(curr_models)

    return models

if __name__ == '__main__':
    import numpy as np
    import time
    import pandas as pd
    import itertools
    import os
    from datetime import datetime
    from tensorflow.keras.models import load_model
    from functions import oneHot

    startTime = time.time()
    aaList = 'ACDEFGHIKLMNPQRSTVWY'
    date = datetime.now().strftime('%Y%m%d')

    folder_name = input("What is the date of the model?\n"
                        "format:YYYYMMDD (in model_YYYYMMDD_100_best_parameters)\n")
    name = f"{date}_100_best_model"
    MODELS_PATH = os.path.join(f"./Saved_models/model_{folder_name}_100_best_parameters")
    directory_path = f"./Predict_all_{name}"
    os.makedirs(directory_path, exist_ok=True)

    model_groups = ['MMP-1 G1', 'MMP-3 G1']
    models = load_models_from_paths(MODELS_PATH, model_groups)
    NUM_MODELS = len(models[0])
    print(f"Initiations: {time.time() - startTime:.2f} sec")

    startTime = time.time()
    for lastAA in aaList:
        for prelastAA in aaList:
            # Generate 20^5 variants
            MMP1_G1_predict = None
            MMP3_G1_predict = None
            num_positions = 5

            # One hot encode the data
            seq_one_hot = generate_all_combinations(aaList, num_positions, prelastAA, lastAA)

            for i in range(NUM_MODELS):
                pred_MMP1 = models[0][i].predict(seq_one_hot,batch_size=32768)
                pred_MMP3 = models[1][i].predict(seq_one_hot,batch_size=32768)

                if MMP1_G1_predict is None:
                    MMP1_G1_predict = pred_MMP1
                    MMP3_G1_predict = pred_MMP3

                else:
                    MMP1_G1_predict = MMP1_G1_predict + pred_MMP1
                    MMP3_G1_predict = MMP3_G1_predict + pred_MMP3

            # calculate the average value
            MMP1_G1_predict = MMP1_G1_predict / NUM_MODELS
            MMP3_G1_predict = MMP3_G1_predict / NUM_MODELS

            file_path_MMP_1 = os.path.join(directory_path, f"{name}_MMP1_G1_{prelastAA}{lastAA}.npy")
            file_path_MMP_3 = os.path.join(directory_path, f"{name}_MMP3_G1_{prelastAA}{lastAA}.npy")

            np.save(file_path_MMP_1, MMP1_G1_predict)
            np.save(file_path_MMP_3, MMP3_G1_predict)
            print(f"prelastAA {prelastAA}, time: {(time.time() - startTime) // 60} min")

        print(f"End cycle: lastAA {lastAA}, time: {(time.time() - startTime) // 60} min")
