
def get_training_choice():
    """
    Prompts the user to select a training methodology and validates the input.

    :return: int - 1 if training on 90% of the data with a test set,
                  2 if training on 100% of the data.
    """
    valid_choices = {"1", "2"}
    while (choice := input(
            "Select a training methodology:\n"
            "1: Train on 90% of the data, using the remaining 10% as a test set.\n"
            "2: Train on 100% of the data without a separate test set.\n"
            "Enter your choice (1 or 2): "
    )) not in valid_choices:
        print("Invalid input. Please enter 1 or 2.")

    return int(choice)

def predict_model_100(train_inputs, train_ER, weight, hyper_parameters, seed):
    """
    Trains a neural network model to predict experimental response (ER) values
    based on all data.

    :param train_inputs: np.ndarray - input features for training
    :param train_ER: np.ndarray - target ER values for training
    :param weight: np.ndarray - sample weights for training (used if 'with_weight' is specified)
    :param hyper_parameters: dict - contains hyperparameters for model training:
        - 'dense_size_1': int - size of the first dense layer
        - 'drop_1': float - dropout rate after the first dense layer
        - 'dense_size_2': int - size of the second dense layer
        - 'drop_2': float - dropout rate after the second dense layer
        - 'lr': float - learning rate for the optimizer
        - 'epochs': int - number of training epochs
        - 'batch': int - batch size for training
        - 'weight': str - whether to use sample weights ('with_weight' or 'no_weight')
    :param seed: int - random seed for reproducibility
    :return: Sequential model - trained Keras model
    """

    from tensorflow.keras.losses import MeanSquaredError

    hyper_parameters['seed'] = seed

    # Shuffle data for training
    train_index, train_value, train_weights = shuffle_data(train_inputs,train_ER, weight,seed)

    # Initialize model and compile
    model = model_initiation(hyper_parameters)
    model.compile(optimizer=tf.optimizers.Adam(learning_rate=hyper_parameters['lr']),
                  loss=MeanSquaredError())

    if hyper_parameters['weight'] == 'with_weight':
        model.fit(train_index, train_value, epochs=hyper_parameters['epochs'],
                  batch_size=hyper_parameters['batch'], verbose=0, shuffle=True,
                  sample_weight=train_weights)

    elif hyper_parameters['weight'] == 'no_weight':
        model.fit(train_index, train_value, epochs=hyper_parameters['epochs'],
                  batch_size=hyper_parameters['batch'], verbose=0, shuffle=True)

    return model

def train_and_evaluate(gate, filter_data, label_data_log2, weight_column, parameters_dic, choice, path, num_models):
    """
    Trains and evaluates models for a specific gate.

    :param gate: str - name of the target gate (e.g., 'MMP-1 G1')
    :param filter_data: pd.DataFrame - filtered dataset for training and testing
    :param label_data_log2: pd.DataFrame - log2-transformed target values
    :param weight_column: str - name of the weight column used for training
    :param parameters_dic: dict - contains hyperparameter configurations for models
    :param choice: int - training choice (1: train-test split, 2: full dataset training)
    :param path: str - directory path to save trained models
    :param num_models: int - number of models to train for ensemble prediction
    :return: None
    """

    if choice == 1: # train 90% test the remaining 10%
        size_10 = int(0.1 * len(filter_data))
        train_indexes = filter_data.index[size_10:]
        test_indexes = filter_data.index[:size_10]

        # Prepare test data
        test_data, test_weight, test_index_encoded = prepare_data(filter_data, test_indexes, label_data_log2,
                                                       weight_column)
    else: # 100% of the data using for training
        train_indexes =  filter_data.index[:]

    # Prepare training data
    train_data, train_weight, train_index_encoded = prepare_data(filter_data, train_indexes, label_data_log2,
                                                                 weight_column)

    total_pred = np.zeros_like(test_data) if choice == 1 else None

    for seed in range(num_models):
        parameters_dic[gate]['seed'] = seed

        # Train the model
        model = predict_model_100(train_index_encoded, train_data, train_weight, parameters_dic[gate], seed)

        # Save the model
        model_dir = os.path.join(path, gate)
        os.makedirs(model_dir, exist_ok=True)  # Ensure directory exists
        model_path = os.path.join(model_dir, f"model_{gate}_{seed}_{now}.h5")
        model.save(model_path)
        print(f"Trained {gate}, iteration {seed + 1}\n"
              f"Model saved at: {model_path}")

        if choice == 1:
            test_pred = model.predict(test_index_encoded).flatten()
            total_pred += test_pred
            print(f"Pearson correlation (iteration {seed + 1}): {person_cal(test_data, test_pred):.4f}")

    if choice == 1:
        total_pred = total_pred / NUM_MODELS
        person = person_cal(test_data, total_pred)

        # Ensure the directory exists before saving
        test_summary_dir = os.path.join(PATH, "10_test_summary")
        os.makedirs(test_summary_dir, exist_ok=True)

        # Save results
        summary_df = pd.DataFrame({'Observed log2 ER': test_data, 'Predicted log2 ER': total_pred})
        summary_path = os.path.join(test_summary_dir, f"{gate}_final_data.csv")
        summary_df.to_csv(summary_path, index=False)

        print(f"Pearson correlation after training {NUM_MODELS} models on {gate}: {person}")

    print(f"finish working on {gate}")


if __name__ == '__main__':
    import time
    import os
    import tensorflow as tf
    import numpy as np
    import pandas as pd
    from datetime import datetime
    from functions import (model_initiation,shuffle_data, create_res_df_filter,
                           generate_label, prepare_data, person_cal)

    # Initiate parameters for training models
    script_start_time = time.time()
    threshold = 1
    NUM_MODELS = 10
    now = datetime.now().strftime("%Y%m%d")
    choice = get_training_choice()

    # Model file path
    if choice == 1:
        file_id = f"{now}_90_best_parameters"
    else:
        file_id = f"{now}_100_best_parameters"

    PATH = f"./Saved_models/model_{file_id}"
    os.makedirs(PATH, exist_ok=True)

    # Load data
    raw_data = "./pre_process_data"
    final_names = ['pre', 'MMP-1 G1', 'MMP-3 G1', 'MMP-1 G4', 'MMP-3 G4']
    final_names_no_pre = final_names[1:]
    res_df_filter = create_res_df_filter(raw_data, final_names)
    print(f"Number of different variants: {res_df_filter.shape[0]}")

    # Optimized hyperparameters for each model
    parameters_dic = {
        'MMP-1 G1': {'dense_size_1': 8, 'drop_1': 0, 'dense_size_2': 2, 'drop_2': 0.3, 'weight': 'no_weight',
                     'batch': 4, 'lr': 0.005, 'epochs': 30},
        'MMP-1 G4': {'dense_size_1': 32, 'drop_1': 0.2, 'dense_size_2': 4, 'drop_2': 0.1, 'weight': 'no_weight',
                     'batch': 8, 'lr': 0.001, 'epochs': 50},
        'MMP-3 G1': {'dense_size_1': 32, 'drop_1': 0.2, 'dense_size_2': 4, 'drop_2': 0.1, 'weight': 'no_weight',
                     'batch': 2, 'lr': 0.001, 'epochs': 40},
        'MMP-3 G4': {'dense_size_1': 8, 'drop_1': 0, 'dense_size_2': 2, 'drop_2': 0.3, 'weight': 'with_weight',
                     'batch': 8, 'lr': 0.005, 'epochs': 40}
    }

    # Iterate over each gate and train models
    for gate in final_names_no_pre:
        df = res_df_filter[['Mutations number', 'pre', gate]]
        weight_column = f"{gate}_sum"
        filter_data = df[df[gate] >= threshold].copy()
        filter_data[weight_column] = np.log2(filter_data['pre'] + filter_data[gate].replace(0, np.nan))

        # Prepare the label data and sort according to size
        label_data = generate_label(df, gate)
        label_data_log2 = np.log2(label_data.replace(0, np.nan))
        filter_data.sort_values(by=weight_column, ascending=False, inplace=True)

        # train and evaluate the model
        train_and_evaluate(gate, filter_data, label_data_log2, weight_column,
                           parameters_dic, choice, PATH, NUM_MODELS)
