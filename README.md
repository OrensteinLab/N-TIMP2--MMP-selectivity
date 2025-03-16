# Predicting the affinity landscape of N-TIMP2/MMP1<sub>CAT</sub> or MMP3<sub>CAT</sub> by combining deep neural networks and deep mutational scans

# Introduction

We present a novel approach using deep neural networks trained on High-Throughput Sequencing (HTS) data for protein mutagenesis library affinity screens.
In this study, we have focused on the experimental raw data from the N-terminal domain of the tissue inhibitor of metalloproteinases 2 (N-TIMP2) with the catalytic domain of matrix metalloproteinase 1 and 3 (MMP1<sub>CAT</sub> and MMP3<sub>CAT</sub>). Our goal is to comprehensively measure Protein-Protein Interactions (PPIs) and accurately predict unobserved affinity-enhancing or affinity-reducing variants to isolate selective N-TIMP2 variant.

# Purpose

This repository contains the code related to our research project, which aims to isolate a variant with:
   	- High affinity for N-TIMP2/MMP9<sub>CAT</sub>.
    	- Medium affinity for N-TIMP2/MMP3<sub>CAT</sub>.
     	- The lowest affinity for N-TIMP2/MMP1<sub>CAT</sub>.
This is achieved by combining deep neural networks with deep mutational scans.

# Data
The experimental raw data used in this study consists of N-TIMP2 variants in complex with MMP9<sub>CAT</sub>, MMP3<sub>CAT</sub>, and MMP1<sub>CAT</sub>.
All variants in this dataset contain mutations only at the following seven positions: 4, 35, 38, 68, 71, 97, and 99 (referred to as relevant positions).
HTS data was used to train deep neural networks, allowing accurate prediction of unobserved affinity-enhancing or affinity-reducing variants.

# Setup Environment Instructions (Windows)

Before you proceed with the setup, make sure to have Python and Anaconda installed on your system.


1. **Download the Code Repository:**
   - Visit the GitHub repository: [https://github.com/OrensteinLab/N-TIMP2--MMP-selectivity](https://github.com/OrensteinLab/N-TIMP2--MMP-selectivity)
   - Download the contents of the "Code" folder.

2. **Inside the "Code" Folder, Add Your Raw Data in folder named "Raw_data":**

   For code validation, include the following datasets:
   	- All_mutations_Pre.csv
   	- All_mutations_MMP1_H.csv
   	- All_mutations_MMP1_L.csv
   	- All_mutations_MMP3_H.csv
   	- All_mutations_MMP3_L.csv


3. **Create a Virtual Conda Environment:**
   - Open a command prompt.
   - Navigate to the directory where you downloaded the "Code" repository.
   - Run the following command to create a virtual conda environment named "my_env" with Python 3.8.3 and the required modules:
     
     ```bash
     conda create --name my_env python=3.8.3 --file Requirements.txt 
     ```

4. **Activate the New Environment and Run the Script:**
   - Activate the environment using the following command:
     ```bash
     conda activate my_env
     pip install tensorflow==2.11.0
     ```
   - Run the scripts according to the provided usage instructions.


# Usage
### 1. MMP9 model:
The script `MMP9_prediction.py` takes the `MMP9_model.h5` located in the code folder. To generate N-TIM2 list of variants with high affinity to MMP9.

### 2.	Pre-processing Script (Optional):
The script variants_count.py processes FASTA files to count valid variants.
A valid variant must:  
* Align with the N-TIMP2 wild-type amino acid sequence.
* Contain mutations only at the relevant positions.
 
The script merges data from four repetitions of each library.

Important Note: This script is for users who wish to train a model based on their own data. To work with our data and prdeict N-TIMP2 variants, you dont need to run this script. Becasue the output files of this script (All_mutations_X.csv, where X represents the library subpopulation) are already provided in the data folder.


### 3.	Train models script:
Run `train_save_model.py`. You will be prompted to choose an action:
1. Split the dataset into training (90%), and test (10%) sets. The model is trained on the training set, to predict the test set.
2. Train the model using all available data without making predictions.
   
Output:
* Trained models will be saved in sub folder in "Saved_models".
* If Option 1 is selected, prediction files will be generated at the end of execution.

### 4.	Get predictions:
After training the model by running: `python train_save_model.py` using option number 2, you can now generate predictions.

There are three available scripts for making predictions.
Each script will prompt you to enter the model date in YYYYMMDD format.

#### 4.1	Predicting All possibole Mutations:
The script `predict_all_library.py` generates predictions for all possible mutations (library size of 20^7).

#### 4.2	User-Specified Variant Predictions:
The script `predict_user_variants.py` predicts the log2 ER of any variant.
To use this script, create a TXT file containing full-length protein sequences (99 amino acids), ensuring mutations are only at the relevant positions.
Enter the sequnce in FASTA format:
```bash
>variant_ID
CSCSPVHPQQAFCNADVVIRAKAVSEKEVDSGNDIYGNPIKRIQYEIKQIKMFKGPEKDIEFIYTAPSSAVCGVSLDVGGKKEYLIAGKAEGDGKMHIT
```
Output:
* The script creates a csv file '{current data}_user_predictions', with the Predicted_log2_ER of each variant in the MMP1 and MMP3 models.

#### 4.3	Variant quantification:
The script `calcaulte_library_above_variants.py` calculates how many variants from the entire library have a higher ER than the five benchmark variants presented in the paper:
```bash
['RTDWWIQ', 'RMDWWIQ', 'RTDWWID', 'RFDWWIQ', 'RFDWWID']
```

Output:
* The data is shown as a percentage of the entire library.
* For MMP1 and MMP3, the script creates a csv files '{current data}_df_sum', summarizing the percentage of variants above the benchmark.
   
During execution, you will be prompted to enter the date of the prediction folder in YYYYMMDD format.


