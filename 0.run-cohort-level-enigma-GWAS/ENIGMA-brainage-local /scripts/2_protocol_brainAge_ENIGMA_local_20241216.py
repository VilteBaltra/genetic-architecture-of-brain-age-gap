import os
import pandas as pd
import numpy as np
from photonai.base import Hyperpipe

def process_sex(sex, modelfolder, rawfolder, predictionfolder):
    """
    Process data for a specific sex by loading the model, making predictions, and saving output.

    Parameters:
    sex (str): 'females' or 'males'
    modelfolder (str): Path to the folder containing photon models
    rawfolder (str): Path to the folder containing raw data
    predictionfolder (str): Path to save the predictions
    """
    print(f"Processing {sex} data...")

    # Define file paths
    best_model_file = os.path.join(modelfolder, f'ENIGMA_MDD_{sex}.photon')
    raw_data_file = os.path.join(rawfolder, f'{sex}_raw.csv')
    output_file = os.path.join(predictionfolder, f'{sex}_raw_out.csv')

    # Load the model
    my_model = Hyperpipe.load_optimum_pipe(best_model_file)

    # Load and prepare data
    df = pd.read_csv(raw_data_file, header=0, delimiter=',')
    y_te = df['AGE'].values
    X_te = df.drop(columns=['AGE', 'SUBJID']).values

    # Make predictions
    preds_ind = my_model.predict(X_te)

    # Calculate and print Mean Absolute Error (MAE)
    mae = np.mean(np.abs(y_te - preds_ind))
    print(f'\n\nTest: MAE for {sex} = {mae}')

    # Calculate and print Pearson correlation coefficient (r)
    corr_coef = np.corrcoef(y_te, preds_ind)[0, 1]
    print(f'Test: Pearson correlation coefficient for {sex} = {corr_coef}')

    # Add predictions to the dataframe
    #df = df.drop(columns=['AGE'])  # Remove AGE column
    df['age_prediction'] = preds_ind  # Add predictions as a new column

    # Save the updated dataframe
    df.to_csv(output_file, sep="\t", index=False)
    print(f"Predictions saved to {output_file}\n")

if __name__ == "__main__":
    import argparse

    # Argument parser for user-defined paths
    parser = argparse.ArgumentParser(description="Run age prediction models for males and females.")
    parser.add_argument('--modelfolder', type=str, required=True, help="Path to the folder containing photon models.")
    parser.add_argument('--rawfolder', type=str, required=True, help="Path to the folder containing raw data.")
    parser.add_argument('--predictionfolder', type=str, required=True, help="Path to save the predictions.")

    args = parser.parse_args()

    # Process data for females and males
    for sex in ['females', 'males']:
        process_sex(sex, args.modelfolder, args.rawfolder, args.predictionfolder)

