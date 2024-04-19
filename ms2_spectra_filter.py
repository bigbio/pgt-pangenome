import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import click

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])


@click.group(context_settings=CONTEXT_SETTINGS)
def cli():
    """
    This is the main tool that gives access to all commands to convert SDRF files into pipeline-specific configuration files
    """
    pass


@click.command("filter-ms2pip", help="Run the ms2pip filtering process to remove low-quality peptides.")
@click.option( "-p", "--peptide_file", type=str, required=True, help="Peptide sequence to be used for the ms2pip")
@click.option("-f", "--folder_output", type=str, required=True, help="Output file after filtering ms2pip threshold")
@click.option("-o", "--output_file", type=str, required=True, help="Output file after filtering ms2pip threshold")
@click.option("--number_aa", type=int, default=8, help="Minimum number of amino acids in the peptide sequence")
def filter_ms2pip(peptide_file: str, folder_output: str, output_file: str, number_aa: int = 8):

    data = pd.read_csv(peptide_file, sep=',')
    print("Number from Previous Step and before filtering:", len(data))

    # Remove rows that the sequence_x length is less than 8 aa
    data = data[data['seq'].str.len() >= number_aa]

    # Remove ratio total_ions/number_peaks > 1 the number of total ions should be less than the number of peaks.
    data = data[data['total_ions'] / data['number_peaks'] < 1]

    # Dynamically set thresholds based on percentiles
    lower_percentile = 5
    upper_percentile = 95

    lower_threshold = data['signal_to_noise'].quantile(lower_percentile / 100)
    upper_threshold = data['signal_to_noise'].quantile(upper_percentile / 100)
    lower_threshold = round(lower_threshold, 2)
    upper_threshold = round(upper_threshold, 2)

    print(f"Lower Threshold for SNR: {lower_threshold}")
    print(f"Upper Threshold for SNR: {upper_threshold}")

    plt.hist(data['signal_to_noise'], bins=500)
    plt.title('Signal to Noise (SNR) Distribution')
    plt.vlines(lower_threshold, 0, 100, color='red', label=f'Lower 5th Percentile threshold {lower_threshold}')
    plt.vlines(upper_threshold, 0, 100, color='red', label=f'Upper 95th Percentile threshold {upper_threshold}')
    plt.xlabel('SNR')
    plt.ylabel('Frequency')
    plt.legend()
    plt.savefig(folder_output + "/signal_to_noise.png")

    # Remove rows with y-ions and b-ions wrongly correlated
    data = data[data['pearsonr_Y'] > 0.1]
    data = data[data['signal_to_noise'] > lower_threshold]
    data = data[data['signal_to_noise'] < upper_threshold]

    print("Number of rows after filtering SNR:", len(data))

    # Plot the distribution of probability scores
    plt.hist(data['corrected_pearsonr_B'], bins=50, alpha=0.6, label='b')
    plt.hist(data['corrected_pearsonr_Y'], bins=50, alpha=0.6, label='y')
    plt.xlabel('B and Y ion Pearson Correlations')
    plt.ylabel('Frequency')
    plt.title('B and Y ion Pearson Correlations')
    plt.show()

    # Choose a percentile threshold (e.g., 10th percentile)
    threshold_percentile = 5

    data_filtered_b = data
    data_filtered_y = data
    data_filtered_dp = data
    threshold_value_b = np.percentile(data_filtered_b['corrected_pearsonr_B'], threshold_percentile)
    threshold_value_y = np.percentile(data_filtered_b['corrected_pearsonr_Y'], threshold_percentile)
    threshold_value_dp = np.percentile(data_filtered_dp['corrected_dot_product'], threshold_percentile)

    # round to 2 decimal places
    threshold_value_b = round(threshold_value_b, 2)
    threshold_value_y = round(threshold_value_y, 2)
    threshold_value_dp = round(threshold_value_dp, 2)

    # Display the threshold value
    print(f"B Threshold Value (at {threshold_percentile}th percentile): {threshold_value_b}")
    print(f"Y Threshold Value (at {threshold_percentile}th percentile): {threshold_value_y}")
    print(f"Dot Product Threshold Value (at {threshold_percentile}th percentile): {threshold_value_dp}")

    # Plot the distribution of probability scores
    plt.hist(data['corrected_pearsonr_B'], bins=50, alpha=0.6, label='b')
    plt.hist(data['corrected_pearsonr_Y'], bins=50, alpha=0.6, label='y')
    plt.hist(data['corrected_dot_product'], bins=50, alpha=0.6, label='MS2PIP Adjusted Score')
    plt.xlabel('B and Y ion Pearson Distribution, and MS2PIP Adjusted Score')
    plt.vlines(threshold_value_dp, 0, 1600, color='red',
               label=f'{threshold_percentile}th Percentile  MS2PIP Adjusted Score {threshold_value_dp}')
    plt.ylabel('Frequency')
    plt.title('B and Y ion Pearson Correlations')
    plt.legend()
    plt.savefig(folder_output + "/ms2pip_filtered_score.png")

    data = data[data['corrected_dot_product'] > threshold_value_dp]
    data.to_csv(folder_output + '/' + output_file, sep=',', index=False)


cli.add_command(filter_ms2pip)

if __name__ == "__main__":
    cli()