import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


data = pd.read_csv('gca_peptides_for_deeplc_95thperc_observations_ms2pip.tsv', sep='\t')

print("Number from DeepLC", len(data))

# Remove rows that the sequence_x length is less than 8 aa
data = data[data['sequence_x'].str.len() >= 8]

# Remove ratio total_ions/number_peaks > 1 the number of total ions should be less than the number of peaks.
data = data[data['total_ions']/data['number_peaks'] < 1]

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
plt.savefig("./ms2pip-plots/signal_to_noise.png")

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

#round to 2 decimal places
threshold_value_b = round(threshold_value_b, 2)
threshold_value_y = round(threshold_value_y, 2)
threshold_value_dp = round(threshold_value_dp, 2)

# Filter entries with scores above the threshold
#good_hits = data_filtered[data_filtered['dot_product'] >= threshold_value]

# Display the threshold value
print(f"B Threshold Value (at {threshold_percentile}th percentile): {threshold_value_b}")
print(f"Y Threshold Value (at {threshold_percentile}th percentile): {threshold_value_y}")
print(f"Dot Product Threshold Value (at {threshold_percentile}th percentile): {threshold_value_dp}")


# Plot the distribution of probability scores
plt.hist(data['corrected_pearsonr_B'], bins=50, alpha=0.6, label='b')
plt.hist(data['corrected_pearsonr_Y'], bins=50, alpha=0.6, label='y')
plt.hist(data['corrected_dot_product'], bins=50, alpha=0.6, label='MS2PIP Adjusted Score')
plt.xlabel('B and Y ion Pearson Distribution, and MS2PIP Adjusted Score')
plt.vlines(threshold_value_dp, 0, 1600, color='red', label=f'{threshold_percentile}th Percentile  MS2PIP Adjusted Score {threshold_value_dp}')
plt.ylabel('Frequency')
plt.title('B and Y ion Pearson Correlations')
plt.legend()
plt.savefig("./ms2pip-plots/ms2pip_filtered_score.png")


data = data[data['corrected_dot_product'] > threshold_value_dp]

## plot the distribution of peptides per number of observations in PeptideAtlas column (PeptideAtlas_observations), the number of observations
# must be binned in 20 bins.

# Remove redundant peptides (same sequence)
data_peptides = data.drop_duplicates(subset=['sequence_x'])
# Define custom bins
bins = [0, 10, 20, 30, 40, 50, np.inf]

# Plotting the histogram
plt.show()
plt.hist(data_peptides['PeptideAtlas_observations'], bins=bins)
plt.title('Distribution of PeptideAtlas observations')
plt.xlabel('PeptideAtlas observations')
plt.ylabel('Frequency')
plt.xticks(bins[:-1])  # Set x-axis ticks to bin edges for clarity
plt.savefig("./peptide-atlas-plots/peptide_atlas.png")

data.to_csv('gca_peptides_for_deeplc_95thperc_observations_ms2pip_filtered.tsv', sep='\t', index=False)

print("Number of rows with PeptideAtlas observations = 0:", len(data_peptides[data_peptides['PeptideAtlas_observations'] == 0]))
print("Number of rows with PeptideAtlas observations > 0:", len(data_peptides[data_peptides['PeptideAtlas_observations'] > 0]))

print("Number of rows after filtering MS2PIP:", len(data))
print("Number of unique sequences in the results:", len(data_peptides))

# Sort the protein accessions divided by comma in descending order in the column protein_accessions.

data_peptides['protein_accessions'] = data_peptides['protein_accessions'].str.split(',').apply(lambda x: sorted(x, reverse=True))

# Count the number of times a protein accession appears in the dataframe.
proteins_count = data_peptides['protein_accessions'].value_counts()
proteins_count = proteins_count.reset_index()
proteins_count.columns = ['protein_accessions', 'count']
proteins_count['number_accessions_group'] = proteins_count['protein_accessions'].str.len()
proteins_count.to_csv('proteins_count.tsv', sep = "\t", index=False)
print(proteins_count)


# Display the entries with good hits
print("Entries with Good Hits:")
#print(good_hits)
