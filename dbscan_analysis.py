import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt
import seaborn as sns


# Assuming 'data' is your DataFrame with columns: PEP, Retention Time, Observations
# Replace 'data.csv' with your actual dataset
cols = ["Posterior error probability", "abserror", "error_percentile", "PeptideAtlas_observations", "GPMDB_observations"]
data = pd.read_csv('gca_peptides_for_deeplc_95thperc_observations.tsv', usecols=cols, sep='\t')
# Feature Engineering: Standardize the numerical features
numerical_features = ['Posterior error probability', 'error_percentile', 'PeptideAtlas_observations']

scaler = StandardScaler()
data_scaled = scaler.fit_transform(data[numerical_features])

# Perform KMeans clustering for each individual feature
for feature in numerical_features:
    feature_data = data[[feature]]

    # Create a DataFrame with the scaled feature and cluster labels
    kmeans = KMeans(n_clusters=3)
    feature_df = pd.DataFrame(data={feature: feature_data[feature], 'Cluster': kmeans.fit_predict(feature_data)})

    # Count the number of points in each cluster
    cluster_sizes = feature_df['Cluster'].value_counts().sort_index()

    # Visualize the clusters
    plt.figure(figsize=(10, 6))
    sns.scatterplot(x=feature, y=feature, hue='Cluster', data=feature_df, palette='viridis')

    # Display cluster sizes
    for cluster, size in cluster_sizes.items():
        plt.text(feature_df[feature][feature_df['Cluster'] == cluster].mean(),
                 feature_df[feature][feature_df['Cluster'] == cluster].mean(),
                 f'Cluster {cluster} (Size: {size})', fontsize=10, color='black')

    plt.title(f'Clustering of Peptides based on {feature}')
    plt.show()