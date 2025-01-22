#!/usr/bin/python

import pandas as pd
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import squareform, pdist
import matplotlib.pyplot as plt
from sklearn.metrics import silhouette_score
import os
import argparse

parser = argparse.ArgumentParser(description="Adding clusters to samples file")
parser.add_argument("-k", type=int, required=True, help="number of clusters desired by the user")
parser.add_argument("-f", type=str, required=True, help="file path to the table with the distance matrix")
args = parser.parse_args()

file_path = args.f
data = pd.read_csv(file_path, sep='\t', header=0)

print(os.getcwd())

labels = data.iloc[:, 0].values
print(labels)
values = data.iloc[:, 1:].values

# Removing the whole path from the labels that will be shown in the dendrogram, for readability
adj_labels = []

for i in labels:
    adj_labels.append(os.path.basename(i))

print(adj_labels)

condensed_distance_matrix = squareform(values)
Z = linkage(condensed_distance_matrix, method='ward')

# Dendrogram
plt.figure(figsize=(10, 7))
dendrogram(Z, labels=adj_labels)
plt.title('Dendrogram')
plt.xlabel('Samples')
plt.ylabel('Distance')
plt.xticks(rotation=90)
plt.savefig('results/clustering/dendrogram.png', dpi=300, bbox_inches='tight')

## DEFINED BY THE USER
# Clustering based on desired number of clusters (DEFINED BY THE USER)
min_k = 2
max_k = [len(adj_labels) / 2]
k = args.k

if k == 0:
    print("Ideal number of clusters being calculated...")
    if len(adj_labels) < 3:
        print("The number of samples is too low to perform clustering.")
        best_k = 1
    # To check what the ideal K is with the Silhouette Method
    # Determine the silhouette score for different values of k
    else:
        max_clusters = len(adj_labels) - 1  # Define the maximum number of clusters to test
        best_k = None  # Initialize the best number of clusters
        best_score = float('-inf')  # Initialize the best silhouette score

        silhouette_scores = []
        for K in range(2, max_clusters + 1):  # Silhouette score is not defined for k=1
            cluster_labels = fcluster(Z, K, criterion='maxclust')
            print(cluster_labels)
            print(f"Testing K={K}...")
            if len(set(cluster_labels)) < 2:
                print(f"Only one cluster was found for K={K}. Skipping silhouette score calculation")
                continue
            score = silhouette_score(squareform(condensed_distance_matrix), cluster_labels, metric='precomputed')
            silhouette_scores.append(score)
            if score > best_score:
                best_score = score
                best_k = K
        if best_k is not None:
            print(f"The optimal number of clusters is {best_k} with a silhouette score of {best_score}")
            print("You can verify this by checking clustering/silhouette_score.png")
            # Plot
            plt.plot(range(2, max_clusters + 1), silhouette_scores, marker='o')
            plt.xlabel('Number of clusters')
            plt.ylabel('Silhouette Score')
            plt.title('Silhouette Method for Optimal k')
            plt.savefig('results/clustering/silhouette_score.png', dpi=300, bbox_inches='tight')
            plt.show()
        else:
            print("The optimal number of clusters could not be determined.")
            best_k = 1

    cluster_labels = fcluster(Z, best_k, criterion='maxclust')
    if len(set(cluster_labels)) < 2: # If the optimal number of clusters is 1, then the user should be warned
        print("The optimal number of clusters is 1. This means that the samples are too similar to each other.")
        print("Please, check the dendrogram to see the similarity between samples.")
    print("Sample labels and their corresponding cluster labels:")
    for sample, cluster in zip(adj_labels, cluster_labels):
        print(f"{sample}: Cluster {cluster}")
    print("Adding cluster information to 'samples.tsv' file...")

else:
    cluster_labels = fcluster(Z, k, criterion='maxclust')
    if len(set(cluster_labels)) < 2: # If the optimal number of clusters is 1, then the user should be warned
        print("The optimal number of clusters is 1. This means that the samples are too similar to each other.")
        print("Please, check the dendrogram to see the similarity between samples.")
    print("Sample labels and their corresponding cluster labels:")
    for sample, cluster in zip(adj_labels, cluster_labels):
        print(f"{sample}: Cluster {cluster}")
    print("Adding cluster information to 'samples.tsv' file...")

# Adding cluster division to "samples.tsv"
df = pd.read_csv('config/samples.tsv', sep='\t')
df.iloc[:, -1] = cluster_labels
df.to_csv('config/samples_supervised.tsv', sep='\t', index=False)

print("'samples_supervised.tsv' file successfully changed.")