################
#This code gets the coordinates of the B-cells and tumor cells and find the aggregates of B-cells.
#Last modified: Mohammad Yosofvand, 09/05/2024
################
import os
import numpy as np
import pandas as pd
from scipy.spatial import ConvexHull, distance
from sklearn.cluster import DBSCAN
import scipy.stats
from sklearn.decomposition import PCA
from math import pi

def calculate_cluster_characteristics(cell_centroids, clusters, tumor_cells):
    cluster_characteristics = []
    unique_clusters = [c for c in np.unique(clusters) if c != -1]
    for cluster_id in unique_clusters:
        cluster_points = cell_centroids[clusters == cluster_id]
        hull = ConvexHull(cluster_points)
        hull_area = hull.volume
        centroid = np.mean(cluster_points, axis=0)
        roundness = (4 * pi * hull_area) / (hull.area ** 2)
        skewness = scipy.stats.skew(cluster_points, axis=0)
        pca = PCA(n_components=2).fit(cluster_points)
        aspect_ratio = pca.explained_variance_ratio_[0] / pca.explained_variance_ratio_[1]
        transformed_points = pca.transform(cluster_points)
        min_caliper = np.max(transformed_points[:, 0]) - np.min(transformed_points[:, 0])
        max_caliper = np.max(transformed_points[:, 1]) - np.min(transformed_points[:, 1])
        min_max_caliper_ratio = min_caliper / max_caliper
        hull_points = cluster_points[hull.vertices]
        distances_to_centroid = distance.cdist([centroid], hull_points)[0]
        average_radius = np.mean(distances_to_centroid)
        distances = distance.cdist([centroid], tumor_cells[['X', 'Y']].values)[0]
        nearest_distances = np.partition(distances, 100)[:100] if len(distances) > 100 else distances
        average_distance = np.mean(nearest_distances)
        
        cluster_characteristics.append({
            'Cluster ID': cluster_id,
            'Size': len(cluster_points),
            'Centroid X': centroid[0],
            'Centroid Y': centroid[1],
            'Area': hull_area,
            'Roundness': roundness,
            'Skewness X': skewness[0],
            'Skewness Y': skewness[1],
            'Aspect Ratio': aspect_ratio,
            'Min Max Caliper Ratio': min_max_caliper_ratio,
            'Distance': average_distance,
            'Radius': average_radius
        })
    return cluster_characteristics

def process_files(directory):
    all_characteristics = []
    for filename in os.listdir(directory):
        if filename.endswith(".csv"):
            file_path = os.path.join(directory, filename)
            df = pd.read_csv(file_path)
            df_tumor = df[(df['region'] == 'Tumor') & (df['majortype'] == 'TumorEpithelial')]
            df = df[df['majortype'].isin(['Bcell'])]
            if df.empty:
                continue  
            cell_centroids = df[['X', 'Y']].values
            dbscan = DBSCAN(eps=225, min_samples=100, metric='euclidean')
            clusters = dbscan.fit_predict(cell_centroids)
            characteristics = calculate_cluster_characteristics(cell_centroids, clusters, df_tumor)
            for char in characteristics:
                char['File Name'] = filename[:-4]
                all_characteristics.append(char)

    df_columns = ['File Name', 'Cluster ID', 'Size', 'Centroid X', 'Centroid Y', 'Area', 'Roundness', 
                  'Skewness X', 'Skewness Y', 'Aspect Ratio', 'Min Max Caliper Ratio', 'Distance',
                  'Radius']
    return pd.DataFrame(all_characteristics, columns=df_columns)

# Replace 'your_directory_path' with the actual directory path containing your CSV files
directory_path = '/path/to/csv/files/'
results_df = process_files(directory_path)

# Save to CSV
results_df.to_csv('cluster_characteristics.csv', index=False)
