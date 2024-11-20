import pandas as pd
import numpy as np
from scipy.spatial.distance import pdist, squareform

# Carica i dati dal file CSV
file_path = 'HabHYG50ly.csv'  # Sostituisci con il percorso corretto
stellar_data = pd.read_csv(file_path)

# Filtra i dati per coordinate valide
stellar_data_filtered = stellar_data.dropna(subset=['Xg', 'Yg', 'Zg', 'Display Name'])

# Estrai dati rilevanti
coordinates = stellar_data_filtered[['Xg', 'Yg', 'Zg']].to_numpy()
names = stellar_data_filtered['Display Name'].tolist()
habitability = stellar_data_filtered['Hab?'].fillna(0).astype(int).tolist()

# Calcola le distanze tra tutte le stelle
distances = squareform(pdist(coordinates))
