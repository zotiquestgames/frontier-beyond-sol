import pandas as pd
import numpy as np
from scipy.spatial.distance import pdist, squareform
import networkx as nx
import matplotlib.pyplot as plt

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

# Crea un grafo vuoto
graph = nx.Graph()

# Aggiungi i nodi con attributi (nome e abitabilità)
for i, name in enumerate(names):
    graph.add_node(i, name=name, habitable=habitability[i])

# Aggiungi archi con priorità ai pianeti abitabili
for i in range(len(names)):
    distances_from_node = distances[i].copy()
    distances_from_node[i] = np.inf  # Escludi la distanza a sé stesso

    # Ordina le stelle per distanza
    sorted_indices = np.argsort(distances_from_node)
    closest_two = sorted_indices[:2]  # Prendi le due più vicine

    # Priorità ai pianeti abitabili
    habitable_neighbors = [idx for idx in closest_two if habitability[idx] == 1]
    if habitable_neighbors:
        habitable_idx = habitable_neighbors[0]
        graph.add_edge(i, habitable_idx, weight=distances_from_node[habitable_idx], highlight=True)

        for idx in closest_two:
            if idx != habitable_idx:
                graph.add_edge(i, idx, weight=distances_from_node[idx], highlight=False)
                break
    else:
        for idx in closest_two:
            graph.add_edge(i, idx, weight=distances_from_node[idx], highlight=False)

# Forza il Sole al centro del layout
sun_index = names.index("Sol")
fixed_positions = {sun_index: (0, 0)}
pos = nx.spring_layout(graph, seed=42, pos=fixed_positions, fixed=fixed_positions.keys())

# Crea una figura
plt.figure(figsize=(12, 12))

# Disegna i nodi: verdi per abitabili, blu per non abitabili
habitable_nodes = [n for n, d in graph.nodes(data=True) if d['habitable'] == 1]
non_habitable_nodes = [n for n, d in graph.nodes(data=True) if d['habitable'] == 0]

nx.draw_networkx_nodes(graph, pos, nodelist=habitable_nodes, node_color="green", label="Habitable", node_size=50)
nx.draw_networkx_nodes(graph, pos, nodelist=non_habitable_nodes, node_color="blue", label="Non-habitable", node_size=30)

# Disegna gli archi: tutti tratteggiati
highlighted_edges = [(u, v) for u, v, d in graph.edges(data=True) if d.get('highlight', False)]
nx.draw_networkx_edges(graph, pos, edgelist=highlighted_edges, edge_color="red", width=1.5, alpha=0.8, style='dashed')
nx.draw_networkx_edges(graph, pos, edge_color="lightgray", alpha=0.3, style='dashed')

# Disegna le etichette
labels = {n: d['name'] for n, d in graph.nodes(data=True)}
nx.draw_networkx_labels(graph, pos, labels, font_size=8, font_color="black")

# Aggiungi legenda e titolo
plt.legend(scatterpoints=1)
plt.title("Graph Representation of Stellar Systems within 50 Light Years", fontsize=16)
plt.axis('off')

# Salva in formato SVG
output_svg_path = "stellar_graph.svg"
plt.savefig(output_svg_path, format="svg", bbox_inches="tight")

# Mostra il grafico
plt.tight_layout()
plt.show()

# Esportazioni opzionali
nx.write_graphml(graph, "stellar_graph.graphml")

# Crea una copia del grafo con nodi semplici
simple_graph = nx.DiGraph()

# Aggiungi nodi con solo l'identificatore come nome
for node, data in graph.nodes(data=True):
    simple_graph.add_node(node, label=data['name'])

# Aggiungi archi con lo stesso peso e attributi necessari
for u, v, data in graph.edges(data=True):
    simple_graph.add_edge(u, v, weight=data['weight'], highlight=data.get('highlight', False))

# Esporta in formato DOT
nx.nx_pydot.write_dot(simple_graph, "stellar_graph.dot")

# Esporta in formato GraphML
nx.write_graphml(graph, "stellar_graph.graphml")


