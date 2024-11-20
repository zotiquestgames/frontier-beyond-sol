import numpy as np
import networkx as nx

# Crea un grafo vuoto
graph = nx.Graph()

# Aggiungi i nodi con attributi (nome e abitabilità)
for i, name in enumerate(names):
    graph.add_node(i, name=name, habitable=habitability[i])

# Aggiungi archi con priorità ai pianeti abitabili
for i in range(len(names)):
    # Copia delle distanze per il nodo corrente
    distances_from_node = distances[i].copy()
    distances_from_node[i] = np.inf  # Escludi la distanza a sé stesso

    # Ordina le stelle per distanza
    sorted_indices = np.argsort(distances_from_node)
    closest_two = sorted_indices[:2]  # Prendi le due più vicine

    # Priorità ai pianeti abitabili
    habitable_neighbors = [idx for idx in closest_two if habitability[idx] == 1]
    if habitable_neighbors:
        # Collega il più vicino abitabile
        habitable_idx = habitable_neighbors[0]
        graph.add_edge(i, habitable_idx, weight=distances_from_node[habitable_idx], highlight=True)

        # Collega anche il successivo nella lista delle distanze
        for idx in closest_two:
            if idx != habitable_idx:
                graph.add_edge(i, idx, weight=distances_from_node[idx], highlight=False)
                break
    else:
        # Collega le due stelle più vicine senza priorità
        for idx in closest_two:
            graph.add_edge(i, idx, weight=distances_from_node[idx], highlight=False)
