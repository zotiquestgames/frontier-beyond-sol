import matplotlib.pyplot as plt

# Layout dei nodi (spring layout per chiarezza)
pos = nx.spring_layout(graph, seed=42)

# Crea una figura
plt.figure(figsize=(12, 12))

# Disegna i nodi: verdi per abitabili, blu per non abitabili
habitable_nodes = [n for n, d in graph.nodes(data=True) if d['habitable'] == 1]
non_habitable_nodes = [n for n, d in graph.nodes(data=True) if d['habitable'] == 0]

nx.draw_networkx_nodes(graph, pos, nodelist=habitable_nodes, node_color="green", label="Habitable", node_size=50)
nx.draw_networkx_nodes(graph, pos, nodelist=non_habitable_nodes, node_color="blue", label="Non-habitable", node_size=30)

# Disegna gli archi: rossi per quelli evidenziati, grigi per gli altri
highlighted_edges = [(u, v) for u, v, d in graph.edges(data=True) if d.get('highlight', False)]
nx.draw_networkx_edges(graph, pos, edgelist=highlighted_edges, edge_color="red", width=1.5, alpha=0.8)
nx.draw_networkx_edges(graph, pos, edge_color="lightgray", alpha=0.3)

# Disegna le etichette per i sistemi abitabili
labels = {n: d['name'] for n, d in graph.nodes(data=True) if d['habitable'] == 1}
nx.draw_networkx_labels(graph, pos, labels, font_size=8, font_color="black")

# Aggiungi legenda e titolo
plt.legend(scatterpoints=1)
plt.title("Graph Representation of Stellar Systems within 50 Light Years", fontsize=16)
plt.axis('off')

# Mostra il grafico
plt.tight_layout()
plt.show()
