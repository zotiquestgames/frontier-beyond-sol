# Crea una copia del grafo con nodi semplici
simple_graph = nx.DiGraph()

# Aggiungi nodi con solo l'identificatore come nome
for node, data in graph.nodes(data=True):
    simple_graph.add_node(node, label=data['name'])

# Aggiungi archi con lo stesso peso e attributi necessari
for u, v, data in graph.edges(data=True):
    simple_graph.add_edge(u, v, weight=data['weight'], highlight=data.get('highlight', False))

# Esporta in formato GraphML
nx.write_graphml(graph, "stellar_graph.graphml")

# Esporta in formato DOT (per Graphviz)
nx.nx_pydot.write_dot(graph, "stellar_graph.dot")
