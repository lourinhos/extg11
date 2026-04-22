import networkx as nx
import matplotlib.pyplot as plt
import torch.nn as nn


def visualize_neural_network(model, input_dim):
    G = nx.DiGraph()
    
    # Input layer
    input_nodes = [f"x{i+1}" for i in range(input_dim)]
    for node in input_nodes:
        G.add_node(node, layer=0)

    # Hidden & Output layers
    prev_nodes = input_nodes
    layer = 1
    for module in model.model:
        if isinstance(module, nn.Linear):
            layer_nodes = [f"layer{layer}_n{i+1}" for i in range(module.out_features)]
            for node in layer_nodes:
                G.add_node(node, layer=layer)
            for prev_node in prev_nodes:
                for node in layer_nodes:
                    G.add_edge(prev_node, node)
            prev_nodes = layer_nodes
            layer += 1
    
    # Output node
    output_node = "y"
    G.add_node(output_node, layer=layer)
    for prev_node in prev_nodes:
        G.add_edge(prev_node, output_node)

    # Layout
    pos = nx.multipartite_layout(G, subset_key="layer")
    plt.figure(figsize=(8, 6))
    nx.draw(G, pos, with_labels=True, node_size=2000, node_color="lightblue", font_size=10, edge_color="gray")
    plt.show()


def visualize_perceptron(model, input_dim):
    """
    Visualizes a single-layer perceptron using NetworkX.

    Parameters:
        model (Perceptron): Trained perceptron model.
        input_dim (int): Number of input features.
    """
    G = nx.DiGraph()  # Directed graph

    # Define node positions
    pos = {}
    input_layer_y = 1  # Height of input layer
    output_layer_y = 0  # Height of output layer

    # Add input neurons
    for i in range(input_dim):
        G.add_node(f"x{i+1}", layer="input")
        pos[f"x{i+1}"] = (i, input_layer_y)

    # Add bias neuron
    G.add_node("bias", layer="input")
    pos["bias"] = (input_dim, input_layer_y)

    G.add_node("Learning Rule", layer="act")
    pos["Learning Rule"] = (input_dim / 2, 0.5)

    # Add output neuron
    G.add_node("Output", layer="output")
    pos["Output"] = (input_dim / 2, output_layer_y)

    # Add edges with weights
    for i in range(input_dim):
        G.add_edge(f"x{i+1}", "Learning Rule", weight=model.weights[i].item())  # Convert tensor to number
    G.add_edge("bias", "Learning Rule", weight=model.weights[-1].item())  # Bias weight
    G.add_edge("Learning Rule", "Output")  # Output weight

    # Draw the network
    plt.figure(figsize=(6, 4))
    nx.draw(G, pos, with_labels=True, node_color=["lightblue"]*input_dim + ["gray", "lightgreen", "lightgreen"],
            node_size=2000, edge_color="black", font_size=10, font_weight="bold")

    # Draw edge labels (weights)
    edge_labels = {(u, v): f"{d['weight']:.2f}" for u, v, d in G.edges(data=True) if v != "Output"}
    nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels, font_size=10)

    plt.title("Perceptron Visualization")
    plt.show()



