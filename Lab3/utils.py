import torch
import torch.nn as nn
from visualize_network import visualize_perceptron
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np

class HardLimit(torch.nn.Module):
    def forward(self, x):
        return torch.where(x >= 0, torch.tensor(1.0, device=x.device), torch.tensor(0.0, device=x.device))

class Perceptron(torch.nn.Module):
    """ Single-layer perceptron"""
    def __init__(self, input_dim):
        super(Perceptron, self).__init__()
        self.weights = torch.randn(input_dim + 1, dtype=torch.float32)  # +1 for bias
        self.activation = HardLimit()

    def forward(self, x):
        x = torch.cat((x, torch.ones(x.shape[0], 1)), dim=1)  # Add bias term
        net_input = torch.matmul(x, self.weights)  # Linear combination
        return self.activation(net_input)  # Apply step function

    def train_perceptron(self, X, y, epochs=10, lr=0.1):
        """ Train the perceptron"""
        X = torch.cat((X, torch.ones(X.shape[0], 1)), dim=1)  # Add bias term
        
        for epoch in range(epochs):
            errors = 0
            for i in range(X.shape[0]):
                output = self.activation(torch.dot(X[i], self.weights))  # Compute output
                error = y[i] - output  # Compute error
                if error != 0:
                    self.weights += lr * error * X[i]  # Update weights
                    errors += 1
            print(f"Epoch {epoch+1}: Misclassified {errors} samples")

            if errors == 0:
                print("Training complete - all samples classified correctly!")
                break
        return {"weights": self.weights, "epochs": epoch+1}

    def visualize(self):
        visualize_perceptron(self, input_dim=2)


class MultiLayerPerceptron(nn.Module):
    def __init__(self, hiddenSizes=2):
        super(MultiLayerPerceptron, self).__init__()

        if isinstance(hiddenSizes, int):
            self.hiddenSizes = [hiddenSizes]
        else:
            self.hiddenSizes = hiddenSizes

        layers = []
        for i, n in enumerate(self.hiddenSizes):
            if i == 0:
                layers.append(nn.Linear(2, n, bias=True))
            else:
                layers.append(nn.Linear(self.hiddenSizes[i-1], n, bias=True))
            layers.append(nn.ReLU())

        self.model = nn.Sequential(*layers, nn.Linear(self.hiddenSizes[-1], 1, bias=True), nn.Sigmoid())

        # Initialize weights
        for layer in self.model:
            if isinstance(layer, nn.Linear):
                nn.init.xavier_uniform_(layer.weight)

    def forward(self, x):
        return self.model(x)
    
    def train_mlp(self, X, y, epochs=10, lr=0.1):
        criterion = nn.BCELoss()
        optimizer = torch.optim.SGD(self.parameters(), lr=lr)

        if len(y.shape) == 1:
            y = y.unsqueeze(1)

        for epoch in range(epochs):
            optimizer.zero_grad()
            output = self.forward(X)
            loss = criterion(output, y)
            error = sum((output > 0.5).int().squeeze() != y.int().squeeze())

            print(f"Epoch: {epoch+1} Loss: {loss.item()} Error: {error}")
            if error == 0:
                print("Training complete - all samples classified correctly!")
                break

            loss.backward()
            optimizer.step()

        return {"epochs": epoch+1}

    def predict(self, x):
        return (self.forward(x) > 0.5).int().squeeze()
    
    def visualize(self):
        G = nx.DiGraph()
        pos = {}
        node_colors = []  # List to store colors for each node

        # Input layer
        input_nodes = [f"x{i+1}" for i in range(2)]  # Assuming 2 input features
        for i, node in enumerate(input_nodes):
            G.add_node(node, layer=0)
            # Center the input nodes horizontally
            pos[node] = (i - (len(input_nodes) - 1) / 2, len(self.hiddenSizes) + 1)
            node_colors.append("lightgreen")  # Input nodes color

        # Hidden layers
        prev_nodes = input_nodes
        for layer_idx, num_neurons in enumerate(self.hiddenSizes):
            current_nodes = [f"h{layer_idx+1}_n{i+1}" for i in range(num_neurons)]
            for i, node in enumerate(current_nodes):
                G.add_node(node, layer=layer_idx + 1)
                # Center the hidden nodes horizontally
                pos[node] = (i - (len(current_nodes) - 1) / 2, len(self.hiddenSizes) - layer_idx)
                node_colors.append("lightblue")  # Hidden nodes color
                for prev_node in prev_nodes:
                    G.add_edge(prev_node, node)
            prev_nodes = current_nodes

        # Output layer
        output_node = "y"
        G.add_node(output_node, layer=len(self.hiddenSizes) + 1)
        # Center the output node horizontally
        pos[output_node] = (0, 0)
        node_colors.append("lightcoral")  # Output node color
        for prev_node in prev_nodes:
            G.add_edge(prev_node, output_node)

        # Draw the graph
        plt.figure(figsize=(10, 6))
        nx.draw(
            G,
            pos,
            with_labels=True,
            node_size=2000,
            node_color=node_colors,  # Use the node_colors list
            font_size=10,
            edge_color="gray",
        )
        plt.show()


def plot_db(model, X_train, y_train):
    # Generate grid points
    x_min, x_max = -0.5, 1.5
    y_min, y_max = -0.5, 1.5
    xx, yy = np.meshgrid(np.linspace(x_min, x_max, 100), np.linspace(y_min, y_max, 100))
    grid = np.c_[xx.ravel(), yy.ravel()]

    grid_tensor = torch.tensor(grid, dtype=torch.float32)
    with torch.no_grad():
        predictions = model(grid_tensor)  
        class_preds = (predictions > 0.5).numpy().reshape(xx.shape)

    # Plot decision boundary
    plt.contourf(xx, yy, class_preds, alpha=0.5, cmap=plt.cm.coolwarm)

    # Plot training points
    plt.scatter(X_train[:, 0], X_train[:, 1], c=y_train, edgecolors='k', cmap=plt.cm.coolwarm)

    plt.xlabel("X1")
    plt.ylabel("X2")
    plt.title("MLP Decision Boundary")
    plt.show()


