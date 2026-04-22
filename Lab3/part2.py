import torch
from utils import Perceptron, MultiLayerPerceptron


X = torch.tensor([[0, 0], [0, 1], [1, 0], [1, 1]], dtype=torch.float32)
y = torch.tensor([0, 1, 1, 1], dtype=torch.float32)  
  
perceptron = Perceptron(input_dim=2)

res = perceptron.train_perceptron(X, y, epochs=100, lr=0.8)

print(f"Final predictions: {perceptron.forward(X)}")

print(f"Total epochs: {res['epochs']}")

   