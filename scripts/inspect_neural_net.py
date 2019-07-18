#!/usr/bin/env python

# PyTorch imports
import torch

# Load the model from file
model = torch.load('model.dat')

# number of layers excluding input layer
params = list(model.parameters())

for p in params:
  print(p)
