#!/usr/bin/env python3

#
# Import a saved model.dat
#

# PyTorch imports
import torch
import torch.nn as nn
import torch.nn.functional as F

# NumPy imports
import numpy as np

# Load the model from file
model = torch.load('model.dat')

# number of layers excluding input layer
params = list(model.parameters())
n_layers = len(params) // 2

with open('weights_biases.txt', "w") as f:
    # write number of layers
    f.write("%d\n" % n_layers)

    # loop over weights
    for layer in range(n_layers):
        p = params[layer*2]
        pp = p.detach().cpu().numpy().transpose()
        # print size
        x = pp.shape[0]
        y = pp.shape[1]
        f.write("%d %d\n" % (x, y))
        for i in range(x):
            np.savetxt(f, pp[i])

    # loop over biases
    for layer in range(n_layers):
        p = params[layer*2+1]
        pp = p.detach().cpu().numpy()
        # print size
        x = p.size()[0]
        f.write("%d\n" % x)
        np.savetxt(f, pp)
