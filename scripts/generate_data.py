#!/usr/bin/env python3

# PyTorch imports
import torch
import torch.nn as nn
import torch.nn.functional as F

# NumPy imports
import numpy as np

# SymPy
from sympy import *

#
# Parameters
#

# output data file with D_in + D_out columns
filename = 'data.txt'

# output dimension (usually one free energy)
D_in = 2

# output dimension (usually one free energy)
D_out = 1

# gradient outputs
D_grad = D_in * D_out

#
# The model
#

# free energy expression
in_vars = symbols("T c_p")
T, cp = in_vars
F = cp*(1-cp) + 0.001*(T*300+400)*(cp*log(cp)+(1-cp)*log(1-cp))

# derivatives
D = [diff(F,var) for var in in_vars]

x = [0] * len(in_vars)
map = {}
with open(filename, "w") as f:
    for x[0] in np.arange(0.15,0.25,0.02):
        for x[1] in np.arange(0.05,0.95,0.001):
            # substitution rule
            for i in range(len(in_vars)):
                map[in_vars[i]] = x[i]

            # input vec and output node
            row = x + [F.evalf(subs=map)]

            # evaluate gradients
            row += [d.evalf(subs=map) for d in D]

            print(row)
            np.savetxt(f, [row])
