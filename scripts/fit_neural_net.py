#!/usr/bin/env python

import os
import sys
import argparse

# PyTorch imports
import torch
import torch.nn as nn
import torch.nn.functional as F

# NumPy imports
import numpy as np

#
# Parameters
#

cmd = ' '.join(sys.argv)

parser = argparse.ArgumentParser(description='Fit a neural network to a Gibbs free energy.')
parser.add_argument('datafile', help='Text file with columns for the inputs, free energy, and (optionally) chemical potential data.')
parser.add_argument('inputs', type=int, help="Number of input nodes")
parser.add_argument('outputs', type=int, help="Number of input nodes")
parser.add_argument('-c', '--use_chemical_potentials', action="store_true")
parser.add_argument('-o', '--output_status', type=int, default=10, help="Epoch interval for outputting the loss function value to screen and disk")
parser.add_argument('-O', '--output_model', type=int, default=1000, help="Epoch interval for outputting the NN torch model to disk")
parser.add_argument('-e', '--epochs', type=int, default=50000, help="Epochs to run for training (50000)")
parser.add_argument('-m', '--mini_batch', type=int, default=128, help="Mini batch size (128)")
parser.add_argument('-H', '--hidden_layer_nodes', type=int, default=20, help="Number of nodes per hidden layer (20)")
parser.add_argument('-n', '--hidden_layer_count', type=int, default=2, help="Number of hidden layers (2)")
parser.add_argument('-r', '--learning_rate', type=float, default=1e-5, help="Learning rate meta parameter (1e-5)")
parser.add_argument('-G', '--disable_gpu', action="store_true", help="Disable GPU detection and learn on the CPU")
parser.add_argument('-a', '--activation', choices=['softsign','sigmoid','tanh'], default='softsign', help="Activation function")
args = parser.parse_args()
print(args)

# use the derivatives of the output w.r.t. the inputs (i.e. chemical potentials)
use_chemical_potentials = args.use_chemical_potentials
if use_chemical_potentials:
    print("Using chemical potential data in loss function")

# use cuda if available
gpu = torch.cuda.is_available() and not args.disable_gpu

# data file with D_in + D_out columns
filename = args.datafile

# Activation function
if args.activation == 'sigmoid' :
    activation = torch.nn.Sigmoid
elif args.activation == 'softsign' :
    activation = torch.nn.Softsign
elif args.activation == 'tanh' :
    activation = torch.nn.Tanh
else :
    print("invalid activation function '%s'" % args.activation)

# input dimension (temperature and concentrations)
D_in = args.inputs

# output dimension (usually one free energy)
D_out = args.outputs

# gradient outputs
D_grad = D_in * D_out

# write model to disk every this many steps
output_model = args.output_model

# write model to disk every this many steps
output_status = args.output_status

# load data file
data = np.loadtxt(filename, unpack=True)
N = len(data[0])
print("Training set size is %d" % N)

# mini batch size
batch_size = args.mini_batch
print("Mini batch size is %d" % batch_size)

# Number of epochs to run
n_epochs = args.epochs
print("Training for %d epochs" % n_epochs)

# learning rate meta parameter
learning_rate = args.learning_rate
print("Learning rate %f" % learning_rate)

# number of hidden layers...
hidden_layer_count = args.hidden_layer_count

# ...and their node count
H = args.hidden_layer_nodes

# open training log
log = open("model.log", "a")
log.write("# %s\n" % cmd)
log.flush()

# Create tensors to hold inputs and outputs (and gradients)
if gpu:
    x = torch.from_numpy(data[0:D_in]).transpose(0,1).float().cuda()
    y = torch.from_numpy(data[D_in:D_in+D_out]).transpose(0,1).float().cuda()
    grad = torch.from_numpy(data[D_in+D_out:]).transpose(0,1).float().cuda()
else:
    x = torch.Tensor(data[0:D_in]).transpose(0,1)
    y = torch.Tensor(data[D_in:D_in+D_out]).transpose(0,1)
    grad = torch.Tensor(data[D_in+D_out:]).transpose(0,1)

# get min and max of input data
in_bounds = [(np.amin(d), np.amax(d)) for d in data[:D_in]]
in_norm = [(b[1]-b[0], b[0]/(b[1]-b[0])) for b in in_bounds]
out_bounds = [(np.amin(d), np.amax(d)) for d in data[D_in:D_in+D_out]]
out_norm = [(b[1]-b[0], b[0]) for b in out_bounds]

print(D_in, D_in+D_out)
print(np.amin(data[D_in:D_in+D_out][0]))
print(data[D_in:D_in+D_out][0])

def weights_init(m):
    classname = m.__class__.__name__
    # for every Linear layer in a model..
    if classname.find('Linear') != -1:
        # get the number of the inputs
        n = m.in_features
        m.weight.data.normal_(std=1.0/n)
        m.bias.data.fill_(0)

# adjust input and output weights and biases for lacking data normalization
def adjust_weights():
    params = list(model.parameters())
    n_layers = int(len(params) / 2)

    # input weights
    p = params[0]
    for h in p:
        for i in range(D_in):
            h[i].data /= in_norm[i][0]

    # input biases
    p = params[1]
    for h in p:
        for i in range(D_in):
            h.data -= in_norm[i][1]

    # output weights
    p = params[n_layers*2-2]
    for i in range(D_out):
      for h in p[i]:
            h.data *= out_norm[i][0]

    # output biases
    p = params[n_layers*2-1]
    for i in range(D_out):
      p[i].data += out_norm[i][1]

# adjust input and output weights and biases for lacking data normalization
def shift_output():
    params = list(model.parameters())
    n_layers = int(len(params) / 2)

    y_pred = model.forward(x)
    if gpu:
        yd = np.transpose(y.cpu().detach().numpy())
        ypd = np.transpose(y_pred.cpu().detach().numpy())
    else:
        yd = np.transpose(y.detach().numpy())
        ypd = np.transpose(y_pred.detach().numpy())

    # output biases
    p = params[n_layers*2-1]
    for i in range(D_out):
      p[i].data += np.mean(yd[i]-ypd[i])

if os.path.exists('model.dat') :
    # Load the model from file
    print("Resuming training of existing model")
    model = torch.load('model.dat')
    # model = model.load_state_dict(torch.load('model.dat'))
else:
    # Use the nn package to define our model and loss function.
    print("Seting up new model with %d hidden layers with %d nodes each." %
          (hidden_layer_count, H))

    layers = [torch.nn.Linear(D_in, H), activation()]
    for i in range(1, hidden_layer_count):
        layers += [torch.nn.Linear(H, H), activation()]
    layers += [torch.nn.Linear(H, D_out)]
    model = torch.nn.Sequential(*layers)

    model.apply(weights_init)
    adjust_weights()

if gpu:
    model = model.float().cuda()

# shift outputs by recalculating output layer bias
shift_output()

# loss_fn = torch.nn.MSELoss(reduction='sum')

# Use the optim package to define an Optimizer that will update the weights of
# the model for us. Here we will use Adam; the optim package contains many other
# optimization algoriths. The first argument to the Adam constructor tells the
# optimizer which Tensors it should update.
optimizer = torch.optim.Adam(model.parameters(), lr=learning_rate)
for epoch in range(n_epochs):

    # generate a new permutation of the training set
    permutation = torch.randperm(x.size()[0])

    # iterate over mini batches
    for i in range(0,x.size()[0], batch_size):

      indices = permutation[i:i+batch_size]
      batch_x = torch.autograd.Variable(x[indices], requires_grad=True)
      batch_y = torch.autograd.Variable(y[indices], requires_grad=False)

      # Before the backward pass, use the optimizer object to zero all of the
      # gradients for the variables it will update (which are the learnable
      # weights of the model). This is because by default, gradients are
      # accumulated in buffers( i.e, not overwritten) whenever .backward()
      # is called. Checkout docs of torch.autograd.backward for more details.
      optimizer.zero_grad()
      y_pred = model.forward(batch_x)

      if use_chemical_potentials:
          ones = torch.ones(y_pred.shape)
          if gpu:
              ones = ones.float().cuda()
          # Forward pass: compute predicted y by passing x to the model.
          y_pred.backward(ones, retain_graph=True)
          grad_pred = torch.autograd.grad(y_pred, batch_x,
                        grad_outputs=y_pred.data.new(y_pred.shape).fill_(1),
                        create_graph=True)[0]

          optimizer.zero_grad()

          batch_grad = torch.autograd.Variable(grad[indices], requires_grad=False)

          # skip temperature derivative
          batch_grad = batch_grad[:,1:]
          grad_pred = grad_pred[:,1:]

          loss = torch.mean((y_pred - batch_y) ** 2 + (grad_pred - batch_grad) ** 2)

      else:
          loss = torch.mean((y_pred - batch_y) ** 2)

      # print learning status
      if epoch % output_status == 0 and i == 0:
          print(epoch, loss.item())
          log.write("%d %E\n" % (epoch, loss.item()))
          log.flush()

      # save model once in a while to allow for resuming of the training
      if epoch > 0 and epoch % output_model == 0 and i == 0:
          torch.save(model, 'model.dat')
          #torch.save(model.state_dict(), 'model.dat')
          print("Model saved.")

      # Backward pass: compute gradient of the loss with respect to model
      # parameters
      loss.backward()

      # Calling the step function on an Optimizer makes an update to its
      # parameters
      optimizer.step()

log.close()

if use_chemical_potentials or True:
    # Forward pass: compute predicted y by passing x to the model.
    xg = torch.autograd.Variable(x, requires_grad=True)
    y_pred = model.forward(xg)
    ones = torch.ones(y_pred.shape)
    if gpu:
        ones = ones.float().cuda()
    y_pred.backward(ones, retain_graph=True)
    grad_pred = torch.autograd.grad(y_pred, xg,
                  grad_outputs=y_pred.data.new(y_pred.shape).fill_(1),
                  create_graph=True)[0]
    if gpu:
        xd = x.cpu().detach().numpy()
        yd = y_pred.cpu().detach().numpy()
        gd = grad_pred.cpu().detach().numpy()
    else:
        xd = x.detach().numpy()
        yd = y_pred.detach().numpy()
        gd = grad_pred.detach().numpy()
    np.savetxt('out.txt', np.concatenate((xd, yd, gd), axis=1))
else:
    y_pred = model.forward(x)
    if gpu:
        xd = x.cpu().detach().numpy()
        yd = y_pred.cpu().detach().numpy()
    else:
        xd = x.detach().numpy()
        yd = y_pred.detach().numpy()
    np.savetxt('out.txt', np.concatenate((xd, yd), axis=1))

torch.save(model, 'model.dat')
#torch.save(model.state_dict(), 'model.dat')
