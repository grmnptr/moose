#!/usr/bin/env python3
# * This file is part of the MOOSE framework
# * https://mooseframework.inl.gov
# *
# * All rights reserved, see COPYRIGHT for full restrictions
# * https://github.com/idaholab/moose/blob/master/COPYRIGHT
# *
# * Licensed under LGPL 2.1, please see LICENSE for details
# * https://www.gnu.org/licenses/lgpl-2.1.html

import os
import sys
import shutil
import numpy as np
import argparse
import importlib.util

import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, TensorDataset

print(torch.__version__)
torch.manual_seed(42)

# This is an example of a neural network that will determine the dislocation density
# based on strain and temperature (might not be physical, this is just an example)
class MyDDNet(nn.Module):
    def __init__(self):
        super(MyDDNet, self).__init__()

        # We save the normalization parameters into the buffer so we can
        # load them in C++
        self.register_buffer('mean', torch.tensor([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]))
        self.register_buffer('std', torch.tensor([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]))

        # We have two linear layers, we can add more if needed
        self.layer1 = nn.Linear(9, 32)
        self.layer2 = nn.Linear(32, 64)

        self.output_layer = nn.Linear(64, 1)

    def update_standardization(self, x):
        self.mean

    def forward(self, x):
        x = (x - self.mean) / self.std
        x = self.layer1(x)
        x = self.output_layer(x)
        return x

# Alright, time to train a neural net
model = MyDDNet()

# # These are the inputs for the learning process, not a lot so far.
# # Feel free to extend.
# strain = [0.0,  0.001, 0.003, 0.004]
# temperature =[300.0, 310.0, 320.0, 330.0]
# dd = [[0.1, 0.2, 0.3, 0.4]]

# # Converting things to tensors
# input_tensor = torch.tensor([strain,temperature]).T
# output_tensor = torch.tensor(dd).T

# # Setting up loss and optimizer
# criterion = nn.MSELoss()
# optimizer = optim.Adam(model.parameters(), lr=0.1)

# # Data loader just in case, we might not need it at this stage
# dataset = TensorDataset(input_tensor, output_tensor)
# dataloader = DataLoader(dataset, batch_size=4, shuffle=True)

# # Main training loop
# num_epochs = 100
# for epoch in range(num_epochs):
#     for inputs, targets in dataloader:
#         # Forward pass
#         outputs = model(inputs)
#         loss = criterion(outputs, targets)

#         # Backward pass and optimization
#         optimizer.zero_grad()
#         loss.backward()
#         optimizer.step()

#         print(f'Epoch [{epoch+1}/{num_epochs}], Loss: {loss.item():.4f}')

# print("Training complete.")

# print("Saving NN.")
# scripted_model = torch.jit.script(model.double())
# scripted_model.save("my_dd_net.pt")

StochasticControl = None
StochasticRunOptions = None

def tryImportStochasticControl(path = None):
    global StochasticControl
    global StochasticRunOptions

    if StochasticControl is not None:
        return True

    append_path = not (path is None or path in sys.path)
    if append_path:
        if not os.path.isdir(path):
            return False
        sys.path.append(path)

    if importlib.util.find_spec("moose_stochastic_tools") is None:
        if append_path:
            sys.path.pop(path)
        return False
    else:
        from moose_stochastic_tools import StochasticControl, StochasticRunOptions
        return True


if not tryImportStochasticControl():
    _moose_dir = os.environ.get("MOOSE_DIR", None)
    if not _moose_dir:
        _moose_dir = os.path.join(os.path.dirname(__file__), *([".."] * 5))
    _stm_python_path = os.path.abspath(
        os.path.join(_moose_dir, "modules", "stochastic_tools", "python")
    )
    tryImportStochasticControl(_stm_python_path)

if __name__ == "__main__":

    num_steps = np.random.randint(3, 10)
    num_cols = np.random.randint(3, 10)
    options = {}  # Options to insert into StochasticRunOptions

    # Arguments gathered from RUNAPP_COMMAND
    cmd = os.environ.get("RUNAPP_COMMAND")
    if cmd is None:
        sys.exit("Missing expected command variable RUNAPP_COMMAND")
    cmd = cmd.split()
    # Gather MPI options if running with MPI
    if "-n" in cmd:
        mpi_index = cmd.index("-n") - 1
        options["mpi_command"] = cmd[mpi_index]
        options["num_procs"] = int(cmd[mpi_index + 2])
    # Get executable and input file based on position of '-i' argument
    exec_index = cmd.index("-i") - 1
    executable = cmd[exec_index]
    input_file = cmd[exec_index + 2]
    # Everything after the input is cli_args
    options["cli_args"] = cmd[(exec_index + 3) :]

    # Import StochasticControl based on executable location (i.e. from conda)
    if StochasticControl is None:
        _exec_dir = os.path.dirname(os.path.abspath(shutil.which(executable)))
        _share_dir = os.path.abspath(os.path.join(_exec_dir, "..", "share"))
        _moose_stm_python = os.path.join(_share_dir, "stochastic_tools", "python")
        if not tryImportStochasticControl(_moose_stm_python):
            raise ModuleNotFoundError("Could not find MOOSE stochastic tools module python utilities.")

    # Arguments from cli
    cli_args = test_options()
    parameters = [
        "Postprocessors/rosenbrock/x["
        + ",".join([str(i) for i in range(num_cols)])
        + "]"
    ]
    qois = ["rosenbrock/value"]
    if cli_args.multi_output:
        parameters.append(f"Postprocessors/eggholder/x[{num_cols},{num_cols + 1}]")
        qois.append("eggholder/value")
        num_cols += 2
    options["multiapp_mode"] = StochasticRunOptions.MultiAppMode(cli_args.mode)

    # Get sampling matrices to test with
    matrices = [None] * num_steps
    for t in range(num_steps):
        # Test at least one matrix with only one entry
        num_rows = 1 if t == 0 else np.random.randint(1, 10)
        matrices[t] = 2.0 * (t + 1) * (np.random.random((num_rows, num_cols)) - 0.5)

    with StochasticControl(
        executable=executable,
        physics_input=input_file,
        parameters=parameters,
        quantities_of_interest=qois,
        options=StochasticRunOptions(**options),
    ) as runner:
        for matrix in matrices:
            result = runner(matrix)
            # Need to reshape since result could be a float or 1-D array
            result = np.array(result).reshape(
                (matrix.shape[0], 2 if cli_args.multi_output else 1)
            )
            compare(matrix, result)
