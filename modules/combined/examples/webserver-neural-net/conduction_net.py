#!/usr/bin/env python3
import os, shlex, sys
from MooseControl import MooseControl
import csv
import math

import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, TensorDataset


import phono3py as p3
import numpy as np

def compute_thermal_conductivity(temperatures):

    # Load the same params file the CLI would read
    ph3 = p3.load("phono3py_params_NaCl.yaml.xz", produce_fc=False, log_level=2)

    # --mesh 50  ⇢ Γ-centered 50×50×50 grid
    ph3.mesh_numbers = 50  # API accepts an int or [n1,n2,n3]. :contentReference[oaicite:1]{index=1}

    temperature = 273

    # Prepare 3-phonon interaction on that grid (the key missing step)
    ph3.init_phph_interaction()  # must precede run_thermal_conductivity. :contentReference[oaicite:2]{index=2}

    # --br  ⇢ RTA (BTE in relaxation-time approximation). The API defaults to RTA.
    ph3.run_thermal_conductivity(temperatures=temperatures, is_LBTE=False)

    # (Optional) quick sanity print
    # ... after ph3.run_thermal_conductivity(...)
    tc = ph3.thermal_conductivity

    # Temperatures array (same as the 'temperature' dataset in the HDF5)
    T = tc.temperatures  # shape (nT,)

    # 6-component vector per T: [κxx, κyy, κzz, κyz, κxz, κxy]
    # Matches the HDF5 layout shown in the docs.
    kappa = tc.kappa

    # Isotropic average = (κxx+κyy+κzz)/3
    kappa = tc.kappa[0]     # shape (nT, 6)
    k_iso = [np.mean(l[:3]) for l in kappa]

    # If you ran only at T=300 K:
    print("Temperatures (K):", T.tolist())
    print("κ_iso [W/mK]:", k_iso)

    return k_iso

# This is an example of a neural network that defines the thermal conductivity
class MyTCNet(nn.Module):
    def __init__(self):
        super(MyTCNet, self).__init__()

        # We save the normalization parameters into the buffer so we can
        # load them in C++
        self.register_buffer('mean', torch.tensor([0.0]))
        self.register_buffer('std', torch.tensor([1.0]))

        # We have two linear layers, we can add more if needed
        self.layer1 = nn.Linear(1, 2)
        self.layer2 = nn.Linear(2, 4)
        self.output_layer = nn.Linear(4, 1)

    def update_standardization(self, x):
        self.mean = np.mean(x)
        self.std = np.std(x)

    def forward(self, x):
        x = (x - self.mean) / self.std
        x = self.layer1(x)
        x = self.layer2(x)
        x = self.output_layer(x)
        return x

def train_nn(model, input, output, num_epochs, learning_rate):

    # Converting things to tensors
    input_tensor = torch.tensor(input).T
    output_tensor = torch.tensor(output).T

    # Setting up loss and optimizer
    criterion = nn.MSELoss()
    optimizer = optim.Adam(model.parameters(), lr=learning_rate)

    # Data loader just in case, we might not need it at this stage
    dataset = TensorDataset(input_tensor, output_tensor)
    dataloader = DataLoader(dataset, batch_size=len(input), shuffle=False)

    # Main training loop
    for epoch in range(num_epochs):
        for inputs, targets in dataloader:
            # Forward pass
            outputs = model(inputs)
            loss = criterion(outputs, targets)

            # Backward pass and optimization
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()

            print(f'Epoch [{epoch+1}/{num_epochs}], Loss: {loss.item():.4f}')

    print("Training complete.")

    print("Saving NN.")
    scripted_model = torch.jit.script(model)
    scripted_model.save("my_tc_net.pt")

# Passed into the base controller to run the MooseControl
def run_moose(control):
    # Get the command we should run
    # You'll hit this if you don't run a RunApp-derived Tester or
    # don't run it with the "command_proxy" option
    RUNAPP_COMMAND = os.environ.get('RUNAPP_COMMAND')
    if RUNAPP_COMMAND is None:
        sys.exit('Missing expected command variable RUNAPP_COMMAND')

    moose_command = shlex.split(RUNAPP_COMMAND)
    control = MooseControl(moose_command=moose_command,
                           moose_control_name=control)
    control.initialize()

    value = []
    # Run the user test
    try:
        # Wait until the server is ready, where it should be at INITIAL
        control.wait('INITIAL')

        # Tell MOOSE to continue with the solve
        control.setContinue()

        # Wait, where we shold be at TIMESTEP_BEGIN
        control.wait('FINAL')

        # Get the current value of the postprocessor and compare to the gold
        value = control.getReporterValue('temps/u')

        # Tell MOOSE to continue with the solve
        control.setContinue()
    except:
        control.kill()
        raise

    # Wait for the webserver to stop listening
    control.finalize()

    return value

def select_temperatures(temperatures, current_temperatures, distance):

    if not current_temperatures:
        return temperatures  # No reference temps; everything qualifies

    min_curr = min(current_temperatures)
    max_curr = max(current_temperatures)

    selected = []
    for t in temperatures:
        # Outside the current range entirely
        if t < min_curr or t > max_curr:
            selected.append(t)
            continue

        # Check if it's far enough from every current temp
        if all(abs(t - ct) > distance for ct in current_temperatures):
            selected.append(t)

    return selected

# This should be called by the test harness with the get_postprocessor.i
# input file to obtain changing postprocessor values from the web server
if __name__ == '__main__':

    # St
    temperatures = [[300, 400]]
    k = [compute_thermal_conductivity(temperatures[0])]

    model = MyTCNet().double()

    train_nn(model, temperatures, k, 10, 1e-2)

    for i in range(3):

        new_temps = run_moose('web_server')

        selected_temps = select_temperatures(new_temps, temperatures[0], 10)

        new_ks = compute_thermal_conductivity(selected_temps)

        temperatures[0] = temperatures[0] + selected_temps
        k[0] = k[0] + new_ks

        train_nn(model, temperatures, k, 10, 1e-2)

        print(temperatures)



