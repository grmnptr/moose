# LibtorchDRLControl

!syntax description /Controls/LibtorchDRLControl

## Overview

This object controls a physical process using a neural network, just like [LibtorchNeuralNetControl.md], 
with an additional functionality of randomizing the action values to avoid overfitting in the control process.
This control object is supposed to be used in conjunction with [LibtorchDRLControlTrainer.md]. In other 
cases when the neural network needs to be simply evaluated, the user is encouraged to use [LibtorchNeuralNetControl.md]. 

!syntax parameters /Controls/LibtorchDRLControl

!syntax inputs /Controls/LibtorchDRLControl

!syntax children /Controls/LibtorchDRLControl