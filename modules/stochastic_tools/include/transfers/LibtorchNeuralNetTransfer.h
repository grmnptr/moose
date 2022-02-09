//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

// files that contain torch includes
#include "LibtorchSimpleNNControlTrainer.h"

// MOOSE includes
#include "MultiAppTransfer.h"
#include "SurrogateModelInterface.h"

class LibtorchNeuralNetTransfer : public MultiAppTransfer, SurrogateModelInterface
{
public:
  static InputParameters validParams();

  LibtorchNeuralNetTransfer(const InputParameters & parameters);

  virtual void execute() override;

protected:
  /**
   * The trainer object which contains the information about the neural net.
   */

  const std::string _control_name;

  const LibtorchSimpleNNControlTrainer & _trainer;
};
