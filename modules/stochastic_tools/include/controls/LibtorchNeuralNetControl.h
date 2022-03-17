//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "LibtorchSimpleNeuralNet.h"
#include "Control.h"

/**
 * A time-dependent, neural network-based control of multiple input parameters.
 * The control strategy depends on the training of the neural net, which is
 * typically done in a trainer object in the master app.
 */
class LibtorchNeuralNetControl : public Control
{
public:
  static InputParameters validParams();

  /**
   * Class constructor
   * @param parameters Input parameters for this Control object
   */
  LibtorchNeuralNetControl(const InputParameters & parameters);

  virtual void execute() override;

#ifdef TORCH_ENABLED
  void
  loadControlNeuralNet(const std::shared_ptr<StochasticTools::LibtorchSimpleNeuralNet> & input_nn);
#endif

private:
  std::vector<Real> _old_response;
  std::vector<Real> _current_response;

  bool _initialized;

  std::vector<std::string> _control_names;

  std::vector<PostprocessorName> _response_names, _postprocessor_names;

#ifdef TORCH_ENABLED
  /// Pointer to the neural net object which is supposed to be used to control
  /// the input values
  std::shared_ptr<StochasticTools::LibtorchSimpleNeuralNet> _nn;
#endif
};
