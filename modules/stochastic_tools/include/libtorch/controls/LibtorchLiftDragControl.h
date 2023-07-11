#ifdef LIBTORCH_ENABLED

#pragma once

#include "LibtorchArtificialNeuralNet.h"
#include "LibtorchNeuralNetControl.h"
#include "LibtorchDRLControl.h"

class LibtorchLiftDragControl : public LibtorchDRLControl
{
public:
  static InputParameters validParams();

  LibtorchLiftDragControl(const InputParameters & parameters);

  virtual void execute() override;

protected:
  std::vector<Real> _current_valid_control_values;

  const Real _time_delay;
  Real _switch_time;
  const Real _ramp;

  torch::Tensor _action_tensor;
  torch::Tensor _log_probability_tensor;
};

#endif
