//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifdef LIBTORCH_ENABLED

#include "LibtorchLiftDragControl.h"
#include "LibtorchTorchScriptNeuralNet.h"
#include "LibtorchArtificialNeuralNet.h"
#include "Transient.h"

registerMooseObject("StochasticToolsApp", LibtorchLiftDragControl);

InputParameters
LibtorchLiftDragControl::validParams()
{
  InputParameters params = LibtorchDRLControl::validParams();
  params.addClassDescription("Bazinga.");
  params.addRequiredParam<Real>("time_delay", "basinga");
  params.addRequiredParam<Real>("ramp", "basinga");

  return params;
}

LibtorchLiftDragControl::LibtorchLiftDragControl(const InputParameters & parameters)
  : LibtorchDRLControl(parameters),
    _time_delay(getParam<Real>("time_delay")),
    _switch_time(_t),
    _ramp(getParam<Real>("ramp"))
{
}

void
LibtorchLiftDragControl::execute()
{
  if (_nn)
  {
    unsigned int n_controls = _control_names.size();
    unsigned int num_old_timesteps = _input_timesteps - 1;

    // Fill a vector with the current values of the responses
    updateCurrentResponse();

    // If this is the first time this control is called and we need to use older values, fill up the
    // needed old values using the initial values
    if (!_initialized)
    {
      _old_responses.clear();
      for (unsigned int step_i = 0; step_i < num_old_timesteps; ++step_i)
        _old_responses.push_back(_current_response);
      _initialized = true;
    }

    // Organize the old an current solution into a tensor so we can evaluate the neural net
    torch::Tensor input_tensor = prepareInputTensor();

    if (_t >= _switch_time + _time_delay || _t_step == 1)
    {
      // Evaluate the neural network to get the expected control value
      torch::Tensor output_tensor = _nn->forward(input_tensor);

      // Sample control value (action) from Gaussian distribution
      _action_tensor = at::normal(output_tensor, _std);

      // Compute log probability
      _log_probability_tensor = computeLogProbability(_action_tensor, output_tensor);

      _switch_time = _t;
    }

    std::vector<Real> action_signal = {_action_tensor.data_ptr<Real>(),
                                       _action_tensor.data_ptr<Real>() + _action_tensor.size(1)};
    // Convert data
    std::transform(_current_control_signals.cbegin(),
                   _current_control_signals.cend(),
                   action_signal.cbegin(),
                   _current_control_signals.begin(),
                   [this](const Real & signal, const Real & action)
                   { return signal + _ramp * (action - signal); });

    _current_control_signal_log_probabilities = {_log_probability_tensor.data_ptr<Real>(),
                                                 _log_probability_tensor.data_ptr<Real>() +
                                                     _log_probability_tensor.size(1)};

    for (unsigned int control_i = 0; control_i < n_controls; ++control_i)
    {
      // We scale the controllable value for physically meaningful control action
      setControllableValueByName<Real>(_control_names[control_i],
                                       _current_control_signals[control_i] *
                                           _action_scaling_factors[control_i]);
    }

    // We add the curent solution to the old solutions and move everything in there one step
    // backward
    std::rotate(_old_responses.rbegin(), _old_responses.rbegin() + 1, _old_responses.rend());
    _old_responses[0] = _current_response;
  }
}

#endif
