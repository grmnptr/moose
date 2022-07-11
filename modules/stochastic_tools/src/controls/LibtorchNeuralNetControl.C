//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "LibtorchNeuralNetControl.h"
#include "Transient.h"

registerMooseObject("MooseApp", LibtorchNeuralNetControl);

InputParameters
LibtorchNeuralNetControl::validParams()
{
  InputParameters params = Control::validParams();
  params.addClassDescription(
      "Sets the value of a 'Real' input parameter (or postprocessor) based on a Proportional "
      "Integral Derivative control of a postprocessor to match a target a target value.");
  params.addRequiredParam<std::vector<std::string>>("parameters",
                                                    "The input parameter(s) to control.");
  params.addRequiredParam<std::vector<PostprocessorName>>(
      "responses", "The responses (prostprocessors) which are used for the control.");
  params.addRequiredParam<std::vector<PostprocessorName>>(
      "postprocessors", "The postprocessors which stores the control values.");

  return params;
}

LibtorchNeuralNetControl::LibtorchNeuralNetControl(const InputParameters & parameters)
  : Control(parameters),
    _initialized(false),
    _control_names(getParam<std::vector<std::string>>("parameters")),
    _response_names(getParam<std::vector<PostprocessorName>>("responses")),
    _postprocessor_names(getParam<std::vector<PostprocessorName>>("postprocessors"))
{
}

void
LibtorchNeuralNetControl::execute()
{
#ifdef LIBTORCH_ENABLED
  if (_nn != NULL)
  {
    unsigned int n_responses = _response_names.size();

    _current_response.clear();
    for (unsigned int resp_i = 0; resp_i < _response_names.size(); ++resp_i)
      _current_response.push_back(getPostprocessorValueByName(_response_names[resp_i]));

    if (!_initialized)
    {
      _old_response = _current_response;
      _initialized = true;
    }

    std::vector<Real> raw_input(_old_response);
    raw_input.insert(raw_input.end(), _current_response.begin(), _current_response.end());

    auto options = torch::TensorOptions().dtype(at::kDouble);
    torch::Tensor input_tensor =
        torch::from_blob(raw_input.data(), {1, 2 * n_responses}, options).to(at::kDouble);

    torch::Tensor output_tensor = _nn->forward(input_tensor);

    std::vector<Real> converted_output = {output_tensor.data_ptr<Real>(),
                                          output_tensor.data_ptr<Real>() + output_tensor.size(0)};

    for (unsigned int control_i = 0; control_i < _control_names.size(); ++control_i)
    {
      setControllableValueByName<Real>(_control_names[control_i], converted_output[control_i]);
      _fe_problem.setPostprocessorValueByName(_postprocessor_names[control_i],
                                              converted_output[control_i]);
    }

    _old_response = _current_response;
  }
#endif
}

#ifdef LIBTORCH_ENABLED
void
LibtorchNeuralNetControl::loadControlNeuralNet(
    const std::shared_ptr<Moose::LibtorchArtificialNeuralNet> & input_nn)
{
  std::vector<std::string> activation_names;
  const MultiMooseEnum & activation_functions = input_nn->activationFunctions();
  for (unsigned int i = 0; i < activation_functions.size(); ++i)
    activation_names.push_back(activation_functions[i]);

  _nn = std::make_shared<Moose::LibtorchArtificialNeuralNet>(input_nn->name(),
                                                             input_nn->numInputs(),
                                                             input_nn->numOutputs(),
                                                             input_nn->numNeuronsPerLayer(),
                                                             activation_names);

  torch::load(_nn, input_nn->name());
}
#endif
