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
  unsigned int n_responses = _response_names.size();
  unsigned int n_controls = _control_names.size();

  _current_response.clear();
  for (unsigned int resp_i = 0; resp_i < _response_names.size(); ++resp_i)
    _current_response.push_back(getPostprocessorValueByName(_respons_names[resp_i]));

  if (!_initialized)
  {
    _old_response = _current_response;
    _initialized = true;
  }

  raw_input.insert(raw_input.end(), _current_response.begin(), _current_response.end());

  auto options = torch::TensorOptions().dtype(at::kDouble);
  torch::Tensor input_tensor =
      torch::from_blob(raw_input.data(), {1, 2 * n_responses}, options).to(at::kDouble);

  torch::Tensor output_tensor = _nn->forward(input_tensor);

  std::vector<Real> converted_output = {output_tensor.data_ptr<Real>(),
                                        output_tensor.data_ptr<Real>() + output_tensor.size(1)};

  for (unsigned int control_i = 0; control_i < _control_names.size(); ++control_i)
  {
    setControllableValueByName<Real>(_control_names[control_i], converted_output[control_i]);
    _fe_problem.setPostprocessorValueByName(_postprocessor_names[control_i],
                                            converted_output[control_i]);
  }

  _old_response = _current_response;
}
