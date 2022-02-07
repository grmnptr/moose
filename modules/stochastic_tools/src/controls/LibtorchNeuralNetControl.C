//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "LibtorchNeuralNetControl.h"
#include "Function.h"
#include "Transient.h"

registerMooseObject("MooseApp", LibtorchNeuralNetControl);

InputParameters
LibtorchNeuralNetControl::validParams()
{
  InputParameters params = Control::validParams();
  params.addClassDescription(
      "Sets the value of a 'Real' input parameter (or postprocessor) based on a Proportional "
      "Integral Derivative control of a postprocessor to match a target a target value.");
  // params.addRequiredParam<PostprocessorName>(
  //     "postprocessor", "The postprocessor to watch for controlling the specified parameter.");
  params.addRequiredParam<FunctionName>("target",
                                        "The target value 1D time function for the postprocessor");
  params.addRequiredParam<std::string>(
      "parameter",
      "The input parameter(s) to control. Specify a single parameter name and all "
      "parameters in all objects matching the name will be updated");
  params.addRequiredParam<PostprocessorName>("postprocessor",
                                             "The postprocessor which stores the control values.");
  // params.addParam<std::string>("parameter_pp",
  //                              "The postprocessor to control. Should be accessed by reference by
  //                              " "the objects depending on its value.");
  // params.addParam<Real>(
  //     "start_time", -std::numeric_limits<Real>::max(), "The time to start the PID controller
  //     at");
  // params.addParam<Real>(
  //     "stop_time", std::numeric_limits<Real>::max(), "The time to stop the PID controller at");
  // params.addParam<bool>(
  //     "reset_every_timestep",
  //     false,
  //     "Reset the PID integral when changing timestep, for coupling iterations within a
  //     timestep");
  // params.addParam<bool>("reset_integral_windup",
  //                       true,
  //                       "Reset the PID integral when the error crosses zero and the integral is "
  //                       "larger than the error.");

  return params;
}

LibtorchNeuralNetControl::LibtorchNeuralNetControl(const InputParameters & parameters)
  : Control(parameters), _target(getFunction("target"))
{
}

void
LibtorchNeuralNetControl::execute()
{
  Point dummy;
  auto value = _target.value(_t, dummy);
  setControllableValue<Real>("parameter", value);
  _fe_problem.setPostprocessorValueByName(getParam<PostprocessorName>("postprocessor"), value);
}
