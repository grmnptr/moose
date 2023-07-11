//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "LiftDragDRLRewardFunction.h"

registerMooseObject("StochasticToolsApp", LiftDragDRLRewardFunction);

InputParameters
LiftDragDRLRewardFunction::validParams()
{
  InputParameters params = Function::validParams();

  params.addClassDescription(
      "Implements R = -<C_drag>-0.2*|<C_lift>| where <> is a time integral over a given interval.");

  params.addRequiredParam<PostprocessorName>(
      "lift_pp", "The name of the postprocessor which computes the lift coefficient.");
  params.addRequiredParam<PostprocessorName>(
      "drag_pp", "The name of the postprocessor which computes the drag coefficient.");

  params.addRequiredParam<unsigned int>("smoothing_interval", "The averaging interval.");

  return params;
}

LiftDragDRLRewardFunction::LiftDragDRLRewardFunction(const InputParameters & parameters)
  : Function(parameters),
    _observed_drag(getPostprocessorValueByName(getParam<PostprocessorName>("drag_pp"))),
    _observed_lift(getPostprocessorValueByName(getParam<PostprocessorName>("lift_pp"))),
    _avg_interval(getParam<unsigned int>("smoothing_interval")),
    _c2(getParam<Real>("c2"))
{
}

Real
LiftDragDRLRewardFunction::value(Real /*t*/, const Point & /*p*/) const
{
  _drag_history.push_back(_observed_drag);
  _lift_history.push_back(_observed_lift);

  auto num_history(std::min(_avg_interval, _drag_history.size()));

  Real sum_drag = std::transform_reduce(_drag_history.cend() - num_history,
                                        _drag_history.cend(),
                                        0.0,
                                        std::plus{},
                                        [](auto val) { return val; });
  Real mean_drag = sum_drag / num_history;

  Real sum_lift = std::transform_reduce(_lift_history.cend() - num_history,
                                        _lift_history.cend(),
                                        0.0,
                                        std::plus{},
                                        [](auto val) { return std::abs(val); });
  Real mean_lift = sum_lift / num_history;

  return -mean_drag - 0.2 * mean_lift;
}

ADReal
LiftDragDRLRewardFunction::value(const ADReal & t, const ADPoint & p) const
{
  _drag_history.push_back(_observed_drag);
  _lift_history.push_back(_observed_lift);

  auto num_history(std::min(_avg_interval, _drag_history.size()));

  Real sum_drag = std::transform_reduce(_drag_history.cend() - num_history,
                                        _drag_history.cend(),
                                        0.0,
                                        std::plus{},
                                        [](auto val) { return val; });
  Real mean_drag = sum_drag / num_history;

  Real sum_lift = std::transform_reduce(_lift_history.cend() - num_history,
                                        _lift_history.cend(),
                                        0.0,
                                        std::plus{},
                                        [](auto val) { return std::abs(val); });
  Real mean_lift = sum_lift / num_history;

  return -mean_drag - 0.2 * mean_lift;
}
