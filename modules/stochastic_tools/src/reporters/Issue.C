//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "Issue.h"
#include "Sampler.h"

registerMooseObject("StochasticToolsApp", Issue);

template <typename T>
InputParameters
IssueTempl<T>::validParams()
{
  InputParameters params = ActiveLearningReporterTempl<T>::validParams();
  params += FunctionParserUtils<false>::validParams();
  params.addClassDescription("Evaluates parsed function to determine if sample needs to be "
                             "evaluated, otherwise data is set to a default value.");

  params.addRequiredCustomTypeParam<std::string>(
      "function",
      "FunctionExpression",
      "function expression that should evaluate to 0 if sample should NOT be evaluated.");
  params.addRequiredParam<T>("default_value",
                             "Value to replace in data when sample is not evaluated.");
  params.addParam<std::vector<std::string>>(
      "sampler_vars", std::vector<std::string>(), "Vector of variables defined by sampler values.");
  params.addParam<std::vector<unsigned int>>("sampler_var_indices",
                                             std::vector<unsigned int>(),
                                             "Vector of indices defining the sampler column index "
                                             "associated with the variables in 'sampler_vars'.");
  params.addParam<bool>(
      "use_time", true, "Make time (t) variable available in the function expression.");
  return params;
}

template <typename T>
IssueTempl<T>::IssueTempl(const InputParameters & parameters)
  : ActiveLearningReporterTempl<T>(parameters),
    FunctionParserUtils<false>(parameters),
    _function(this->template getParam<std::string>("function")),
    _default_value(this->template getParam<T>("default_value")),
    _sampler_vars(
        this->template getParam<std::string, unsigned int>("sampler_vars", "sampler_var_indices")),
    _use_time(this->template getParam<bool>("use_time")),
    _step(this->template getCheckedPointerParam<FEProblemBase *>("_fe_problem_base")->timeStep())
{
  // Gather variables and see if any of the indexes are out-of-bounds
  std::stringstream vars;
  std::string sep = "";
  for (const auto & it : _sampler_vars)
  {
    vars << sep << it.first;
    if (it.second >= this->sampler().getNumberOfCols())
      this->paramError("sampler_var_indices",
                       "Provided index ",
                       it.second,
                       " must be smaller than the sampler number of columns (",
                       this->sampler().getNumberOfCols(),
                       ").");
    sep = ",";
  }
  if (_use_time)
    vars << sep << "t";

  // Parse the inputted function
  _func_F = std::make_shared<SymFunction>();
  setParserFeatureFlags(_func_F);
  if (_func_F->Parse(_function, vars.str()) >= 0)
    this->paramError("function", "Invalid function:\n", _function, "\n", _func_F->ErrorMsg());
  if (!_disable_fpoptimizer)
    _func_F->Optimize();
  if (_enable_jit)
    _func_F->JITCompile();
  _func_params.resize(_sampler_vars.size() + (_use_time ? 1 : 0));
}

template <typename T>
bool
IssueTempl<T>::needSample(const std::vector<Real> & row,
                          dof_id_type local_ind,
                          dof_id_type,
                          T & val)
{
  bool decision = false;
  std::cout << "Step " << _step << std::endl;
  std::cout << "Local index: " << local_ind << std::endl;

  if (_step == 0)
    val = _default_value;

  if (_step == 1)
  {
    if (local_ind == 0)
      decision = true;
    // else if (local_ind == 1)
    //   decision = true;
    else
      val = _default_value;
  }

  if (_step == 2)
  {
    if (local_ind == 1)
      decision = true;
    // else if (local_ind == 2)
    //   decision = true;
    else
      val = _default_value;
  }

  if (_step == 3)
  {
    if (local_ind == 0)
      decision = true;
    // else if (local_ind == 2)
    //   decision = true;
    else
      val = _default_value;
  }

  return decision;
}

// Explicit instantiation (more types can easily be added)
template class IssueTempl<Real>;
