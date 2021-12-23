//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "MCT.h"
#include "Distribution.h"

registerMooseObject("StochasticToolsApp", MCT);

InputParameters
MCT::validParams()
{
  InputParameters params = Sampler::validParams();
  params.addClassDescription("Monte Carlo Sampler.");
  params.addRequiredParam<dof_id_type>("num_rows", "The number of rows per matrix to generate.");
  params.addRequiredParam<std::vector<DistributionName>>(
      "distributions",
      "The distribution names to be sampled, the number of distributions provided defines the "
      "number of columns per matrix.");
  return params;
}

MCT::MCT(const InputParameters & parameters)
  : Sampler(parameters),
    _distribution_names(getParam<std::vector<DistributionName>>("distributions")),
    _step(getCheckedPointerParam<FEProblemBase *>("_fe_problem_base")->timeStep())
{
  for (const DistributionName & name : _distribution_names)
    _distributions.push_back(&getDistributionByName(name));

  setNumberOfRows(getParam<dof_id_type>("num_rows"));
  setNumberOfCols(_distributions.size());
  setNumberOfRandomSeeds(100000);
  // _inputs_sto.resize(_distributions.size());
  // _check_step = 0;
}

Real
MCT::computeSample(dof_id_type /*row_index*/, dof_id_type col_index)
{
  // const bool sample = _step > 1 && col_index == 0 && _check_step != _step;
  // _check_step = _step;
  // return _inputs_sto[col_index];
  return _distributions[col_index]->quantile(getRand(_step));
}
