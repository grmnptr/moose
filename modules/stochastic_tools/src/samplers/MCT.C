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
  // params.addRequiredParam<ReporterName>("flag_sample", "Flag samples.");
  return params;
}

MCT::MCT(const InputParameters & parameters)
  : Sampler(parameters),
    ReporterInterface(this),
    _distribution_names(getParam<std::vector<DistributionName>>("distributions")),
    // _flag_sample(getReporterValue<std::vector<bool>>("flag_sample")),
    _step(getCheckedPointerParam<FEProblemBase *>("_fe_problem_base")->timeStep())
{
  for (const DistributionName & name : _distribution_names)
    _distributions.push_back(&getDistributionByName(name));

  setNumberOfRows(getParam<dof_id_type>("num_rows"));
  setNumberOfCols(_distributions.size());
  setNumberOfRandomSeeds(100000);
  _inputs_sto.resize(getParam<dof_id_type>("num_rows"));
  for (unsigned int i = 0; i < _inputs_sto.size(); ++i)
    _inputs_sto[i].resize(_distributions.size());
  _check_step = 0;
}

Real
MCT::computeSample(dof_id_type row_index, dof_id_type col_index)
{
  if (_step == 1)
  {
    if (row_index == 0)
    {
      if (col_index == 0)
        _inputs_sto[row_index][col_index] = 0.05;
      else
        _inputs_sto[row_index][col_index] = 0.95;
    }

    if (row_index == 1)
    {
      if (col_index == 0)
        _inputs_sto[row_index][col_index] = 0.25;
      else
        _inputs_sto[row_index][col_index] = 0.75;
    }

    // if (row_index == 2)
    // {
    //   if (col_index == 0)
    //     _inputs_sto[row_index][col_index] = 0.15;
    //   else
    //     _inputs_sto[row_index][col_index] = 0.55;
    // }
  }

  if (_step == 2)
  {
    if (row_index == 0)
    {
      if (col_index == 0)
        _inputs_sto[row_index][col_index] = 0.0;
      else
        _inputs_sto[row_index][col_index] = 0.0;
    }

    if (row_index == 1)
    {
      if (col_index == 0)
        _inputs_sto[row_index][col_index] = 0.30;
      else
        _inputs_sto[row_index][col_index] = 0.70;
    }

    // if (row_index == 2)
    // {
    //   if (col_index == 0)
    //     _inputs_sto[row_index][col_index] = 0.35;
    //   else
    //     _inputs_sto[row_index][col_index] = 0.85;
    // }
  }

  if (_step == 3)
  {
    if (row_index == 0)
    {
      if (col_index == 0)
        _inputs_sto[row_index][col_index] = 0.1;
      else
        _inputs_sto[row_index][col_index] = 0.1;
    }

    if (row_index == 1)
    {
      if (col_index == 0)
        _inputs_sto[row_index][col_index] = 0.25;
      else
        _inputs_sto[row_index][col_index] = 0.45;
    }

    // if (row_index == 2)
    // {
    //   if (col_index == 0)
    //     _inputs_sto[row_index][col_index] = 0.35;
    //   else
    //     _inputs_sto[row_index][col_index] = 0.55;
    // }
  }

  std::cout << " In sampler: step " << _step << " row " << row_index << " Params "
            << Moose::stringify(_inputs_sto) << std::endl;

  if (_step > 3 && col_index == 0 && _check_step != _step) //
  {
    for (dof_id_type j = 0; j < _distributions.size(); ++j)
      _inputs_sto[row_index][j] = _distributions[j]->quantile(getRand(_step));
  }
  _check_step = _step;
  return _inputs_sto[row_index][col_index];
}
