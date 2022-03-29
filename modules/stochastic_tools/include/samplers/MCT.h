//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Sampler.h"

/**
 * A class used to perform Monte Carlo Sampling
 */
class MCT : public Sampler, public ReporterInterface
{
public:
  static InputParameters validParams();

  MCT(const InputParameters & parameters);

protected:
  /// Return the sample for the given row and column
  virtual Real computeSample(dof_id_type row_index, dof_id_type col_index) override;

  /// Storage for distribution objects to be utilized
  std::vector<Distribution const *> _distributions;

  /// Distribution names
  const std::vector<DistributionName> & _distribution_names;

  const bool & _flag_sample;
  bool _last_flag;

private:
  /// Track the current step of the main App
  const int & _step;

  // /// Ensure that the MCMC algorithm proceeds in a sequential fashion
  // int _check_step;
  //
  // /// Storage for previously accepted samples by the decision reporter system
  std::vector<Real> _inputs_sto;
};
