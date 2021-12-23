//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ActiveLearningReporterBase.h"
// #include "GeneralReporter.h"
#include <Eigen/Dense>

class ActiveLearningBayesGaussBasis : public ActiveLearningReporterTempl<Real>
{
public:
  static InputParameters validParams();
  ActiveLearningBayesGaussBasis(const InputParameters & parameters);

protected:

  std::vector<Real> Normalize(const std::vector<Real> & inputs, const std::vector<Real> & ref);

  Real InvNormalize(const Real & inputs, const std::vector<Real> & ref);

  void computeParams(const Real & userSig, const Real & priorSig);

  void computePrediction(const std::vector<Real> & inputs, const Real & userSig);

  /**
   * This evaluates the inputted function to determine whether a multiapp solve is
   * necessary/allowed, otherwise it replaces the "transferred" quantity with a
   * default value.
   */
  virtual bool needSample(const std::vector<Real> & row,
                          dof_id_type local_ind,
                          dof_id_type global_ind,
                          Real & val) override;

  /// Model output value from SubApp
  const std::vector<Real> & _output_value;

private:
  /// Track the current step of the main App
  const int & _step;

  /// The adaptive Monte Carlo sampler
  Sampler & _sampler;

  std::vector<std::vector<Real>> _inputs_sto;

  std::vector<Real> _outputs_sto;

  RealEigenMatrix _SigP;

  RealEigenVector _muP;

  Real _muPredict;

  Real _SigPredict;

  /// store the default value
  // Real _default_value;

  std::vector<Real> _inputs_prev;

  bool _decision;
};
