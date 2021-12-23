//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ActiveLearningBayesGaussBasis.h"
#include "Sampler.h"
#include "AdaptiveMonteCarloUtils.h"

registerMooseObject("StochasticToolsApp", ActiveLearningBayesGaussBasis);

InputParameters
ActiveLearningBayesGaussBasis::validParams()
{
  InputParameters params = ActiveLearningReporterTempl<Real>::validParams();
  params.addClassDescription("Evaluates parsed function to determine if sample needs to be "
                             "evaluated, otherwise data is set to a default value.");
  params.addRequiredParam<ReporterName>("output_value", "Value of the model output from the SubApp.");
  params.addRequiredParam<SamplerName>("sampler", "The sampler object.");
  return params;
}

ActiveLearningBayesGaussBasis::ActiveLearningBayesGaussBasis(
    const InputParameters & parameters)
  : ActiveLearningReporterTempl<Real>(parameters),
  _output_value(getReporterValue<std::vector<Real>>("output_value")),
  _step(getCheckedPointerParam<FEProblemBase *>("_fe_problem_base")->timeStep()),
  _sampler(getSampler("sampler"))
{
  _inputs_sto.resize(_sampler.getNumberOfCols());
  _decision = true;
}

std::vector<Real> ActiveLearningBayesGaussBasis::Normalize(const std::vector<Real> & inputs, const std::vector<Real> & ref)
{
  Real mu1 = AdaptiveMonteCarloUtils::computeMean(ref, 1);
  Real std1 = AdaptiveMonteCarloUtils::computeSTD(ref, 1);
  // std::cout << mu1 << std::endl;
  // std::cout << std1 << std::endl;
  std::vector<Real> result;
  result.resize(inputs.size());
  for (unsigned int k = 0; k < inputs.size(); ++k)
    result[k] = (inputs[k] - mu1) / std1;
  return result;
  // return inputs;
}

Real ActiveLearningBayesGaussBasis::InvNormalize(const Real & inputs, const std::vector<Real> & ref)
{
  Real mu1 = AdaptiveMonteCarloUtils::computeMean(ref, 1);
  Real std1 = AdaptiveMonteCarloUtils::computeSTD(ref, 1);
  return (inputs * std1 + mu1);
  // return inputs;
}

void ActiveLearningBayesGaussBasis::computeParams(const Real & userSig, const Real & priorSig)
{
  std::vector<std::vector<Real>> inputs_norm;
  inputs_norm.resize(_inputs_sto.size());
  for (unsigned int k = 0; k < _inputs_sto.size(); ++k)
    inputs_norm[k] = Normalize(_inputs_sto[k], _inputs_sto[k]);
  std::vector<Real> outputs_norm;
  outputs_norm.resize(_outputs_sto.size());
  outputs_norm = Normalize(_outputs_sto, _outputs_sto);
  RealEigenMatrix Sig = 1/priorSig * RealEigenMatrix::Identity(_step-1,_step-1);
  RealEigenMatrix psi;
  psi.resize(_step-1,_step-1);
  for (int i = 0; i < _step-1; ++i)
  {
    for (int j = 0; j < _step-1; ++j)
    {
      Real tmp = 0.0;
      for (unsigned int k = 0; k < inputs_norm.size(); ++k)
        tmp += std::pow((inputs_norm[k][j] - inputs_norm[k][i]), 2);
      psi(i,j) = std::exp(tmp/2);
    }
  }
  // std::cout << Moose::stringify(psi) << std::endl;
  RealEigenMatrix SigP_inv = std::pow(userSig, -2) * psi.transpose() * psi + Sig;
  _SigP = SigP_inv.inverse();
  RealEigenVector output_tmp;
  output_tmp.resize(outputs_norm.size());
  for (unsigned int i = 0; i < outputs_norm.size(); ++i)
    output_tmp(i) = outputs_norm[i];
  _muP = std::pow(userSig, -2) * _SigP * psi.transpose() * output_tmp;
}

void ActiveLearningBayesGaussBasis::computePrediction(const std::vector<Real> & inputs, const Real & userSig)
{
  std::vector<std::vector<Real>> inputs_norm;
  inputs_norm.resize(_inputs_sto.size());
  for (unsigned int k = 0; k < _inputs_sto.size(); ++k)
    inputs_norm[k] = Normalize(_inputs_sto[k], _inputs_sto[k]);
  RealEigenVector psi_pred;
  psi_pred.resize(_muP.size());
  for (unsigned int j = 0; j < _muP.size(); ++j)
  {
    Real tmp = 0.0;
    for (unsigned int k = 0; k < inputs.size(); ++k)
      tmp += std::pow((inputs[k] - inputs_norm[k][j]), 2);
    psi_pred(j) = std::exp(tmp/2);
  }
  _muPredict = _muP.transpose() * psi_pred;
  _SigPredict = std::pow((psi_pred.transpose() * _SigP * psi_pred + userSig * userSig), 0.5);
}

bool
ActiveLearningBayesGaussBasis::needSample(const std::vector<Real> & row,
                                              dof_id_type,
                                              dof_id_type,
                                              Real & val)
{
  int N = 12;
  if (_step < N)
  {
    if (_step > 1)
    {
      _outputs_sto.push_back(_output_value[0]);
      for (unsigned int k = 0; k < _inputs_sto.size(); ++k)
        _inputs_sto[k].push_back(_inputs_prev[k]);
    }
    _decision = true;
  } else if (_step == N)
  {
    _outputs_sto.push_back(_output_value[0]);
    for (unsigned int k = 0; k < _inputs_sto.size(); ++k)
      _inputs_sto[k].push_back(_inputs_prev[k]);
    computeParams(1.0, 1000.0);
    // std::cout << "Here inps " << Moose::stringify(_inputs_sto) << std::endl;
    // std::cout << "Here outs " << Moose::stringify(_outputs_sto) << std::endl;
    // for (unsigned int k = 0; k < _muP.size(); ++k)
    //   std::cout << "Here ****** " << Moose::stringify(_muP(k)) << std::endl;
    // std::cout << "Here test inp " << Moose::stringify(row) << std::endl;
    computePrediction(row, 1.0);
    // std::cout << "Mean is " << _muPredict << std::endl;
    // std::cout << "Std is " << _SigPredict << std::endl;
    // std::cout << "InvNorm Mean is " << InvNormalize(_muPredict, _outputs_sto) << std::endl;
    // computePrediction(row, 1.0);
    if (_SigPredict  <= 0.01) // / _muPredict
    {
      // _decision = true;
      _decision = false;
      val = InvNormalize(_muPredict, _outputs_sto);
    } else
      _decision = true;
  } else
  {
    if (_decision == true)
    {
      _outputs_sto.push_back(_output_value[0]);
      for (unsigned int k = 0; k < _inputs_sto.size(); ++k)
        _inputs_sto[k].push_back(_inputs_prev[k]);
      computeParams(1.0, 1000.0);
    }
    computePrediction(row, 1.0);
    std::cout << "Mean is " << _muPredict << std::endl;
    std::cout << "Std is " << _SigPredict << std::endl;
    std::cout << "InvNorm Mean is " << InvNormalize(_muPredict, _outputs_sto) << std::endl;

    if (_SigPredict  <= 0.01) // / _muPredict
    {
      // _decision = true;
      _decision = false;
      val = InvNormalize(_muPredict, _outputs_sto);
    } else
      _decision = true;
  }
  _inputs_prev = row;
  return _decision;
}
