//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "AL_ADAM.h"
#include "Sampler.h"
#include "AdaptiveMonteCarloUtils.h"

#include <petsctao.h>
#include <petscdmda.h>

#include "libmesh/petsc_vector.h"
#include "libmesh/petsc_matrix.h"

#include <math.h>

registerMooseObject("StochasticToolsApp", AL_ADAM);

InputParameters
AL_ADAM::validParams()
{
  InputParameters params = ActiveLearningReporterTempl<Real>::validParams();
  params += CovarianceInterface::validParams();
  params.addClassDescription("Evaluates parsed function to determine if sample needs to be "
                             "evaluated, otherwise data is set to a default value.");
  // params.addRequiredParam<ReporterName>("output_value", "Value of the model output from the SubApp.");
  params.addRequiredParam<SamplerName>("sampler", "The sampler object.");
  params.addRequiredParam<UserObjectName>("covariance_function", "Name of covariance function.");
  params.addParam<ReporterValueName>("flag_sample", "flag_sample", "Flag samples.");
  params.addParam<bool>(
      "standardize_params", true, "Standardize (center and scale) training parameters (x values)");
  params.addParam<bool>(
      "standardize_data", true, "Standardize (center and scale) training data (y values)");
  params.addParam<std::string>(
      "tao_options", "", "Command line options for PETSc/TAO hyperparameter optimization");
  params.addParam<bool>("show_tao", false, "Switch to show TAO solver results");
  params.addParam<std::vector<std::string>>("tune_parameters",
                                            "Select hyperparameters to be tuned");
  params.addParam<std::vector<Real>>("tuning_min", "Minimum allowable tuning value");
  params.addParam<std::vector<Real>>("tuning_max", "Maximum allowable tuning value");
  return params;
}

AL_ADAM::AL_ADAM(
    const InputParameters & parameters)
  : ActiveLearningReporterTempl<Real>(parameters),
  CovarianceInterface(parameters),
  // _output_value(getReporterValue<std::vector<Real>>("output_value")),
  _step(getCheckedPointerParam<FEProblemBase *>("_fe_problem_base")->timeStep()),
  _sampler(getSampler("sampler")),
  _param_standardizer(declareModelData<StochasticTools::Standardizer>("_param_standardizer")),
  _data_standardizer(declareModelData<StochasticTools::Standardizer>("_data_standardizer")),
  _covariance_function(
      getCovarianceFunctionByName(getParam<UserObjectName>("covariance_function"))),
  _flag_sample(declareValue<bool>("flag_sample")),
  _do_tuning(isParamValid("tune_parameters")),
  _tao_options(getParam<std::string>("tao_options")),
  _show_tao(getParam<bool>("show_tao")),
  _tao_comm(MPI_COMM_SELF),
  _covar_type(declareModelData<std::string>("_covar_type", _covariance_function->type())),
  _training_params(declareModelData<RealEigenMatrix>("_training_params")),
  _K(declareModelData<RealEigenMatrix>("_K")),
  _K_results_solve(declareModelData<RealEigenMatrix>("_K_results_solve")),
  // _K_cho_decomp(declareModelData<Eigen::LLT<RealEigenMatrix>>("_K_cho_decomp")),
  _standardize_params(getParam<bool>("standardize_params")),
  _standardize_data(getParam<bool>("standardize_data")),
  _hyperparam_map(declareModelData<std::unordered_map<std::string, Real>>("_hyperparam_map")),
  _hyperparam_vec_map(declareModelData<std::unordered_map<std::string, std::vector<Real>>>(
      "_hyperparam_vec_map"))
{
  _inputs_sto.resize(_sampler.getNumberOfCols());
  _decision = true;
  _covar_type = _covariance_function->type();

  _amp_sto = _covariance_function->getSignalVariance();
  _len_sto = _covariance_function->getLengthFactor();

  _num_tunable = 0;
  std::vector<std::string> tune_parameters(getParam<std::vector<std::string>>("tune_parameters"));
  for (unsigned int ii = 0; ii < tune_parameters.size(); ++ii)
  {
    const auto & hp = tune_parameters[ii];
    if (_covariance_function->isTunable(hp))
    {
      unsigned int size;
      Real min;
      Real max;
      // Get size and default min/max
      _covariance_function->getTuningData(hp, size, min, max);
      // Check for overridden min/max
      min = isParamValid("tuning_min") ? getParam<std::vector<Real>>("tuning_min")[ii] : min;
      max = isParamValid("tuning_max") ? getParam<std::vector<Real>>("tuning_max")[ii] : max;
      // Save data in tuple
      _tuning_data[hp] = std::make_tuple(_num_tunable, size, min, max);
      _num_tunable += size;
    }
  }
  _flag_sample = false;
}

void
AL_ADAM::SetupData(const std::vector<std::vector<Real>> & inputs, const std::vector<Real> & outputs)
{
  _training_params.setZero(outputs.size(), inputs.size());
  _training_data.setZero(outputs.size(), 1);

  for (unsigned int i = 0; i < outputs.size(); ++i)
  {
    _training_data(i, 0) = outputs[i];
    for (unsigned int j = 0; j < inputs.size(); ++j)
      _training_params(i,j) = inputs[j][i];
  }
}

void
AL_ADAM::Train_ADAM(const unsigned int & iters)
{
  SetupData(_inputs_sto, _outputs_sto);
  // std::cout << Moose::stringify(_training_params) << std::endl;
  // std::cout << Moose::stringify(_training_data) << std::endl;
  for (unsigned int ii = 0; ii < _training_params.rows(); ++ii)
  {
    for (unsigned int jj = 0; jj < _training_params.cols(); ++jj)
      gatherSum(_training_params(ii, jj));
    gatherSum(_training_data(ii, 0));
  }

  // Standardize (center and scale) training params
  if (_standardize_params)
  {
    _param_standardizer.computeSet(_training_params);
    _param_standardizer.getStandardized(_training_params);
  }
  // if not standardizing data set mean=0, std=1 for use in surrogate
  else
    _param_standardizer.set(0, 1, _sampler.getNumberOfCols());

  // Standardize (center and scale) training data
  if (_standardize_data)
  {
    _data_standardizer.computeSet(_training_data);
    _data_standardizer.getStandardized(_training_data);
  }
  // if not standardizing data set mean=0, std=1 for use in surrogate
  else
    _param_standardizer.set(0, 1, _sampler.getNumberOfCols());

  _K.resize(_training_params.rows(), _training_params.rows());

  libMesh::PetscVector<Number> theta(_tao_comm, _num_tunable);
  _covariance_function->buildHyperParamMap(_hyperparam_map, _hyperparam_vec_map);
  mapToVec(theta);
  // std::cout << theta(0) << std::endl;
  // std::cout << theta(1) << std::endl;
  // std::cout << theta(2) << std::endl;
  // unsigned int Iters = 10000;
  std::vector<Real> grad1;
  Real b1;
  Real b2;
  Real eps;
  Real alpha;
  b1 = 0.9; b2 = 0.999; alpha = 0.01; eps = 1e-7;
  std::vector<Real> m0;
  std::vector<Real> v0;
  for (unsigned int ii = 0; ii < _num_tunable; ++ii)
  {
    m0.push_back(0.0);
    v0.push_back(0.0);
  }

  // for (unsigned int ii = 0; ii < _num_tunable; ++ii)
  //   theta.set(ii,1.0);
  // vecToMap(theta);
  // _covariance_function->loadHyperParamMap(_hyperparam_map, _hyperparam_vec_map);

  Real new_val;
  Real m_hat;
  Real v_hat;
  std::cout << "Starting loss " << getLoss() << std::endl;
  Real diff;
  diff = 100;
  Real prev;
  prev = getLoss();
  int jj;
  jj = 0;
  while (diff > 1e-3) // for (unsigned int jj = 1; jj <= iters; ++jj) //
  {
    ++jj;
    grad1 = getGradient();
    // std::cout << Moose::stringify(grad1) << std::endl;
    for (unsigned int ii = 0; ii < _num_tunable; ++ii)
    {
      m0[ii] = b1 * m0[ii] + (1 - b1) * grad1[ii];
      v0[ii] = b2 * v0[ii] + (1 - b2) * grad1[ii] * grad1[ii];
      m_hat = m0[ii] / (1 - std::pow(b1, jj));
      v_hat = v0[ii] / (1 - std::pow(b2, jj));
      new_val = theta(ii) - alpha * m_hat / (std::pow(v_hat, 0.5) + eps);
      if (new_val < 0.0)
        new_val = 0.001;
      theta.set(ii,new_val);
    }
    // std::cout << Moose::stringify(m0) << std::endl;
    // std::cout << Moose::stringify(v0) << std::endl;
    // theta.print();
    vecToMap(theta);
    _covariance_function->loadHyperParamMap(_hyperparam_map, _hyperparam_vec_map);
    diff = std::abs(prev - getLoss()); // / (prev);
    prev = getLoss();
  }
  std::cout << "Ending loss " << getLoss() << std::endl;
  theta.print();
  // std::cout << getLoss() << std::endl;
  // std::cout << Moose::stringify(getGradient()) << std::endl;

  // theta.set(0,100.0);
  // vecToMap(theta);
  // _covariance_function->loadHyperParamMap(_hyperparam_map, _hyperparam_vec_map);
  // std::cout << _covariance_function->getSignalVariance() << std::endl;




  // _covariance_function->computeCovarianceMatrix(_K, _training_params, _training_params, true);
  // _K_cho_decomp = _K.llt();
  // _K_results_solve = _K_cho_decomp.solve(_training_data);
  //
  // _covariance_function->buildHyperParamMap(_hyperparam_map, _hyperparam_vec_map);
}

Real
AL_ADAM::getLoss()
{
  // libMesh::PetscVector<Number> theta(theta_vec, _tao_comm);
  // vecToMap(theta);
  // _covariance_function->loadHyperParamMap(_hyperparam_map, _hyperparam_vec_map);
  _covariance_function->computeCovarianceMatrix(_K, _training_params, _training_params, false);
  _K_cho_decomp = _K.llt();
  _K_results_solve = _K_cho_decomp.solve(_training_data);
  Real log_likelihood = 0;
  log_likelihood += -(_training_data.transpose() * _K_results_solve)(0, 0);
  log_likelihood += -std::log(_K.determinant());

  // log_likelihood += -std::pow(((std::log(_covariance_function->getSignalVariance()) - 0.0) / 1.0) , 2) - 2 * std::log(_covariance_function->getSignalVariance());
  // std::vector<Real> len1;
  // len1 = _covariance_function->getLengthFactor();
  // for (unsigned int ii = 0; ii < len1.size(); ++ii)
  //   log_likelihood += -std::pow(((std::log(len1[ii]) - 0.0) / 1.0) , 2) - 2 * std::log(len1[ii]);

  log_likelihood += -_training_data.rows() * std::log(2 * M_PI);
  log_likelihood = -log_likelihood / 2;
  return log_likelihood;
}

std::vector<Real>
AL_ADAM::getGradient()
{
  // std::cout << _covariance_function->getSignalVariance() << std::endl;
  // _covariance_function->loadHyperParamMap(_hyperparam_map, _hyperparam_vec_map);
  // _covariance_function->computeCovarianceMatrix(_K, _training_params, _training_params, false);
  // _K_cho_decomp = _K.llt();
  // _K_results_solve = _K_cho_decomp.solve(_training_data);

  // testing auto tuning
  RealEigenMatrix dKdhp(_training_params.rows(), _training_params.rows());
  RealEigenMatrix alpha = _K_results_solve * _K_results_solve.transpose();
  std::vector<Real> grad_vec;
  grad_vec.resize(_num_tunable);
  int count;
  count = 1;
  // Real sto1;
  for (auto iter = _tuning_data.begin(); iter != _tuning_data.end(); ++iter)
  {
    std::string hyper_param_name = iter->first;
    for (unsigned int ii = 0; ii < std::get<1>(iter->second); ++ii)
    {
      // std::cout << hyper_param_name << std::endl;
      _covariance_function->computedKdhyper(dKdhp, _training_params, hyper_param_name, ii);
      RealEigenMatrix tmp = alpha * dKdhp - _K_cho_decomp.solve(dKdhp);
      Real grad1 = -tmp.trace() / 2.0;
      // if (hyper_param_name.compare("length_factor") == 0)
      // {
      //   // std::cout << _covariance_function->getLengthFactor()[ii] << std::endl;
      //   grad1 += (std::log(_covariance_function->getLengthFactor()[ii])-0.0) * 1/_covariance_function->getLengthFactor()[ii] + 1/_covariance_function->getLengthFactor()[ii];
      // } else
      // {
      //   grad1 += (std::log(_covariance_function->getSignalVariance())-0.0) * 1/_covariance_function->getSignalVariance() + 1/_covariance_function->getSignalVariance();
      // }
      // std::cout << hyper_param_name << std::endl;
      // grad.set(std::get<0>(iter->second) + ii, grad1);
      if (hyper_param_name.compare("length_factor") == 0)
      {
        grad_vec[count] = grad1;
        ++count;
      } else
      {
        grad_vec[0] = grad1;
      }
      // std::cout << hyper_param_name << std::endl;
      // grad_vec.push_back(grad1);
    }
  }
  return grad_vec;
}

void
AL_ADAM::mapToVec(libMesh::PetscVector<Number> & theta) const
{
  for (auto iter = _tuning_data.begin(); iter != _tuning_data.end(); ++iter)
  {
    std::string hyper_param_name = iter->first;
    if (_hyperparam_map.find(hyper_param_name) != _hyperparam_map.end())
    {
      theta.set(std::get<0>(iter->second), _hyperparam_map[hyper_param_name]);
    }
    else if (_hyperparam_vec_map.find(hyper_param_name) != _hyperparam_vec_map.end())
    {
      for (unsigned int ii = 0; ii < std::get<1>(iter->second); ++ii)
      {
        theta.set(std::get<0>(iter->second) + ii, _hyperparam_vec_map[hyper_param_name][ii]);
      }
    }
  }
}

void
AL_ADAM::vecToMap(libMesh::PetscVector<Number> & theta)
{
  for (auto iter = _tuning_data.begin(); iter != _tuning_data.end(); ++iter)
  {
    std::string hyper_param_name = iter->first;
    if (_hyperparam_map.find(hyper_param_name) != _hyperparam_map.end())
    {
      _hyperparam_map[hyper_param_name] = theta(std::get<0>(iter->second));
    }
    else if (_hyperparam_vec_map.find(hyper_param_name) != _hyperparam_vec_map.end())
    {
      for (unsigned int ii = 0; ii < std::get<1>(iter->second); ++ii)
      {
        _hyperparam_vec_map[hyper_param_name][ii] = theta(std::get<0>(iter->second) + ii);
      }
    }
  }
}

std::vector<Real>
AL_ADAM::Predict_ADAM(const std::vector<Real> & inputs)
{
  RealEigenMatrix test_points(1, inputs.size());
  for (unsigned int ii = 0; ii < inputs.size(); ++ii)
  {
    test_points(0, ii) = inputs[ii];
  }
  // test_points(1, 0) = 5.0;
  // test_points(1, 1) = 10000.0;

  // _param_standardizer.computeSet(_training_params);
  _param_standardizer.getStandardized(test_points);
  // std::cout << Moose::stringify(test_points) << std::endl;

  // std::cout << Moose::stringify(test_points) << std::endl;
  // _covariance_function->loadHyperParamMap(_hyperparam_map, _hyperparam_vec_map);
  // _covariance_function->computeCovarianceMatrix(_K, _training_params, _training_params, true);
  // _K_cho_decomp = _K.llt();
  // _K_results_solve = _K_cho_decomp.solve(_training_data);
  // std::cout << test_points.rows() << std::endl;
  // std::cout << _training_params.rows() << std::endl;
  // RealEigenMatrix K_test_train(test_points.rows(), _training_params.rows());
  // _covariance_function->computeCovarianceMatrix(K_test_train, test_points, _training_params, false);
  // std::cout << "K_test_train" << Moose::stringify(K_test_train) << std::endl;


  RealEigenMatrix K_train_test(_training_params.rows(), test_points.rows());
  _covariance_function->computeCovarianceMatrix(K_train_test, _training_params, test_points, false);
  // std::cout << "K_train_test" << Moose::stringify(K_train_test) << std::endl;
  RealEigenMatrix K_test(test_points.rows(), test_points.rows());
  _covariance_function->computeCovarianceMatrix(K_test, test_points, test_points, false);
  // std::cout << "K_test" << Moose::stringify(K_test) << std::endl;
  // std::cout << "K_train_test.transpose() " << Moose::stringify(K_train_test) << std::endl;

  // Compute the predicted mean value (centered)
  RealEigenMatrix pred_value = (K_train_test.transpose() * _K_results_solve);
  // std::cout << Moose::stringify(pred_value(0,0)) << std::endl;
  // De-center/scale the value and store for return
  // _data_standardizer.computeSet(_training_data);
  _data_standardizer.getDestandardized(pred_value);

  RealEigenMatrix pred_var =
      K_test - (K_train_test.transpose() * _K_cho_decomp.solve(K_train_test));

  // Vairance computed, take sqrt for standard deviation, scale up by training data std and store
  RealEigenMatrix std_dev_mat = pred_var.array().sqrt();
  // std::cout << Moose::stringify(std_dev_mat) << std::endl;
  // _data_standardizer.computeSet(_training_data);
  _data_standardizer.getDescaled(std_dev_mat);
  // std_dev = std_dev_mat(0, 0);

  std::vector<Real> result;
  result.resize(2);
  result[0] = pred_value(0, 0);
  result[1] = std_dev_mat(0, 0);

  return result;
}

bool
AL_ADAM::needSample(const std::vector<Real> & row,
                                              dof_id_type,
                                              dof_id_type,
                                              Real & val)
{
  int N = 13;
  if (_step < N)
  {
    // std::cout << "Here ***** " << val << std::endl;
    if (_step > 1)
    {
      _outputs_sto.push_back(val);
      for (unsigned int k = 0; k < _inputs_sto.size(); ++k)
        _inputs_sto[k].push_back(row[k]); // _inputs_prev[k]
    }
    // std::cout << "Inputs 1 " << Moose::stringify(_inputs_sto[0]) << std::endl;
    // std::cout << "Inputs 2 " << Moose::stringify(_inputs_sto[1]) << std::endl;
    // std::cout << "Outputs " << Moose::stringify(_outputs_sto) << std::endl;
    _decision = true;
  } else if (_step == N)
  {
    _outputs_sto.push_back(val);
    for (unsigned int k = 0; k < _inputs_sto.size(); ++k)
      _inputs_sto[k].push_back(row[k]); // _inputs_prev[k]
    // Train();
    // std::cout << Moose::stringify(_inputs_sto) << std::endl;
    // std::cout << Moose::stringify(_outputs_sto) << std::endl;
    Train_ADAM(10000);
    // std::cout << Moose::stringify(row) << std::endl;

    std::vector<Real> result = Predict_ADAM(row);
    std::cout << Moose::stringify(result) << std::endl;
    _decision = false;
    val = result[0];

    // if (std::abs(result[0]-0.0)/result[1] > 2.0)
    // {
    //   _decision = false;
    //   val = result[0];
    // } else
    //   _decision = true;
  } else
  {
    // std::vector<Real> result = Predict_ADAM(row);
    // std::cout << Moose::stringify(result) << std::endl;
    // _decision = false;
    // val = result[0];

    if (_decision == true && _flag_sample == true)
    {
      _outputs_sto.push_back(val);
      for (unsigned int k = 0; k < _inputs_sto.size(); ++k)
        _inputs_sto[k].push_back(row[k]); // _inputs_prev[k]
      Train_ADAM(1000);
    }
    // std::cout << Moose::stringify(_inputs_sto) << std::endl;
    // std::cout << Moose::stringify(_outputs_sto) << std::endl;
    std::vector<Real> result = Predict_ADAM(row);
    // std::cout << Moose::stringify(result) << std::endl;
    Real U_val;
    if (_flag_sample == false)
      U_val = result[1] / std::abs(result[0]); // std::abs(result[0]-0.0)/result[1]; //
    else
      U_val = 0.0001; // 100; //
    // std::cout << "U function " << U_val << std::endl;
    if (_flag_sample == true)
      _flag_sample = false;
    if (U_val < 0.025) // > 2.0
    {
      // std::cout << "Here" << std::endl;
      val = result[0];
      _decision = false;
    } else
    {
      _decision = true;
      _flag_sample = true;
    }
  }
  _inputs_prev = row;
  return _decision;
}
