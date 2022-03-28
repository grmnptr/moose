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
#include "Standardizer.h"
#include <Eigen/Dense>

#include "CovarianceFunctionBase.h"
#include "CovarianceInterface.h"

class ActiveLearningGP : public ActiveLearningReporterTempl<Real>, public CovarianceInterface
{
public:
  static InputParameters validParams();
  ActiveLearningGP(const InputParameters & parameters);

  CovarianceFunctionBase * getCovarPtr() const { return _covariance_function; }

  ///@{
  /**
   * Declare model data for loading from file as well as restart
   */
  // MOOSEDOCS_BEGIN
  template <typename T>
  T & declareModelData(const std::string & data_name);

  template <typename T>
  T & declareModelData(const std::string & data_name, const T & value);
  // MOOSEDOCS_END
  ///@}

protected:

  void SetupData(const std::vector<std::vector<Real>> & inputs, const std::vector<Real> & outputs);

  void Train();

  // write stored hyperparam_vecs to PetscVec
  void mapToVec(libMesh::PetscVector<Number> & theta) const;

  // loads PetscVec to stored hyperparam_vecs
  void vecToMap(libMesh::PetscVector<Number> & theta);

  // Sets bounds for hyperparameters
  void buildHyperParamBounds(libMesh::PetscVector<Number> & theta_l,
                             libMesh::PetscVector<Number> & theta_u) const;

  // Wrapper for PETSc function callback
  static PetscErrorCode
  FormFunctionGradientWrapper(Tao tao, Vec theta, PetscReal * f, Vec Grad, void * ptr);

  // Computes Gradient of the loss function
  void FormFunctionGradient(Tao tao, Vec theta, PetscReal * f, Vec Grad);

  PetscErrorCode FormInitialGuess(ActiveLearningGP * GP_ptr, Vec theta);

  // Routine to perform hyperparameter tuning
  PetscErrorCode hyperparamTuning();

  std::vector<Real> Predict(const std::vector<Real> & inputs);

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
  // const std::vector<Real> & _output_value;

private:
  /// Name for the meta data associated with training
  const std::string _model_meta_data_name;

  /**
   * Internal function used by public declareModelData methods.
   */
  template <typename T>
  RestartableData<T> & declareModelDataHelper(const std::string & data_name);

  /// Track the current step of the main App
  const int & _step;

  /// The adaptive Monte Carlo sampler
  Sampler & _sampler;

  /// Standardizer for use with params (x)
  StochasticTools::Standardizer & _param_standardizer;

  /// Standardizer for use with data (y)
  StochasticTools::Standardizer & _data_standardizer;

  CovarianceFunctionBase * _covariance_function;

  bool & _flag_sample;

  bool _do_tuning;

  std::string _tao_options;

  bool _show_tao;

  Parallel::Communicator _tao_comm;

  std::string & _covar_type;

  std::unordered_map<std::string, std::tuple<unsigned int, unsigned int, Real, Real>> _tuning_data;

  unsigned int _num_tunable;

  /// Paramaters (x) used for training, along with statistics
  RealEigenMatrix & _training_params;

  /// Data (y) used for training
  RealEigenMatrix _training_data;

  /// An _n_sample by _n_sample covariance matrix constructed from the selected kernel function
  RealEigenMatrix & _K;

  /// A solve of Ax=b via Cholesky.
  RealEigenMatrix & _K_results_solve;

  /// Cholesky decomposition Eigen object
  Eigen::LLT<RealEigenMatrix> _K_cho_decomp;

  /// Switch for training param (x) standardization
  bool _standardize_params;

  /// Switch for training data(y) standardization
  bool _standardize_data;

  /// Scalar hyperparameters. Stored for use in surrogate
  std::unordered_map<std::string, Real> & _hyperparam_map;

  /// Vector hyperparameters. Stored for use in surrogate
  std::unordered_map<std::string, std::vector<Real>> & _hyperparam_vec_map;

  std::vector<std::vector<Real>> _inputs_sto;

  std::vector<Real> _outputs_sto;

  /// store the default value
  // Real _default_value;

  std::vector<Real> _inputs_prev;

  bool _decision;

  std::vector<Real> _len_sto;

  Real _amp_sto;

  Real _threshold;
};

template <typename T>
T &
ActiveLearningGP::declareModelData(const std::string & data_name)
{
  RestartableData<T> & data_ref = declareModelDataHelper<T>(data_name);
  return data_ref.set();
}

template <typename T>
T &
ActiveLearningGP::declareModelData(const std::string & data_name, const T & value)
{
  RestartableData<T> & data_ref = declareModelDataHelper<T>(data_name);
  data_ref.set() = value;
  return data_ref.set();
}

template <typename T>
RestartableData<T> &
ActiveLearningGP::declareModelDataHelper(const std::string & data_name)
{
  auto data_ptr = std::make_unique<RestartableData<T>>(data_name, nullptr);
  RestartableDataValue & value =
      _app.registerRestartableData(data_name, std::move(data_ptr), 0, false, _model_meta_data_name);
  RestartableData<T> & data_ref = static_cast<RestartableData<T> &>(value);
  return data_ref;
}
