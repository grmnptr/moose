//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "MappingBase.h"
#include "ParallelSolutionStorage.h"
#include <slepcsvd.h>
#include "libmesh/parallel_object.h"
#include "libmesh/petsc_vector.h"

/**
 * Class which provides a Proper Orthogonal Decomposition (POD)-based mapping between
 * full-order and reduced-order spaces.
 */
class PODMapping : public MappingBase
{
public:
  static InputParameters validParams();
  PODMapping(const InputParameters & parameters);

  ~PODMapping();

  virtual void buildMapping(const VariableName & vname) override;

  ///@{
  /**
   * Methods used for mapping full-order solutions for a given variable
   * onto a latent space
   */
  void map(const VariableName & vname,
           const DenseVector<Real> & full_order_vector,
           std::vector<Real> & reduced_order_vector) const override;

  void map(const VariableName & vname,
           const unsigned int global_sample_i,
           std::vector<Real> & reduced_order_vector) const override;
  ///@}

  /**
   * Method used for mapping reduced-order solutions for a given variable
   * onto the full-order space
   */
  void inverse_map(const VariableName & vname,
                   const std::vector<Real> & reduced_order_vector,
                   DenseVector<Real> & full_order_vector) const override;

  /**
   * Return all of the left basis functions for a given variable
   * @param vname The name of the variable.
   */
  const std::vector<DenseVector<Real>> & leftBasis(const VariableName & vname)
  {
    mooseAssert(_left_basis_functions.find(vname) != _left_basis_functions.end(),
                "We should have the requested variable!");
    return _left_basis_functions[vname];
  }

  /**
   * Return all of the right basis functions for a given variable
   * @param vname The name of the variable.
   */
  const std::vector<DenseVector<Real>> & rightBasis(const VariableName & vname)
  {
    mooseAssert(_right_basis_functions.find(vname) != _right_basis_functions.end(),
                "We should have the requested variable!");
    return _right_basis_functions[vname];
  }

  /**
   * Return all of the singular values for a given variable
   * @param vname The name of the variable.
   */
  const std::vector<Real> & singularValues(const VariableName & vname)
  {
    mooseAssert(_singular_values.find(vname) != _singular_values.end(),
                "We should have the requested variable!");
    return _singular_values[vname];
  }

  /**
   * Get the `base_i`-th left basis function for a given variable
   * @param vname The name of the variable
   * @param base_i The index of the basis function
   */
  const DenseVector<Real> & leftBase(const VariableName & vname, const unsigned int base_i);

  /**
   * Get the `base_i`-th right basis function for a given variable
   * @param vname The name of the variable
   * @param base_i The index of the basis function
   */
  const DenseVector<Real> & rightBase(const VariableName & vname, const unsigned int base_i);

protected:
  /**
   * Determine the number of basis functions needed for a given variable based on the information
   * on the eigenvalues.
   * @param vname The name of the variable
   * @param converged_evs Vector of converged eigenvalues
   */
  dof_id_type determineNumberOfModes(const VariableName & vname,
                                     const std::vector<Real> & converged_evs);

  /// The number of modes requested by the user
  const std::vector<dof_id_type> & _num_modes;

  /// The energy thresholds for truncation of the number of modes, defined by the user
  const std::vector<Real> & _energy_threshold;

  /// Restartable container holding the basis functions for each variable
  std::map<VariableName, std::vector<DenseVector<Real>>> & _left_basis_functions;

  /// Restartable container holding the basis functions for each variable
  std::map<VariableName, std::vector<DenseVector<Real>>> & _right_basis_functions;

  /// Restartable container holding the singular values
  std::map<VariableName, std::vector<Real>> & _singular_values;

  /// Variable holding additional petsc options for the singular value solve
  const std::string & _extra_slepc_options;

private:
  /// Link to the parallel storage which holds the solution fields that are used for the SVD
  ParallelSolutionStorage * _parallel_storage;

#if !PETSC_VERSION_LESS_THAN(3, 14, 0)
  /// Storage for SLEPC's SVD objects for each variable.
  std::map<VariableName, SVD> _svds;
#endif

  /// Bool to decide if we already have the SVD or not to make sure it is
  /// not computed multiple times unless the user requests it
  std::map<VariableName, bool> _computed_svd;
};
