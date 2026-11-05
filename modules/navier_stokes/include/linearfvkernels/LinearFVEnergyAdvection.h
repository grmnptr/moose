//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "LinearFVFluxKernel.h"
#include "RhieChowMassFlux.h"
#include "LinearFVAdvectionDiffusionBC.h"
#include "FVAdvectedInterpolationMethod.h"
#include "FVInterpolationMethodInterface.h"

class LinearFVGradientReader;

/**
 * An advection kernel that implements the advection term for the enthalpy in the
 * energy equation.
 */
class LinearFVEnergyAdvection : public LinearFVFluxKernel, public FVInterpolationMethodInterface
{
public:
  static InputParameters validParams();
  LinearFVEnergyAdvection(const InputParameters & params);

  virtual Real computeElemMatrixContribution() override;

  virtual Real computeNeighborMatrixContribution() override;

  virtual Real computeElemRightHandSideContribution() override;

  virtual Real computeNeighborRightHandSideContribution() override;

  virtual Real computeBoundaryMatrixContribution(const LinearFVBoundaryCondition & bc) override;

  virtual Real computeBoundaryRHSContribution(const LinearFVBoundaryCondition & bc) override;

  virtual void setupFaceData(const FaceInfo * face_info) override;

protected:
  /// The advected quantity options
  enum class AdvectedQuantityEnum
  {
    ENTHALPY = 0,
    TEMPERATURE = 1
  };

  /// The advected quantity
  AdvectedQuantityEnum _advected_quantity;

  /// The constant specific heat value. Only needed for temperature advection, with
  /// the assumption of constant cp.
  const Real _cp;

  /// The Rhie-Chow user object that provides us with the face velocity
  const RhieChowMassFlux & _mass_flux_provider;

private:
  /// Container for the current advected interpolation coefficients on the face to make sure
  /// we don't compute it multiple times for different terms.
  std::pair<Real, Real> _advected_interp_coeffs;

  /// Optional FVInterpolationMethod object for the advected quantity (e.g. a limited
  /// second-order scheme). When null, the legacy advected_interp_method enum is used.
  const FVAdvectedInterpolationMethod * const _adv_interp_method_object;

  /// Gradient field used by advected interpolation methods that require gradients
  const LinearFVGradientReader * const _gradient_field;

  /// Right hand side part of the advected face value for the current face (deferred
  /// correction methods); zero for weight-based methods and the legacy enum
  Real _advected_rhs_face_value;

  /// Reusable gradient storage used when the advected interpolation requires gradients
  VectorValue<Real> _elem_grad_storage;
  VectorValue<Real> _neighbor_grad_storage;

  /// Container for the mass flux on the face which will be reused in the advection term's
  /// matrix and right hand side contribution
  Real _face_mass_flux;

  /// The interpolation method to use for the advected quantity
  Moose::FV::InterpMethod _advected_interp_method;
};
