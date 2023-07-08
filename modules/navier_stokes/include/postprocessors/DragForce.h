//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

// MOOSE includes
#include "SideIntegralPostprocessor.h"
#include "MathFVUtils.h"

/**
 * This postprocessor which computes the drag force on a surface
 * defined as:
 *
 * $C_D = \int\limits_S (\sigma \vec{n}) \cdot \vec{e} dS$
 *
 * where $\sigma$ is the Cauchy stress tensor
 */
class DragForce : public SideIntegralPostprocessor
{
public:
  static InputParameters validParams();

  DragForce(const InputParameters & parameters);

protected:
  Real computeFaceInfoIntegral(const FaceInfo * fi) override;

  Real computeQpIntegral() override;

  /// Velocity components
  const Moose::Functor<Real> & _vel_x;
  const Moose::Functor<Real> * _vel_y;
  const Moose::Functor<Real> * _vel_z;

  /// The dynamic viscosity
  const Moose::Functor<Real> & _mu;

  /// Pressure field
  const Moose::Functor<Real> & _pressure;

  /// The direction in which the drag/lift is measured
  const RealVectorValue _direction;
};
