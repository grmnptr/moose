//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Material.h"
#include "SolidMaterialProperties.h"
#include "TorchScriptUserObject.h"

/**
 * A class to define materials for the solid structures in the THM application.
 */
class ADTorchScriptMaterial : public Material
{
public:
  ADTorchScriptMaterial(const InputParameters & parameters);

protected:
  virtual void computeQpProperties();

  /// The solid material properties
  ADMaterialProperty<Real> & _thermal_conductivity;

  /// Temperature in the solid structure
  const ADVariableValue & _temp;

  /// The user object that holds the neural network
  const TorchScriptUserObject & _torchscript_uo;

public:
  static InputParameters validParams();
};
