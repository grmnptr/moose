//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "FVDirichletBCBase.h"

class Function;

/**
 * A class for no slip velocity boundary condtions
 */
class INSFVCylinderJetControlBC : public FVDirichletBCBase
{
public:
  static InputParameters validParams();
  INSFVCylinderJetControlBC(const InputParameters & params);

  ADReal boundaryValue(const FaceInfo & fi) const override;

protected:
  const RealVectorValue _origin;
  const Real & _mass_flow;
  const Real _angle;
  const std::vector<Real> _locations;
  const Real _radius;
  const unsigned int _index;
};
