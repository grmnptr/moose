//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "INSFVCylinderJetControlBC.h"
#include "Function.h"
#include "SubProblem.h"
#include "MooseMesh.h"
#include "FaceInfo.h"

#include "libmesh/vector_value.h"

registerMooseObject("NavierStokesApp", INSFVCylinderJetControlBC);

InputParameters
INSFVCylinderJetControlBC::validParams()
{
  InputParameters params = FVDirichletBCBase::validParams();
  params.addRequiredParam<RealVectorValue>("origin", "The origin of the cylinder.");
  params.addRequiredParam<Real>("mass_flow", "The massflow for the cylinder.");
  params.addRequiredParam<Real>("angle", "The angle of the actuator in degrees.");
  params.addRequiredParam<std::vector<Real>>(
      "locations", "The locations of the midpoints of the holes in degrees.");
  params.addRequiredParam<Real>("radius", "radius of the cylinder.");
  MooseEnum momentum_component("x=0 y=1 z=2");
  params.addRequiredParam<MooseEnum>(
      "momentum_component",
      momentum_component,
      "The component of the momentum equation that this bc applies to.");

  params.declareControllable("mass_flow");
  params.addClassDescription("Implements a jet control BC.");
  return params;
}

INSFVCylinderJetControlBC::INSFVCylinderJetControlBC(const InputParameters & params)
  : FVDirichletBCBase(params),
    _origin(getParam<RealVectorValue>("origin")),
    _mass_flow(getParam<Real>("mass_flow")),
    _angle(getParam<Real>("angle")),
    _locations(getParam<std::vector<Real>>("locations")),
    _radius(getParam<Real>("radius")),
    _index(getParam<MooseEnum>("momentum_component"))
{
}

ADReal
INSFVCylinderJetControlBC::boundaryValue(const FaceInfo & fi) const
{
  const Real omega = _angle / 180 * libMesh::pi;
  RealVectorValue face_center = fi.faceCentroid() - _origin;
  face_center /= face_center.norm();
  Real theta = acos(face_center(0));
  Real value = 0.0;
  if (face_center(1) > 0)
  {
    const Real theta_0 = _locations[0] * libMesh::pi / 180;
    value = _mass_flow * libMesh::pi / (2 * _radius * _radius * omega) *
            cos(libMesh::pi / omega * (theta_0 - theta));
  }
  else
  {
    const Real theta_0 = _locations[1] * libMesh::pi / 180;
    theta = 2 * libMesh::pi - theta;
    value = -_mass_flow * libMesh::pi / (2 * _radius * _radius * omega) *
            cos(libMesh::pi / omega * (theta_0 - theta));
  }

  return value * face_center(_index);
}
