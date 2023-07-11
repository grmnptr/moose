//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "DragForce.h"
#include "MathFVUtils.h"
#include "NSFVUtils.h"
#include "NS.h"
#include "MooseFunctorArguments.h"

#include <cmath>

registerMooseObject("NavierStokesApp", DragForce);

InputParameters
DragForce::validParams()
{
  InputParameters params = SideIntegralPostprocessor::validParams();
  params.addClassDescription(
      "Computes the drag force on a surface. Depending on the principal direction, one can "
      "use this object for the computation of the lift coefficient as well.");
  params.addRequiredParam<MooseFunctorName>("vel_x", "The velocity in direction x.");
  params.addParam<MooseFunctorName>("vel_y", "The velocity in direction y.");
  params.addParam<MooseFunctorName>("vel_z", "The velocity in direction z.");
  params.addRequiredParam<MooseFunctorName>(NS::mu, "The dynamic viscosity.");
  params.addRequiredParam<MooseFunctorName>(NS::pressure, "The pressure functor.");
  params.addRequiredParam<RealVectorValue>("principal_direction",
                                           "The direction in which the drag is computed. This can "
                                           "be used to compute the lift coefficient as well.");
  return params;
}

DragForce::DragForce(const InputParameters & parameters)
  : SideIntegralPostprocessor(parameters),
    _vel_x(getFunctor<Real>("vel_x")),
    _vel_y(isParamValid("vel_y") ? &getFunctor<Real>("vel_y") : nullptr),
    _vel_z(isParamValid("vel_z") ? &getFunctor<Real>("vel_z") : nullptr),
    _mu(getFunctor<Real>(NS::mu)),
    _pressure(getFunctor<Real>(NS::pressure)),
    _direction(getParam<RealVectorValue>("principal_direction"))
{
  // _qp_integration = !getFieldVar("vel_x", 0)->isFV();
  // if (_qp_integration)
  //   mooseError("Drag coefficient computation is only supported for finite volume variables!");
  _qp_integration = false;
}

Real
DragForce::computeFaceInfoIntegral(const FaceInfo * fi)
{
  mooseAssert(fi, "We should have a face info in " + name());
  const auto state = determineState();
  const auto face_arg =
      Moose::FaceArg({fi, Moose::FV::LimiterType::CentralDifference, true, false, nullptr});

  const auto elem_arg = Moose::ElemArg({fi->elemPtr(), false});

  RealTensorValue velocity_gradient;
  RealTensorValue pressure_term;
  RealTensorValue normal_term;
  Real pressure = _pressure(face_arg, state);
  Real mu = _mu(face_arg, state);
  RealVectorValue cell_velocity = 0;
  RealVectorValue face_velocity = 0;
  if (_mesh.dimension() == 1)
  {
    const auto & grad_u = _vel_x.gradient(face_arg, state);
    velocity_gradient = RealTensorValue(grad_u);
    pressure_term(0, 0) = -pressure;
    normal_term(0, 0) = 1;
    cell_velocity(0) = _vel_x(elem_arg, state);
    face_velocity(0) = _vel_x(face_arg, state);
  }
  else if (_mesh.dimension() == 2)
  {
    // const auto & grad_v = _v_var->adGradSln(face_info, skewness_correction);
    const auto & grad_u = _vel_x.gradient(face_arg, state);
    const auto & grad_v = _vel_y->gradient(face_arg, state);
    velocity_gradient = RealTensorValue(grad_u, grad_v);
    pressure_term(0, 0) = -pressure;
    pressure_term(1, 1) = -pressure;
    normal_term(0, 0) = 1;
    normal_term(1, 1) = 1;
    cell_velocity(0) = _vel_x(elem_arg, state);
    cell_velocity(1) = (*_vel_y)(elem_arg, state);
    face_velocity(0) = _vel_x(face_arg, state);
    face_velocity(1) = (*_vel_y)(face_arg, state);
  }
  else // if (_dim == 3)
  {
    const auto & grad_u = _vel_x.gradient(face_arg, state);
    const auto & grad_v = _vel_y->gradient(face_arg, state);
    const auto & grad_w = _vel_z->gradient(face_arg, state);
    velocity_gradient = RealTensorValue(grad_u, grad_v, grad_w);
    pressure_term(0, 0) = -pressure;
    pressure_term(1, 1) = -pressure;
    pressure_term(2, 2) = -pressure;
    normal_term(0, 0) = 1;
    normal_term(1, 1) = 1;
    normal_term(2, 2) = 1;
    cell_velocity(0) = _vel_x(elem_arg, state);
    cell_velocity(1) = (*_vel_y)(elem_arg, state);
    cell_velocity(2) = (*_vel_z)(elem_arg, state);
    face_velocity(0) = _vel_x(face_arg, state);
    face_velocity(1) = (*_vel_y)(face_arg, state);
    face_velocity(2) = (*_vel_z)(face_arg, state);
  }

  auto sheer_force = -mu *
                     (cell_velocity - face_velocity -
                      (cell_velocity - face_velocity) * fi->normal() * fi->normal()) /
                     std::abs(fi->dCN() * fi->normal());

  // auto normal_tensor = normal_term - libMesh::outer_product(fi->normal(), fi->normal());
  // auto normal_gradient = ((velocity_gradient + velcoity_gradient.transpose()) * normal_tensor);

  // Real value = 0;
  // if (std::abs(_direction(0)) == 1)
  //   value = (mu * (velocity_gradient) + pressure_term) * fi->normal() * _direction;
  // else
  //   value = (mu * velocity_gradient.transpose() + pressure_term) * fi->normal() * _direction;

  // return (mu * (velocity_gradient + velocity_gradient.transpose()) + pressure_term) *
  // fi->normal() *
  //        _direction;

  return (sheer_force * _direction + pressure_term * fi->normal() * _direction);
}

Real
DragForce::computeQpIntegral()
{
  return 0.0;
}
