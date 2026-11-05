//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "LinearFVEnergyAdvection.h"
#include "LinearFVGradientManager.h"
#include "MooseLinearVariableFV.h"
#include "NSFVUtils.h"
#include "NS.h"

registerMooseObject("NavierStokesApp", LinearFVEnergyAdvection);

InputParameters
LinearFVEnergyAdvection::validParams()
{
  InputParameters params = LinearFVFluxKernel::validParams();
  params.addClassDescription("Represents the matrix and right hand side contributions of an "
                             "advection term for the energy e.g. h=int(cp dT). A user may still "
                             "override what quantity is advected, but the default is temperature.");
  MooseEnum advected_quantity("enthalpy temperature", "enthalpy");
  params.addParam<MooseEnum>("advected_quantity", advected_quantity, "The advected quantity");
  params.addParam<Real>("cp", "Constant specific heat value");
  params.addRequiredParam<UserObjectName>(
      "rhie_chow_user_object",
      "The rhie-chow user-object which is used to determine the face velocity.");
  params += Moose::FV::advectedInterpolationParameter();
  params.addParam<InterpolationMethodName>(
      "advected_interp_method_name",
      "Name of an FVInterpolationMethod (e.g. FVAdvectedVanLeerWeightBased) to use for the "
      "advected quantity. If given, it replaces 'advected_interp_method'.");
  return params;
}

LinearFVEnergyAdvection::LinearFVEnergyAdvection(const InputParameters & params)
  : LinearFVFluxKernel(params),
    FVInterpolationMethodInterface(this),
    _advected_quantity(getParam<MooseEnum>("advected_quantity").getEnum<AdvectedQuantityEnum>()),
    _cp(isParamValid("cp") ? getParam<Real>("cp") : 1.0),
    _mass_flux_provider(getUserObject<RhieChowMassFlux>("rhie_chow_user_object")),
    _advected_interp_coeffs(std::make_pair<Real, Real>(0, 0)),
    _adv_interp_method_object(
        isParamValid("advected_interp_method_name")
            ? &getFVAdvectedInterpolationMethod(
                  getParam<InterpolationMethodName>("advected_interp_method_name"))
            : nullptr),
    _gradient_field(_adv_interp_method_object && _adv_interp_method_object->needsGradients()
                        ? &_var.requestCellGradients(_adv_interp_method_object->gradientMethodName())
                        : nullptr),
    _advected_rhs_face_value(0.0),
    _face_mass_flux(0.0)
{
  if (_adv_interp_method_object && isParamSetByUser("advected_interp_method"))
    paramError("advected_interp_method_name",
               "Only one of 'advected_interp_method' and 'advected_interp_method_name' may be "
               "set.");

  Moose::FV::setInterpolationMethod(*this, _advected_interp_method, "advected_interp_method");

  if (isParamValid("cp") && _advected_quantity == AdvectedQuantityEnum::ENTHALPY)
    paramError("cp", "cp should not be specified for enthalpy advection");

  if (!isParamValid("cp") && _advected_quantity == AdvectedQuantityEnum::TEMPERATURE)
    paramError("cp", "cp should be specified for temperature advection");
}

Real
LinearFVEnergyAdvection::computeElemMatrixContribution()
{
  return _cp * _advected_interp_coeffs.first * _face_mass_flux * _current_face_area;
}

Real
LinearFVEnergyAdvection::computeNeighborMatrixContribution()
{
  return _cp * _advected_interp_coeffs.second * _face_mass_flux * _current_face_area;
}

Real
LinearFVEnergyAdvection::computeElemRightHandSideContribution()
{
  // rhs_face_value = phi_matrix - phi_high (deferred correction), so the explicit part of the
  // face flux, -m*rhs_face_value, moves to the right hand side with a positive sign
  return _cp * _advected_rhs_face_value * _face_mass_flux * _current_face_area;
}

Real
LinearFVEnergyAdvection::computeNeighborRightHandSideContribution()
{
  return -_cp * _advected_rhs_face_value * _face_mass_flux * _current_face_area;
}

Real
LinearFVEnergyAdvection::computeBoundaryMatrixContribution(const LinearFVBoundaryCondition & bc)
{
  const auto * const adv_bc = cast_ptr<const LinearFVAdvectionDiffusionBC *>(&bc);
  mooseAssert(adv_bc, "This should be a valid BC!");

  const auto boundary_value_matrix_contrib = adv_bc->computeBoundaryValueMatrixContribution();

  // We support internal boundaries too so we have to make sure the normal points always outward
  const auto factor = (_current_face_type == FaceInfo::VarFaceNeighbors::ELEM) ? 1.0 : -1.0;

  return _cp * boundary_value_matrix_contrib * factor * _face_mass_flux * _current_face_area;
}

Real
LinearFVEnergyAdvection::computeBoundaryRHSContribution(const LinearFVBoundaryCondition & bc)
{
  const auto * const adv_bc = cast_ptr<const LinearFVAdvectionDiffusionBC *>(&bc);
  mooseAssert(adv_bc, "This should be a valid BC!");

  // We support internal boundaries too so we have to make sure the normal points always outward
  const auto factor = (_current_face_type == FaceInfo::VarFaceNeighbors::ELEM ? 1.0 : -1.0);

  const auto boundary_value_rhs_contrib = adv_bc->computeBoundaryValueRHSContribution();
  return -_cp * boundary_value_rhs_contrib * factor * _face_mass_flux * _current_face_area;
}

void
LinearFVEnergyAdvection::setupFaceData(const FaceInfo * face_info)
{
  LinearFVFluxKernel::setupFaceData(face_info);

  // Caching the mass flux on the face which will be reused in the advection term's matrix and right
  // hand side contributions
  _face_mass_flux = _mass_flux_provider.getMassFlux(*face_info);

  // Caching the interpolation coefficients so they will be reused for the matrix and right hand
  // side terms
  if (!_adv_interp_method_object)
  {
    _advected_interp_coeffs =
        interpCoeffs(_advected_interp_method, *_current_face_info, true, _face_mass_flux);
    return;
  }

  // Only internal faces need advected interpolation results; boundary contributions are handled
  // through the linear FV boundary conditions.
  _advected_rhs_face_value = 0.0;
  if (_current_face_type != FaceInfo::VarFaceNeighbors::BOTH)
    return;

  const auto state = determineState();
  const auto & elem_info = *_current_face_info->elemInfo();
  const auto & neighbor_info = *_current_face_info->neighborInfo();
  const Real elem_value = _var.getElemValue(elem_info, state);
  const Real neighbor_value = _var.getElemValue(neighbor_info, state);
  if (_gradient_field)
  {
    _elem_grad_storage = _gradient_field->gradient(elem_info);
    _neighbor_grad_storage = _gradient_field->gradient(neighbor_info);
  }

  const auto result = _adv_interp_method_object->advectedInterpolate(*_current_face_info,
                                                                     elem_value,
                                                                     neighbor_value,
                                                                     &_elem_grad_storage,
                                                                     &_neighbor_grad_storage,
                                                                     _face_mass_flux);
  _advected_interp_coeffs = result.weights_matrix;
  _advected_rhs_face_value = result.rhs_face_value;
}
