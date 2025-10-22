//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ADTorchScriptMaterial.h"
#include "HeatConductionModel.h"
#include "LibtorchUtils.h"

registerMooseObject("ThermalHydraulicsApp", ADTorchScriptMaterial);

InputParameters
ADTorchScriptMaterial::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription("Computes solid thermal properties as a function of temperature");
  // Coupled variables
  params.addRequiredCoupledVar("u", "Coupled variable");
  params.addRequiredParam<UserObjectName>("torch_uo", "TC UO.");
  return params;
}

ADTorchScriptMaterial::ADTorchScriptMaterial(const InputParameters & parameters)
  : Material(parameters),
    _thermal_conductivity(declareADProperty<Real>("thermal_conductivity")),
    _temp(adCoupledValue("u")),
    _torchscript_uo(getUserObject<TorchScriptUserObject>("torch_uo"))
{
}

void
ADTorchScriptMaterial::computeQpProperties()
{
  std::vector<Real> input_vector(1, raw_value(_temp[_qp]));

  torch::Tensor input_tensor;
  LibtorchUtils::vectorToTensor(input_vector, input_tensor);
  input_tensor = input_tensor.transpose(0, 1);

  _thermal_conductivity[_qp] = _torchscript_uo.evaluate(input_tensor).item<Real>();
}
