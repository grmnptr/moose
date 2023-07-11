velocity_interp_method = 'rc'

[GlobalParams]
  rhie_chow_user_object = 'rc'
[]

[UserObjects]
  [rc]
    type = INSFVRhieChowInterpolator
    u = vel_x
    v = vel_y
    pressure = pressure
  []
[]

[Variables]
  [vel_x]
    type = INSFVVelocityVariable
    two_term_boundary_expansion = true
  []
  [vel_y]
    type = INSFVVelocityVariable
    two_term_boundary_expansion = true
  []
  [pressure]
    type = INSFVPressureVariable
    two_term_boundary_expansion = true
  []
[]

[FVKernels]
  [mass]
    type = INSFVMassAdvection
    variable = pressure
    advected_interp_method = ${advected_interp_method}
    velocity_interp_method = ${velocity_interp_method}
    rho = ${rho}
  []
  [u_time]
    type = INSFVMomentumTimeDerivative
    variable = vel_x
    rho = ${rho}
    momentum_component = 'x'
  []
  [u_advection]
    type = INSFVMomentumAdvection
    variable = vel_x
    advected_interp_method = ${advected_interp_method}
    velocity_interp_method = ${velocity_interp_method}
    rho = ${rho}
    momentum_component = 'x'
  []
  [u_viscosity]
    type = INSFVMomentumDiffusion
    variable = vel_x
    mu = ${mu}
    momentum_component = 'x'
  []
  [u_pressure]
    type = INSFVMomentumPressure
    variable = vel_x
    pressure = pressure
    momentum_component = 'x'
  []
  [v_time]
    type = INSFVMomentumTimeDerivative
    variable = vel_y
    rho = ${rho}
    momentum_component = 'y'
  []
  [v_advection]
    type = INSFVMomentumAdvection
    variable = vel_y
    advected_interp_method = ${advected_interp_method}
    velocity_interp_method = ${velocity_interp_method}
    rho = ${rho}
    momentum_component = 'y'
  []
  [v_viscosity]
    type = INSFVMomentumDiffusion
    variable = vel_y
    mu = ${mu}
    momentum_component = 'y'
  []
  [v_pressure]
    type = INSFVMomentumPressure
    variable = vel_y
    pressure = pressure
    momentum_component = 'y'
  []
[]

[FVBCs]
  [inlet_x]
    type = INSFVInletVelocityBC
    variable = vel_x
    boundary = 'inlet'
    function = 'inlet_function'
  []
  [inlet_y]
    type = INSFVInletVelocityBC
    variable = vel_y
    boundary = 'inlet'
    function = 0
  []
  [circle_x]
    type = INSFVNoSlipWallBC
    variable = vel_x
    boundary = 'circle'
    function = 0
  []
  [circle_y]
    type = INSFVNoSlipWallBC
    variable = vel_y
    boundary = 'circle'
    function = 0
  []
  [jets_x]
    type = INSFVCylinderJetControlBC
    variable = vel_x
    boundary = 'jets'
    mass_flow = 0.0
    angle = 10
    locations = '90 270'
    origin = '0 0 0'
    radius = 0.5
    momentum_component = x
  []
  [jets_y]
    type = INSFVCylinderJetControlBC
    variable = vel_y
    boundary = 'jets'
    mass_flow = 0.0
    angle = 10
    locations = '90 270'
    origin = '0 0 0'
    radius = 0.5
    momentum_component = y
  []
  [walls_x]
    type = INSFVNoSlipWallBC
    variable = vel_x
    boundary = 'wall'
    function = 0
  []
  [walls_y]
    type = INSFVNoSlipWallBC
    variable = vel_y
    boundary = 'wall'
    function = 0
  []

  [outlet_p]
    type = INSFVOutletPressureBC
    variable = pressure
    boundary = 'outlet'
    function = 0
  []
[]
