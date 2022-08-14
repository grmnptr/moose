refinement=6
mu=3.2215e-3
rho=1
inlet_velocity=1
velocity_interp_method='rc'
advected_interp_method='upwind'

[GlobalParams]
  rhie_chow_user_object = 'rc'
[]

[Mesh]
  [ccmg]
    type = ConcentricCircleMeshGenerator
    num_sectors = ${fparse refinement*2}
    radii = '0.2546 0.3368'
    rings = '4 ${fparse 2*refinement} ${refinement}'
    has_outer_square = on
    pitch = 1
    preserve_volumes = off
    smoothing_max_it = 3
  []
  [in_between]
    type = SideSetsBetweenSubdomainsGenerator
    input = ccmg
    primary_block = 2
    paired_block = 1
    new_boundary = 'circle'
  []
  [delete]
    type = BlockDeletionGenerator
    input = in_between
    block = '1'
  []
  [final_ccmg]
    type = RenameBlockGenerator
    input = delete
    old_block = '2 3'
    new_block = '0 0'
  []
  [left]
    type = GeneratedMeshGenerator
    dim = 2
    xmin = -1.5
    xmax = -0.5
    ymin = -0.5
    ymax = 0.5
    nx = ${fparse refinement*5}
    ny = ${fparse refinement*4+2}
  []
  [right]
    type = GeneratedMeshGenerator
    dim = 2
    xmin = 0.5
    xmax = 4.5
    ymin = -0.5
    ymax = 0.5
    nx = ${fparse refinement*20}
    ny = ${fparse refinement*4+2}
  []
  [stitch_left]
    type = StitchedMeshGenerator
    inputs = 'final_ccmg left'
    stitch_boundaries_pairs = 'left right'
  []
  [combined_middle]
    type = StitchedMeshGenerator
    inputs = 'stitch_left right'
    stitch_boundaries_pairs = 'right left'
  []

  [middle_top_sideset]
    input = combined_middle
    type = ParsedGenerateSideset
    combinatorial_geometry = 'y > 0.49999999'
    normal = '0 1 0'
    new_sideset_name = 'middle_top'
  []
  [middle_bottom_sideset]
    input = middle_top_sideset
    type = ParsedGenerateSideset
    combinatorial_geometry = 'y < -0.4999999'
    normal = '0 -1 0'
    new_sideset_name = 'bottom_boundary'
  []

  [top_left_block]
    type = GeneratedMeshGenerator
    dim = 2
    xmin = -1.5
    xmax = -0.5
    ymin = 0.5
    ymax = ${fparse 0.5 + 2. / 16.}
    nx = ${fparse refinement*5}
    ny = ${fparse ceil(refinement/2)}
  []
  [top_middle_block]
    type = GeneratedMeshGenerator
    dim = 2
    xmin = -0.5
    xmax = 0.5
    ymin = 0.5
    ymax = ${fparse 0.5 + 2. / 16.}
    nx = ${fparse refinement*4+2}
    ny = ${fparse ceil(refinement/2)}
  []
  [top_right_block]
    type = GeneratedMeshGenerator
    dim = 2
    xmin = 0.5
    xmax = 4.5
    ymin = 0.5
    ymax = ${fparse 0.5 + 2. / 16.}
    nx = ${fparse refinement*20}
    ny = ${fparse ceil(refinement/2)}
  []
  [stitch_top_left]
    type = StitchedMeshGenerator
    inputs = 'top_middle_block top_left_block'
    stitch_boundaries_pairs = 'left right'
  []
  [combined_top]
    type = StitchedMeshGenerator
    inputs = 'stitch_top_left top_right_block'
    stitch_boundaries_pairs = 'right left'
  []
  [top_bottom_sideset]
    input = combined_top
    type = ParsedGenerateSideset
    combinatorial_geometry = 'y < 0.5000001'
    normal = '0 -1 0'
    new_sideset_name = 'top_bottom'
  []
  [combined_middle_top]
    type = StitchedMeshGenerator
    inputs = 'top_bottom_sideset middle_bottom_sideset'
    stitch_boundaries_pairs = 'top_bottom middle_top'
  []
  [create_fused_top_sideset]
    input = combined_middle_top
    type = ParsedGenerateSideset
    combinatorial_geometry = 'y > ${fparse 0.999999 * 0.5 + 2. / 16.}'
    normal = '0 1 0'
    new_sideset_name = 'top_boundary'
  []
  [create_fused_left_sideset]
    input = create_fused_top_sideset
    type = ParsedGenerateSideset
    combinatorial_geometry = 'x < -1.499999'
    normal = '-1 0 0'
    new_sideset_name = 'left_boundary'
  []
  [create_fused_right_sideset]
    input = create_fused_left_sideset
    type = ParsedGenerateSideset
    combinatorial_geometry = 'x > 4.4999999'
    normal = '1 0 0'
    new_sideset_name = 'right_boundary'
  []
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
    boundary = 'left_boundary'
    function = ${inlet_velocity}
  []
  [inlet_y]
    type = INSFVInletVelocityBC
    variable = vel_y
    boundary = 'left_boundary'
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
  [walls_x]
    type = INSFVNaturalFreeSlipBC
    variable = vel_x
    boundary = 'top_boundary bottom_boundary'
    momentum_component = 'x'
  []
  [walls_y]
    type = INSFVNaturalFreeSlipBC
    variable = vel_y
    boundary = 'top_boundary bottom_boundary'
    momentum_component = 'y'
  []

  [outlet_p]
    type = INSFVOutletPressureBC
    variable = pressure
    boundary = 'right_boundary'
    function = 0
  []
[]

[Executioner]
  type = Transient
  solve_type = 'NEWTON'
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  petsc_options = '-snes_converged_reason -ksp_converged_reason'
  line_search = 'none'
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10
  nl_max_its = 10
  end_time = 15
  dtmax = 1e-1
  dtmin = 1e-5
  scheme = 'bdf2'
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 1e-4
    optimal_iterations = 6
    growth_factor = 1.5
  []
[]

[Outputs]
  exodus = true
  csv = true
  checkpoint = true
[]

[Postprocessors]
  [Re]
    type = ParsedPostprocessor
    function = 'rho * U * D / mu'
    constant_names = 'rho U D mu'
    constant_expressions = '${rho} ${inlet_velocity} ${fparse 2 * .2546} ${mu}'
    pp_names = ''
  []
  [Strouhal_Number]
    type = ParsedPostprocessor
    function = '0.198*(1-(19.7/(rho*U*D/mu)))'
    constant_names = 'rho U D mu'
    constant_expressions = '${rho} ${inlet_velocity} ${fparse 2 * .2546} ${mu}'
    pp_names = ''
  []
  [St_PVS]
    type = ParsedPostprocessor
    function = '.285 + (-1.3897/sqrt(rho*U*D/mu)) + 1.8061/(rho*U*D/mu)'
    constant_names = 'rho U D mu'
    constant_expressions = '${rho} ${inlet_velocity} ${fparse 2 * .2546} ${mu}'
    pp_names = ''
  []
  [Frequency]
    type = ParsedPostprocessor
    function = '((0.198*U)/D)*(1-(19.7/(rho*U*D/mu)))'
    constant_names = 'rho U D mu'
    constant_expressions = '${rho} ${inlet_velocity} ${fparse 2 * .2546} ${mu}'
    pp_names = ''
  []
  [Frequency_PVS]
    type = ParsedPostprocessor
    function = 'St_PVS*U/D'
    constant_names = 'U D'
    constant_expressions = '${inlet_velocity} ${fparse 2 * .2546}'
    pp_names = 'St_PVS'
  []
  [element_44146_x]
    type = ElementalVariableValue
    variable = 'vel_x'
    elementid = 21082
  []
  [element_44146_y]
    type = ElementalVariableValue
    variable = 'vel_y'
    elementid = 21082
  []
[]
