pi = 3.14159265359
period = 0.25

[Functions]
  [src_func]
    type = ParsedFunction
    value = "200*sin(${pi}/${period}*t)"
  []
[]

[Mesh]
  [msh]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 100
    xmin = -0.5
    xmax = 0.5
    ny = 100
    ymin = -0.5
    ymax = 0.5
  []
  [source_domain]
    type = ParsedSubdomainMeshGenerator
    input = msh
    combinatorial_geometry = '(x<0.2 & x>-0.2) & (y<0.2 & y>-0.2)'
    block_id=1
  []
[]

[Variables]
  [T]
    initial_condition = 0
  []
[]

[Kernels]
  [diffusion]
    type = MatDiffusion
    variable = T
    diffusivity = diff_coeff
  []
  [source]
    type = BodyForce
    variable = T
    function = src_func
    block = 1
  []
  [anti_source]
    type = BodyForce
    variable = T
    value = 0
    block = 1
  []
  [time_deriv]
    type = TimeDerivative
    variable = T
  []
[]

[Materials]
  [diff_coeff]
    type = ParsedMaterial
    f_name = diff_coeff
    function = '3'
  []
[]

[BCs]
  [neumann_rest]
    type = NeumannBC
    variable = T
    boundary = 'left right top bottom'
    value = 0
  []
[]

[Executioner]
  type = Transient
  num_steps = 100
  dt = 0.0025
  solve_type = NEWTON
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  timestep_tolerance = 1e-6
  line_search = 'none'
[]

[Postprocessors]
  [T_max]
    type = NodalExtremeValue
    variable = T
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [T_min]
    type = NodalExtremeValue
    variable = T
    value_type = min
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [control_value]
    type = Receiver
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[Controls]
  [src_control]
    type = LibtorchNeuralNetControl
    parameters = "Kernels/anti_source/value"
    postprocessors = "control_value"
    responses = 'T_min T_max'
  []
[]

[Reporters]
  [T_reporter]
    type = AccumulateReporter
    reporters = 'T_min/value T_max/value control_value/value'
  []
[]

[Outputs]
  [csv]
    type = CSV
    execute_on = FINAL
  []
  exodus = true
[]
