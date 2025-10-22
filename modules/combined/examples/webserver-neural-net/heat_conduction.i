[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 1
    nx = 20
    xmin = 0
    xmax = 1
  []
[]

[Variables]
  [u]
    initial_condition = 300
  []
[]

[UserObjects]
  [my_net]
    type = TorchScriptUserObject
    filename = "my_tc_net.pt"
    load_during_construction = true
    execute_on = INITIAL
  []
[]

[BCs]
  [left]
    type = DirichletBC
    variable = u
    boundary = 'left'
    value = 300
  []
  [right]
    type = DirichletBC
    variable = u
    boundary = 'right'
    value = 400
  []
[]

[Kernels]
  [diff]
    type = ADMatDiffusion
    variable = u
    diffusivity = thermal_conductivity
  []
  [source]
    type = BodyForce
    variable = u
    function = src
  []
[]

[Materials]
  [tc]
    type = ADTorchScriptMaterial
    torch_uo = my_net
    u = u
  []
[]

[Functions]
  [src]
    type = ParsedFunction
    expression = '5000*exp(x)'
  []
[]

[VectorPostprocessors]
  [temps]
    type = LineValueSampler
    warn_discontinuous_face_values = false
    start_point = '0 0 0'
    end_point = '1 0 0'
    num_points = 20
    sort_by = x
    variable = 'u'
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]


[Controls/web_server]
  type = WebServerControl
  execute_on = 'INITIAL FINAL'
[]

[Executioner]
  type = Steady
  solve_type = 'PJFNK'
[]

[Outputs]
  execute_on = 'timestep_end'
  csv = true
[]
