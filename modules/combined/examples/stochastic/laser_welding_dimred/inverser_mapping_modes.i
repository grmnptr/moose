thickness=0.9e-4 # m
xmin=-0.1e-3 # m
xmax=0.75e-3 # m
ymin=${fparse -thickness}

[Mesh]
  [cmg]
    type = GeneratedMeshGenerator
    dim = 2
    xmin = ${xmin}
    xmax = ${xmax}
    ymin = ${fparse ymin}
    ymax = 0
    nx = 161
    ny = 50
  []
[]

[Variables]
  [T]
  []
  [T_pod]
  []
[]

[Problem]
  solve = false
[]

[Executioner]
  type = Steady
  solve_type = NEWTON
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre boomeramg'
[]

[UserObjects]
  [im_sol]
    type = InverseMapping
    mapping = pod_mapping_sol
    variable_to_fill = "T_pod"
    variable_to_reconstruct = "T"
    parameters = '0 0 0 0 0 0 0 0 0 1'
    execute_on = TIMESTEP_END
  []
[]

[VariableMappings]
  [pod_mapping_sol]
    type = PODMapping
    filename = pod_mapping_train_mapping_sol_pod_mapping_sol.rd
    num_modes_to_compute = 10
  []
[]

[Outputs]
  exodus = true
  execute_on = 'FINAL'
[]
