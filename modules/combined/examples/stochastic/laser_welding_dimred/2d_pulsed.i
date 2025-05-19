!include parameters.i

!include mesh.i

[Variables]
  [vel]
    family = LAGRANGE_VEC
  []
  [T]
  []
  [p]
  []
  [disp_x]
  []
  [disp_y]
  []
[]

!include physics_objects_pulsed.i

[Executioner]
  type = Transient
  end_time = ${endtime}
  dtmin = 1e-10
  dtmax = 5e-6
  petsc_options_iname = '-pc_type -pc_factor_shift_type'
  petsc_options_value = 'lu       NONZERO'
  petsc_options = '-snes_converged_reason -ksp_converged_reason -options_left'
  solve_type = 'NEWTON'
  line_search = 'none'
  nl_max_its = 5
  l_max_its = 100
  [TimeStepper]
    type = IterationAdaptiveDT
    optimal_iterations = 5
    iteration_window = 1
    dt = ${timestep}
    linear_iteration_ratio = 1e6
    growth_factor = 1.25
  []
[]

[Outputs]
  exodus = true
  csv = true
[]

[Postprocessors]
  [const]
    type = ConstantPostprocessor
    value = 0
  []
[]

[Debug]
  show_var_residual_norms = true
[]
