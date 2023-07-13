[Functions]
  [inlet_function]
    type = ParsedFunction
    expression = '4*U*(y-ymin)*(ymax-y)/(ymax-ymin)/(ymax-ymin)'
    symbol_names = 'U ymax ymin'
    symbol_values = '${inlet_velocity} ${y_max} ${y_min}'
  []
  #   [inlet_function]
  #     type = ParsedFunction
  #     expression = '6*((ymax-ymin)/2-y)*((ymax-ymin)/2+y)/((ymax-ymin)*(ymax-ymin))'
  #     symbol_names = 'ymax ymin'
  #     symbol_values = '${y_max} ${y_min}'
  #   []
  [reward_function]
    type = LiftDragDRLRewardFunction
    lift_pp = lift_coeff
    drag_pp = drag_coeff
    smoothing_interval = 336
  []
[]

[Problem]
  restart_file_base = executioner_postprocessor_out_cp_restart/LATEST
[]

[Preconditioning]
  active = SMP
  [FSP]
    type = FSP
    # It is the starting point of splitting
    topsplit = 'up' # 'up' should match the following block name
    [up]
      splitting = 'u p' # 'u' and 'p' are the names of subsolvers
      splitting_type = schur
      petsc_options_iname = '-pc_fieldsplit_schur_fact_type  -pc_fieldsplit_schur_precondition -ksp_gmres_restart -ksp_rtol -ksp_atol -ksp_type'
      petsc_options_value = 'full                            selfp                             300                1e-4       8e-8 fgmres'
    []
    [u]
      vars = 'vel_x vel_y'
      petsc_options_iname = '-pc_type -pc_hypre_type -ksp_type -ksp_rtol -ksp_gmres_restart -ksp_pc_side'
      petsc_options_value = 'hypre    boomeramg      gmres    5e-1      300                 right'
    []
    [p]
      vars = 'pressure'
      petsc_options_iname = '-ksp_type -ksp_gmres_restart -ksp_rtol -pc_type -ksp_pc_side'
      petsc_options_value = 'gmres    300                5e-1      jacobi    right'
    []
  []
  [SMP]
    type = SMP
    full = true
    petsc_options_iname = '-pc_type -pc_factor_mat_solver_package -pc_factor_shift_type -pc_factor_mat_ordering_type'
    petsc_options_value = 'lu  superlu_dist NONZERO nd'
  []
[]

[Executioner]
  type = Transient
  solve_type = 'NEWTON'
  # petsc_options = '-snes_converged_reason -ksp_converged_reason'
  # petsc_options_iname = '-pc_type'
  # petsc_options_value = 'lu'
  line_search = 'none'
  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-6
  nl_max_its = 10
  end_time = 20
  start_time = 0
  dtmax = 1e-2
  timestep_tolerance = 1e-10
  scheme = 'bdf2'
[]

[Outputs]
  csv = true
  [exo]
    type = Exodus
    interval = 10
    overwrite = false
  []
[]

[Postprocessors]
  [reward]
    type = FunctionValuePostprocessor
    function = reward_function
    execute_on = 'INITIAL TIMESTEP_END'
    indirect_dependencies = 'lift_coeff drag_coeff'
  []
  [p_probe_0]
    type = PointValue
    point = '0.75 1 0'
    variable = 'pressure'
  []
  [p_probe_1]
    type = PointValue
    point = '0.75 -1 0'
    variable = 'pressure'
  []
  [p_probe_2]
    type = PointValue
    point = '0.75 0 0'
    variable = 'pressure'
  []
  [p_probe_3]
    type = PointValue
    point = '0 0.7 0'
    variable = 'pressure'
  []
  [p_probe_4]
    type = PointValue
    point = '0 -0.7 0'
    variable = 'pressure'
  []

  [drag_force]
    type = DragForce
    vel_x = vel_x
    vel_y = vel_y
    mu = ${mu}
    pressure = pressure
    principal_direction = '-1 0 0'
    boundary = 'circle jets'
    outputs = none
  []
  [drag_coeff]
    type = ParsedPostprocessor
    function = '2*drag_force/rho/(avgvel*avgvel)/D'
    constant_names = 'rho avgvel D'
    constant_expressions = '${rho} ${fparse 2/3*inlet_velocity} ${fparse 2*circle_radius}'
    pp_names = 'drag_force'
  []
  [lift_force]
    type = DragForce
    vel_x = vel_x
    vel_y = vel_y
    mu = ${mu}
    pressure = pressure
    principal_direction = '0 -1 0'
    boundary = 'circle jets'
    outputs = none
  []
  [lift_coeff]
    type = ParsedPostprocessor
    function = '2*lift_force/rho/(avgvel*avgvel)/D'
    constant_names = 'rho avgvel D'
    constant_expressions = '${rho} ${fparse 2/3*inlet_velocity} ${fparse 2*circle_radius}'
    pp_names = 'lift_force'
  []
[]
