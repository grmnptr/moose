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
[]

# [Problem]
#   restart_file_base = executioner_postprocessor_out_cp/LATEST
# []

[Preconditioning]
  active = SMP
  [FSP]
    type = FSP
    # It is the starting point of splitting
    topsplit = 'up' # 'up' should match the following block name
    [up]
      splitting = 'u p' # 'u' and 'p' are the names of subsolvers
      splitting_type  = schur
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
  end_time = 60
  dtmax = 1e-2
  timestep_tolerance = 1e-10
  scheme = 'bdf2'
[]

[Outputs]
  csv = true
  checkpoint = true
  [exo]
    type = Exodus
    interval = 10
    overwrite = false
  []
[]

[Postprocessors]
  [Re]
    type = ParsedPostprocessor
    function = 'rho * U * D / mu'
    constant_names = 'rho U D mu'
    constant_expressions = '${rho} ${fparse 2/3*inlet_velocity} ${fparse 2*circle_radius} ${mu}'
    pp_names = ''
  []
  [point_vel_x]
    type = PointValue
    point = '${fparse (x_max-x_min)/2} ${fparse (y_max-y_min)/2} 0'
    variable = 'vel_x'
  []
  [point_vel_y]
    type = PointValue
    point = '${fparse (x_max-x_min)/2} ${fparse (y_max-y_min)/2} 0'
    variable = 'vel_y'
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
