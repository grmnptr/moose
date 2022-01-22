[Mesh]
  type = GeneratedMesh
  dim = 2
  xmin = 0
  xmax = 3
  ymin = 1
  ymax = 4
  nx = 100
  ny = 100
[]

# [Problem]
#   coord_type = RZ
# []

[Variables]
  [./p]
  [../]
[]

# [AuxVariables]
#   [./vel_p]
#   [../]
#   [./accel_p]
#   [../]
#   [./accel_x]
#   [../]
#   [./accel_y]
#   [../]
# []

[Kernels]
  # [./inertia_p]
  #   type = InertialForce
  #   variable = 'p'
  # [../]
  [./inertia_p]
    type = AcousticInertia
    variable = p
  [../]
  [./diff_p]
    type = ADMatDiffusion
    variable = 'p'
    diffusivity = 'Diff'
  [../]
[]

[DiracKernels]
  [point_source]
    type = FunctionDiracSource
    variable = p
    function = switch_off
    point = '0.0 1.75 0.0'
  []
[]

[Functions]
  [./switch_off]
    type = ParsedFunction
    value = 'if(t < 0.1, 5e-4, 0)'
  [../]
[]

[Materials]
  # [./density]
  #   type = GenericConstantMaterial
  #   prop_names = 'density'
  #   prop_values = '444.44'
  # [../]
  [./co_sq]
    type = GenericConstantMaterial
    prop_names = inv_co_sq
    prop_values = 4.44e-7
  [../]
  [./diff]
    type = ADGenericConstantMaterial
    prop_names = 'Diff'
    prop_values = '1000'
  [../]
[]

[BCs]
  [axial]
    type = DiffusionFluxBC
    variable = p
    boundary = left
  []
[]


[Postprocessors]
  [./p_1]
    type = PointValue
    variable = p
    # point = '0.5 1.01 0.0'
    point = '0 1.75 0'
  [../]
[]

[Executioner]
  type = Transient
  solve_type = 'NEWTON'
  petsc_options_iname = '-pc_type -ksp_grmres_restart -sub_ksp_type -sub_pc_type -pc_asm_overlap'
  petsc_options_value = 'asm      31                  preonly       lu           1'
  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-8
  # automatic_scaling = true
  end_time = 0.995 # 2.1
  dt = 1.5e-3
  # end_time = 1
  # dt = 1e-3
  [TimeIntegrator]
    type = NewmarkBeta
  []
[]

[Outputs]
  [./csv]
    type = CSV
    delimiter = ','
    file_base = 'no_rz_press'
  [../]
  [./exodus]
    type = Exodus
    interval = 100
    file_base = no_rz
  [../]
[]
