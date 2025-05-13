# Process parameters
scanning_speed=1.0 # m/s
power=25 # W (this is the effective power so multiplied by eta)
R=70e-6 # m (this is the effective radius)

# Geometric parameters
thickness=70e-6 # m
ymin=-80e-6
ymax=80e-6
xmin=-80e-6 # m
xmax=240e-6 # m
surfacetemp=300 # K (temperature at the other side of the plate)
backtemp=300

# Time stepping parameters
endtime=240e-6 # s
timestep=${fparse endtime/40} # s

[Mesh]
  [cmg]
    type = GeneratedMeshGenerator
    dim = 3
    xmin = ${xmin}
    xmax = ${xmax}
    ymin = ${ymin}
    ymax = ${ymax}
    zmin = ${fparse -thickness}
    zmax = 0
    nx = 300
    ny = 120
    nz = 60
  []
[]

[Variables]
  [T]
  []
[]

[ICs]
  [T]
    type = FunctionIC
    variable = T
    function = '(${surfacetemp} - ${backtemp}) / ${thickness} * z + ${surfacetemp}'
  []
[]

[Kernels]
  [temperature_time]
    type = ADHeatConductionTimeDerivative
    variable = T
    use_displaced_mesh = true
    density_name = 'rho'
    specific_heat = 'cp'
  []
  [temperature_conduction]
    type = ADHeatConduction
    variable = T
    thermal_conductivity = 'k'
    use_displaced_mesh = true
  []
[]

[Functions]
  [fn]
    type = PiecewiseMulticonstant
    direction = 'right'
    data_file = pulsed-source.txt
  []
[]

[BCs]
  [T_cold]
    type = DirichletBC
    variable = T
    boundary = 'back'
    value = ${backtemp}
  []
  [radiation_flux]
    type = FunctionRadiativeBC
    variable = T
    boundary = 'front'
    emissivity_function = '1'
    Tinfinity = 300
    stefan_boltzmann_constant = 5.67e-8
  []
  [weld_flux]
    type = GaussianEnergyFluxBC
    variable = T
    boundary = 'front'
    P0 = ${power}
    R = ${R}
    x_beam_coord = 'fn'
    y_beam_coord = '0'
    z_beam_coord = '0'
  []
[]

[Materials]
  [steel]
    type = LaserWeld316LStainlessSteel
    temperature = T
    use_constant_density = true
  []
[]

[Executioner]
  type = Transient
  end_time = ${endtime}
  dtmin = 1e-10
  dtmax = 5e-6
  # petsc_options_iname = '-pc_type -pc_factor_shift_type'
  # petsc_options_value = 'lu       NONZERO'
  petsc_options_iname = '-pc_type -pc_hypre_type -pc_factor_shift_type'
  petsc_options_value = 'hypre    boomeramg       NONZERO'
  petsc_options = '-snes_converged_reason -ksp_converged_reason -options_left'
  solve_type = 'NEWTON'
  # line_search = 'none'
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

[Debug]
  show_var_residual_norms = true
[]

[Outputs]
  [exodus]
    type = Exodus
    output_material_properties = true
  []
[]
