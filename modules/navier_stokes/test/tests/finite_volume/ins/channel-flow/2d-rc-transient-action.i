# Fluid properties
mu = 1.1
rho = 1.1
cp = 1.1
k = 1e-3

# Operating conditions
u_inlet = 1
T_inlet = 200
T_solid = 190
p_outlet = 10
h_fs = 0.01

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    xmin = 0
    xmax = 5
    ymin = -1
    ymax = 1
    nx = 50
    ny = 20
  []
[]

[Modules]
  [NavierStokesFV]
    simulation_type = 'transient'
    compressibility = 'incompressible'
    add_energy_equation = true

    density = 'rho'
    dynamic_viscosity = 'mu'
    thermal_conductivity = 'k'
    specific_heat = 'cp'

    initial_velocity = '${u_inlet} 1e-12 0'
    initial_pressure = 0.0
    initial_temperature = '${T_inlet}'

    inlet_boundaries = 'left'
    momentum_inlet_types = 'fixed-velocity'
    momentum_inlet_function = '${u_inlet} 0'
    energy_inlet_types = 'heatflux'
    energy_inlet_function = '${fparse u_inlet * rho * cp * T_inlet}'
    wall_boundaries = 'bottom top'
    momentum_wall_types = 'symmetry noslip'
    energy_wall_types = 'symmetry heatflux'
    energy_wall_function = '0'

    outlet_boundaries = 'right'
    momentum_outlet_types = 'fixed-pressure-zero-gradient'
    pressure_function = '${p_outlet}'

    ambient_convection_alpha = 'h_cv'
    ambient_temperature = 'T_solid'
  []
[]

[Materials]
  [constants]
    type = ADGenericFunctorMaterial
    prop_names = 'h_cv T_solid rho mu cp k dcp_dt'
    prop_values = '${h_fs} ${T_solid} ${rho} ${mu} ${cp} ${k} 0'
  []
[]

[Executioner]
  type = Transient
  solve_type = 'NEWTON'
  petsc_options_iname = '-pc_type -ksp_gmres_restart -sub_pc_type -sub_pc_factor_shift_type'
  petsc_options_value = 'asm      100                lu           NONZERO'
  line_search = 'none'
  nl_rel_tol = 7e-13
  dt = 0.4
  end_time = 0.8
[]

[Outputs]
  exodus = true
  csv = true
[]