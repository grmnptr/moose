[StochasticTools]
[]

[Distributions]
  [k_dist]
    type = Uniform
    lower_bound = 0
    upper_bound = 20
  []
  [q_dist]
    type = Uniform
    lower_bound = 7000
    upper_bound = 13000
  []
  [L_dist]
    type = Uniform
    lower_bound = 0.0
    upper_bound = 0.1
  []
  [Tinf_dist]
    type = Uniform
    lower_bound = 270
    upper_bound = 330
  []
[]

[Samplers]
  [mc]
    type = MCT
    num_rows = 1
    distributions = 'k_dist q_dist L_dist Tinf_dist'
    seed = 10
    execute_on = PRE_MULTIAPP_SETUP
  []
[]

[MultiApps]
  [sub]
    type = SamplerFullSolveMultiApp
    sampler = mc
    input_files = 'Sub4D.i'
  []
[]

[Controls]
  [cmdline]
    type = MultiAppCommandLineControl
    multi_app = sub
    sampler = mc
    param_names = 'Materials/conductivity/prop_values Kernels/source/value Mesh/xmax BCs/right/value'
  []
[]

[Transfers]
  # [param]
  #   type = SamplerParameterTransfer
  #   multi_app = sub
  #   sampler = mc
  #   parameters = 'BCs/left/value BCs/right/value'
  #   to_control = 'stochastic'
  # []
  [reporter_transfer]
    type = SamplerReporterTransfer
    from_reporter = 'avg/value'
    stochastic_reporter = 'constant'
    multi_app = sub
    sampler = mc
  []
[]

[Reporters]
  [constant]
    type = StochasticReporter
  []
[]

[Executioner]
  type = Transient
  num_steps = 100
[]

[Outputs]
  [out]
    type = JSON
    execute_system_information_on = none
  []
[]
