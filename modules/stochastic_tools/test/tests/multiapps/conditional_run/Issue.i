[StochasticTools]
[]

[Distributions]
  [left_bc]
    type = Uniform
    lower_bound = 0
    upper_bound = 0.1
  []
  [right_bc]
    type = Uniform
    lower_bound = 0.9
    upper_bound = 1.0
  []
[]

[Samplers]
  [mc]
    type = MCT
    num_rows = 2 # 1
    distributions = 'left_bc right_bc'
    seed = 5
    execute_on = PRE_MULTIAPP_SETUP
  []
[]

[MultiApps]
  [sub]
    type = SamplerFullSolveMultiApp
    sampler = mc
    input_files = 'sub.i' 
    mode = batch-reset
    should_run_reporter = conditional/need_sample
    execute_on = TIMESTEP_END
  []
[]

[Controls]
  [cmdline]
    type = MultiAppCommandLineControl
    multi_app = sub
    sampler = mc
    param_names = 'BCs/left/value BCs/right/value'
  []
[]

[Transfers]
  [reporter_transfer]
    type = SamplerReporterTransfer
    from_reporter = 'average/value'
    stochastic_reporter = 'conditional'
    from_multi_app = sub
    sampler = mc
  []
[]

[Reporters]
  [conditional]
    type = Issue # ConditionalSampleReporter
    sampler = mc
    default_value = 1
    function = 'val >= t'
    sampler_vars = 'val'
    sampler_var_indices = '0'
    parallel_type = ROOT
    execute_on = 'initial timestep_begin'
  []
[]

[Executioner]
  type = Transient
  num_steps = 2
[]

[Outputs]
  execute_on = timestep_end
  [out]
    type = JSON
    execute_system_information_on = none
  []
[]
