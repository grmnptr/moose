[StochasticTools]
[]

[Distributions]
  [mu1]
    type = Normal
    mean = 0.0
    standard_deviation = 0.5
  []
  [mu2]
    type = Normal
    mean = 1
    standard_deviation = 0.5
  []
[]

[Samplers]
  [mc]
    type = MCT
    num_rows = 1
    distributions = 'mu1 mu2'
    seed = 1
  []
[]

[MultiApps]
  [sub]
    type = SamplerFullSolveMultiApp
    sampler = mc
    input_files = 'Sub.i'
    mode = batch-reset
    should_run_reporter = conditional/need_sample
  []
[]

[Transfers]
  [param]
    type = SamplerParameterTransfer
    multi_app = sub
    sampler = mc
    parameters = 'BCs/left/value BCs/right/value'
    to_control = 'stochastic'
  []
  [reporter_transfer]
    type = SamplerReporterTransfer
    from_reporter = 'average/value'
    stochastic_reporter = 'constant'
    multi_app = sub
    sampler = mc
  []
[]

[Reporters]
  [constant]
    type = StochasticReporter
  []
  [conditional]
    type = ActiveLearningBayesGaussBasis
    sampler = mc
    output_value = constant/reporter_transfer:average:value
    parallel_type = ROOT
    execute_on = 'initial timestep_begin'
  []
[]

[Executioner]
  type = Transient
  num_steps = 200
[]

[Outputs]
  # execute_on = timestep_end
  [out]
    type = JSON
    execute_system_information_on = none
  []
[]
