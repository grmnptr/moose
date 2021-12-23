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
  [sample]
    type = MCT
    num_rows = 1
    distributions = 'mu1 mu2'
  []
[]

[MultiApps]
  [sub]
    type = SamplerFullSolveMultiApp
    input_files = sub.i
    sampler = sample
    mode = batch-reset
    should_run_reporter = al_blin/need_sample
  []
[]

[Transfers]
  [param]
    type = SamplerParameterTransfer
    multi_app = sub
    sampler = sample
    parameters = 'BCs/left/value BCs/right/value'
    to_control = 'stochastic'
  []
  [data]
    type = SamplerReporterTransfer
    multi_app = sub
    sampler = sample
    from_reporter = 'average/value'
    stochastic_reporter = al_blin
  []
  # [data_inp]
  #   type = SamplerData
  #   sampler = sample
  #   execute_on = 'initial timestep_end'
  # []
[]

[Reporters]
  [constant]
    type = StochasticReporter
  []
  [al_blin]
    type = ActiveLearningBayesLinReg
    sampler = sample
    parallel_type = ROOT
    execute_on = 'initial timestep_begin' #
  []
[]

[Executioner]
  type = Transient
  num_steps = 10
[]

[Outputs]
  [out]
    type = JSON
    execute_system_information_on = NONE
  []
[]
