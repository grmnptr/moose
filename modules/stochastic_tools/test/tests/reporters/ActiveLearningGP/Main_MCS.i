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
    seed = 10
  []
[]

[MultiApps]
  [sub]
    type = SamplerFullSolveMultiApp
    sampler = mc
    input_files = 'Sub.i'
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
[]

[Executioner]
  type = Transient
  num_steps = 1000
[]

[Outputs]
  [out]
    type = JSON
    execute_system_information_on = none
  []
[]
