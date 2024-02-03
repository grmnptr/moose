[StochasticTools]
[]

[Samplers]
  [sample]
    type = CartesianProduct
    linear_space_items = '0.0 0.1 10'
  []
[]

[MultiApps]
  [sub]
    type = SamplerFullSolveMultiApp
    input_files = sub_modified.i
    sampler = sample
    mode = batch-reset
  []
[]

[Transfers]
  [data]
    type = SamplerReporterTransfer
    from_multi_app = sub
    sampler = sample
    stochastic_reporter = storage
    from_reporter = 'pp/value constant/str'
  []
  [data_to]
    type = SamplerReporterTransfer
    to_multi_app = sub
    source_reporter = constant/myreal
    to_reporter = constant_receiver/receiver
    sampler = sample
    stochastic_reporter = storage
  []
  [runner]
    type = SamplerParameterTransfer
    to_multi_app = sub
    sampler = sample
    parameters = 'BCs/left/value'
  []
[]

[Reporters]
  [storage]
    type = StochasticReporter
    parallel_type = ROOT
  []
  [constant]
    type = ConstantReporter
    real_names = "myreal"
    real_values = "0.5"
  []
[]

[Outputs]
  [out]
    type = JSON
    execute_on = timestep_end
  []
[]
