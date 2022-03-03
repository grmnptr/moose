[StochasticTools]
[]

[Distributions]
  [mu1]
    type = Normal
    mean = '-0.4238230959932244'
    standard_deviation = 1e-5
  []
  [mu2]
    type = Normal
    mean = '2.2865480177108193'
    standard_deviation = 1e-5
  []
[]

[Samplers]
  [mc]
    type = MCT
    num_rows = 1
    distributions = 'mu1 mu2'
    # seed = 10
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
  [data]
    type = SamplerPostprocessorTransfer
    multi_app = sub
    sampler = mc
    to_vector_postprocessor = storage
    from_postprocessor = average
    # execute_on = TIMESTEP_BEGIN
    check_multiapp_execute_on = false
  []
[]

[VectorPostprocessors]
  [storage]
    type = StochasticResults
    # execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[Executioner]
  type = Transient
  num_steps = 1
[]

[Outputs]
  csv = true
[]
