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
    # seed = 1
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
    type = ActiveLearningGP
    sampler = mc
    output_value = constant/reporter_transfer:average:value
    parallel_type = ROOT
    execute_on = 'initial timestep_begin'
    covariance_function = 'covar'             #Choose a squared exponential for the kernel
    standardize_params = 'true'               #Center and scale the training params
    standardize_data = 'true'                 #Center and scale the training data
    tao_options = '-tao_bncg_type ssml_bfgs'
    tune_parameters = 'signal_variance length_factor'
    tuning_min = ' 1e-9 1e-9'
    tuning_max = ' 1e16  1e16'
    show_tao = 'true'
  []
[]

[Covariance]
  [covar]
    type=SquaredExponentialCovariance
    signal_variance = 0.01                       #Use a signal variance of 1 in the kernel
    noise_variance = 1e-3                     #A small amount of noise can help with numerical stability
    length_factor = '0.01 0.01'         #Select a length factor for each parameter (k and q)
  []
[]

[Executioner]
  type = Transient
  num_steps = 12
[]

[Outputs]
  # execute_on = timestep_end
  [out]
    type = JSON
    execute_system_information_on = none
  []
[]
