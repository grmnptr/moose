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
  # [mu1]
  #   type = Uniform
  #   lower_bound = 1
  #   upper_bound = 10
  # []
  # [mu2]
  #   type = Uniform
  #   lower_bound = 9000
  #   upper_bound = 11000
  # []
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
    input_files = 'Sub.i' # 'Sub_new.i'
    mode = batch-reset
    should_run_reporter = conditional/need_sample
  []
[]

[Transfers]
  [param]
    type = SamplerParameterTransfer
    multi_app = sub
    sampler = mc
    parameters = 'BCs/left/value BCs/right/value' # 'Materials/conductivity/prop_values Kernels/source/value' #
    to_control = 'stochastic'
  []
  [reporter_transfer]
    type = SamplerReporterTransfer
    from_reporter = 'average/value'
    stochastic_reporter = 'conditional' # constant
    multi_app = sub
    sampler = mc
  []
[]

[Reporters]
  # [constant]
  #   type = StochasticReporter
  #   # execute_on = 'initial timestep_begin' #
  # []
  [conditional]
    type = AL_ADAM # ActiveLearningGP # 
    sampler = mc
    # output_value = constant/reporter_transfer:average:value
    parallel_type = ROOT
    execute_on = 'initial timestep_begin'
    covariance_function = 'covar'             #Choose a squared exponential for the kernel
    standardize_params = 'true'               #Center and scale the training params
    standardize_data = 'true'                 #Center and scale the training data
    tao_options = '-tao_max_it 100000 -tao_max_funcs 100000 -tao_fatol 1e-10' #  -tao_cg_type hs -tao_bncg_type gd
    tune_parameters = 'signal_variance length_factor' #
    tuning_min = '1e-3 1e-3'
    tuning_max = '1000 1000'
    show_tao = 'true'
  []
[]

[Covariance]
  [covar]
    type= SquaredExponentialCovariance # MaternHalfIntCovariance #
    # p = 1.0
    signal_variance = 1.0                 #Use a signal variance of 1 in the kernel
    noise_variance = 1e-4                    #A small amount of noise can help with numerical stability
    length_factor = '1.0 1.0'         #Select a length factor for each parameter (k and q)
  []
[]

[Executioner]
  type = Transient
  num_steps = 13
[]

[Outputs]
  # execute_on = timestep_end
  [out]
    type = JSON
    execute_system_information_on = none
  []
[]
