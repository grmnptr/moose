[StochasticTools]
[]

[Distributions]
  [k_dist]
    type = Uniform
    lower_bound = 5
    upper_bound = 20
  []
  [q_dist]
    type = Uniform
    lower_bound = 7000
    upper_bound = 13000
  []
  [L_dist]
    type = Uniform
    lower_bound = 0.05
    upper_bound = 0.15
  []
  [Tinf_dist]
    type = Uniform
    lower_bound = 250
    upper_bound = 350
  []
[]

[Samplers]
  [mc]
    type = MCT
    num_rows = 1
    distributions = 'k_dist q_dist L_dist Tinf_dist'
    flag_sample = 'conditional/flag_sample'
    seed = 10
    execute_on = PRE_MULTIAPP_SETUP
  []
[]

[MultiApps]
  [sub]
    type = SamplerFullSolveMultiApp
    sampler = mc
    input_files = 'Sub4D.i' # 'Sub_new.i'
    mode = batch-reset
    should_run_reporter = conditional/need_sample
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
  #   parameters = 'Materials/conductivity/prop_values Kernels/source/value Mesh/xmax BCs/right/value'
  #   to_control = 'stochastic'
  # []
  [reporter_transfer]
    type = SamplerReporterTransfer
    from_reporter = 'avg/value'
    stochastic_reporter = 'conditional' # constant
    multi_app = sub
    sampler = mc
  []
[]

[Reporters]
  [constant]
     type = ConstantReporter
     execute_on = 'final'
  []
  [conditional]
    type =  ActiveLearningGP # AL_ADAM #
    sampler = mc
    parallel_type = ROOT
    execute_on = 'initial timestep_begin'
    covariance_function = 'covar'             #Choose a squared exponential for the kernel
    standardize_params = 'true'               #Center and scale the training params
    standardize_data = 'true'                 #Center and scale the training data
    tao_options = '-tao_max_it 100000 -tao_max_funcs 100000 -tao_fatol 1e-6 -tao_bncg_type prp -tao_bncg_xi=0.7'
    tune_parameters = 'signal_variance length_factor' #
    tuning_min = '1e-3 1e-3'
    tuning_max = '10000 10000'
    show_tao = 'true'
    flag_sample = 'flag_sample'
    threshold = 349.345
    outputs = cout
  []
  [acc]
    type = AccumulateReporter
    reporters = 'conditional/value'
    outputs = aout
  []
[]

[Covariance]
  [covar]
    type= SquaredExponentialCovariance # MaternHalfIntCovariance #
    # p = 1.0
    signal_variance = 1.0                 #Use a signal variance of 1 in the kernel
    noise_variance = 1e-4                    #A small amount of noise can help with numerical stability
    length_factor = '1.0 1.0 1.0 1.0'         #Select a length factor for each parameter (k and q)
  []
[]

[Executioner]
  type = Transient
  num_steps = 10000
[]

[Outputs]
  [cout]
    type = CSV
    execute_on = timestep_end
  []
  [aout]
    type = CSV
    execute_on = final
  []
[]
