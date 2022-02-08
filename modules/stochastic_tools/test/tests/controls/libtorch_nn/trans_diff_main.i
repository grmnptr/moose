[StochasticTools]
[]

[Samplers]
  [dummy]
    type = CartesianProduct
    linear_space_items = '0 0.01 1'
  []
[]

[MultiApps]
  [runner]
    type = FullSolveMultiApp
    input_files = 'trans_diff_sub.i'
  []
[]

[Transfers]
  [r_transfer]
    type = MultiAppReporterTransfer
    multi_app = runner
    direction = 'from_multiapp'
    to_reporters = 'results/T_max results/T_min results/control_value'
    from_reporters = 'T_reporter/T_max:value T_reporter/T_min:value T_reporter/control_value:value'
  []
[]

[Trainers]
  [train]
    type = LibtorchSimpleNNControlTrainer
    sampler = dummy
    response_reporter = 'results/T_max results/T_min'
    response_constraints ='T_max_constraint T_min_constraint'
    control_reporter = 'results/control_value'
    no_epochs = 4000
    no_batches = 5
    no_hidden_layers = 2
    no_neurons_per_layer = '32 16'
    learning_rate = 0.0001
    control_learning_rate = 0.01
    filename = mynet.pt
    read_from_file = false
  []
[]

[Reporters]
  [results]
    type = ConstantReporter
    real_vector_names = 'T_max T_min control_value'
    real_vector_values = '0; 0; 0;'
    outputs = csv
    execute_on = timestep_begin
  []
[]

[Functions]
  [T_max_constraint]
    type = ParsedFunction
    value = 'if(t > 320, 320, t)'
  []
  [T_min_constraint]
    type = ParsedFunction
    value = 'if(t < 300, 300, t)'
  []
[]

# [Outputs]
#   csv = true
#   execute = asdasda
# []
