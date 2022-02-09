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
#    reset_apps = '0'
#    reset_time = 1
     no_backup_and_restore = true
  []
[]

[Transfers]
#  [nn_transfer]
#    type = LibtorchNeuralNetTransfer
#    multi_app = runner
#    direction = 'to_multiapp'
#    trainer_name = 'nn_trainer'
#    control_name = bd_control
#  []
  [r_transfer]
    type = MultiAppReporterTransfer
    multi_app = runner
    direction = 'from_multiapp'
    to_reporters = 'results/T_max results/T_min results/control_value'
    from_reporters = 'T_reporter/T_max:value T_reporter/T_min:value T_reporter/control_value:value'
  []
[]

# [Trainers]
#  [nn_trainer]
#    type = LibtorchSimpleNNControlTrainer
#    sampler = dummy
#    response_reporter = 'results/T_max results/T_min'
#    response_constraints ='T_max_constraint T_min_constraint'
#    control_reporter = 'results/control_value'
#    no_epocs = 8000
#    no_control_epocs = 400
#    no_control_loops = 2
#    no_batches = 5
#    no_hidden_layers = 3
#    no_neurons_per_layer = '64 32 16'
#    learning_rate = 0.0001
#    control_learning_rate = 0.001
#    filename = mynet.pt
#    read_from_file = false
#  []
#[]

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
    value = 'if(t > 310, 310, t)'
  []
  [T_min_constraint]
    type = ParsedFunction
    value = 'if(t < 300, 300, t)'
  []
[]

[Executioner]
  type = Transient
  num_steps = 2
[]

# [Outputs]
#   csv = true
# []
