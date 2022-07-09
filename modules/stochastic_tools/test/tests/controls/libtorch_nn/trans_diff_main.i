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
    reinit_after_solve = true
  []
[]

[Transfers]
  [nn_transfer]
    type = LibtorchNeuralNetTransfer
    multi_app = runner
    direction = 'to_multiapp'
    trainer_name = 'nn_trainer'
    control_name = src_control
  []
  [r_transfer]
    type = MultiAppReporterTransfer
    multi_app = runner
    direction = 'from_multiapp'
    to_reporters = 'results/T_max results/T_min results/control_value'
    from_reporters = 'T_reporter/T_max:value T_reporter/T_min:value T_reporter/control_value:value'
  []
[]

[Trainers]
  [nn_trainer]
    type = LibtorchSimpleNNControlTrainer
    sampler = dummy
    response_reporter = 'results/T_max results/T_min'
    response_constraints ='T_max_constraint T_min_constraint'
    control_reporter = 'results/control_value'

    # Parameters for the emulator neural net
    no_emulator_epocs = 4000
    no_emulator_batches = 5
    no_emulator_neurons_per_layer = '48 24'
    emulator_learning_rate = 0.0001

    # Parameters for the control neural net
    no_control_neurons_per_layer = '10 5'
    control_learning_rate = 0.0002
    no_control_epocs = 500
    no_control_loops = 1

    # General data
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
  # These are the constraints, t is T_max/T_min here (time hack)
  [T_max_constraint]
    type = ParsedFunction
    value = 'if(t > 2.5, 2.5, t)'
    # This constraint will penalize maximum temperatures below 2.5 C
  []
  [T_min_constraint]
    type = ParsedFunction
    value = 'if(t < 0, 0, t)'
    # This constraint will penalize minimum temperatures below 0 C
  []
[]

[Executioner]
  type = Transient
  num_steps = 200 # Number of training iterations
[]

[Outputs]
   csv = true
[]
