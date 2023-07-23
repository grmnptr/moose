[StochasticTools]
[]

[Samplers]
  [vpp]
    type = VectorPostprocessorSampler
    vectors_names = 'csv_reader/data'
    execute_on = 'initial timestep_end'
  []
[]

[MultiApps]
  [runner]
    type = SamplerTransientMultiApp
    sampler = vpp
    input_files = 'diffusion_time.i'
    mode = batch-restore
  []
[]

[Transfers]
  [parameters]
    type = SamplerParameterTransfer
    to_multi_app = runner
    sampler = vpp
    parameters = 'Kernels/source/value'
  []
[]

[VectorPostprocessors]
  [csv_reader]
    type = CSVReader
    csv_file = 'data.csv'
  []
[]

[Executioner]
  type = Transient
  num_steps = 4
  dt = 0.25
[]

[Outputs]
  execute_on = timestep_end
  [out]
    type = JSON
  []
[]
