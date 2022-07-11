//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "LibtorchSimpleNNControlTrainer.h"
#include "Sampler.h"
#include "Function.h"

registerMooseObject("StochasticToolsApp", LibtorchSimpleNNControlTrainer);

InputParameters
LibtorchSimpleNNControlTrainer::validParams()
{
  InputParameters params = SurrogateTrainer::validParams();

  params.addClassDescription("Trains a simple neural network controller using libtorch.");

  params.addRequiredParam<std::vector<ReporterName>>("response_reporter", "Response reporters.");
  params.addRequiredParam<std::vector<ReporterName>>("control_reporter", "Control reporters.");
  params.addRequiredParam<std::vector<FunctionName>>(
      "response_constraints", "Constraints on the postprocessor values in the reporters.");

  params.addParam<unsigned int>("no_emulator_batches", 1, "Number of batches.");
  params.addParam<unsigned int>("no_emulator_epocs", 1, "Number of epochs.");
  params.addParam<unsigned int>("no_control_epocs", 1, "Number of epochs for the control.");
  params.addParam<unsigned int>("no_control_loops", 1, "Number of loops for training the control.");
  params.addParam<std::vector<unsigned int>>(
      "no_emulator_neurons_per_layer", std::vector<unsigned int>(), "Number of neurons per layer.");
  params.addParam<std::vector<unsigned int>>(
      "no_control_neurons_per_layer", std::vector<unsigned int>(), "Number of neurons per layer.");
  params.addParam<std::vector<std::string>>("emulator_activation_functions",
                                            std::vector<std::string>({"relu"}),
                                            "Number of neurons per layer.");
  params.addParam<std::vector<std::string>>("control_activation_functions",
                                            std::vector<std::string>({"relu"}),
                                            "Number of neurons per layer.");
  params.addParam<std::string>(
      "filename", "net.pt", "Filename used to output the neural net parameters.");
  params.addParam<bool>("read_from_file",
                        false,
                        "Switch to allow reading old trained neural nets for further training.");
  params.addParam<Real>("emulator_learning_rate", 0.001, "Learning rate (relaxation).");
  params.addParam<Real>("control_learning_rate", 0.001, "Control learning rate (relaxation).");

  params.addParam<unsigned int>(
      "seed", 11, "Random number generator seed for stochastic optimizers.");

  return params;
}

LibtorchSimpleNNControlTrainer::LibtorchSimpleNNControlTrainer(const InputParameters & parameters)
  : SurrogateTrainer(parameters),
    _response_names(getParam<std::vector<ReporterName>>("response_reporter")),
    _response_constraints(getParam<std::vector<FunctionName>>("response_constraints")),
    _control_names(getParam<std::vector<ReporterName>>("control_reporter")),
    _no_emulator_batches(getParam<unsigned int>("no_emulator_batches")),
    _no_emulator_epocs(getParam<unsigned int>("no_emulator_epocs")),
    _no_emulator_neurons_per_layer(
        getParam<std::vector<unsigned int>>("no_emulator_neurons_per_layer")),
    _emulator_activation_functions(
        getParam<std::vector<std::string>>("emulator_activation_functions")),
    _emulator_learning_rate(getParam<Real>("emulator_learning_rate")),
    _no_control_epocs(getParam<unsigned int>("no_control_epocs")),
    _no_control_loops(getParam<unsigned int>("no_control_loops")),
    _no_control_neurons_per_layer(
        getParam<std::vector<unsigned int>>("no_control_neurons_per_layer")),
    _control_activation_functions(
        getParam<std::vector<std::string>>("control_activation_functions")),
    _control_learning_rate(getParam<Real>("control_learning_rate")),
    _filename(getParam<std::string>("filename")),
    _read_from_file(getParam<bool>("read_from_file"))
{
  if (_response_names.size() == 0)
    mooseError("The number of reponses reporters should be more than 0!");

  if (_control_names.size() == 0)
    mooseError("The number of control reporters should be more than 0!");

  if (_response_names.size() != _response_constraints.size())
    paramError("response_constraints",
               "The number of responses is not equal to the number of response constraints!");

#ifdef LIBTORCH_ENABLED
  // Fixing the RNG seed to make sure every experiment is the same.
  // Otherwise sampling / stochastic gradient descent would be different.
  torch::manual_seed(getParam<unsigned int>("seed"));

  // Initializing and saving the control neural net so that the control can grab it right away
  _control_nn = std::make_shared<Moose::LibtorchArtificialNeuralNet>(_filename,
                                                                     2 * _response_names.size(),
                                                                     _control_names.size(),
                                                                     _no_control_neurons_per_layer,
                                                                     _control_activation_functions);
  torch::save(_control_nn, _control_nn->name());
#endif
}

void
LibtorchSimpleNNControlTrainer::preTrain()
{
}

void
LibtorchSimpleNNControlTrainer::train()
{
}

void
LibtorchSimpleNNControlTrainer::postTrain()
{
  // We collect the results from the resporters
  auto no_steps = getReporterValueByName<std::vector<Real>>(_response_names[0]).size();

  unsigned int no_rows = no_steps - 1;
  unsigned int no_cols = 2 * _response_names.size() + _control_names.size();
  unsigned int no_responses = _response_names.size();

  // Libtorch needs a 1Draw vector to be able to convert to torch::Tensors
  _flattened_data.resize(no_rows * no_cols);
  _flattened_response.resize(no_rows * _response_names.size());

  // Fill the 1D containers with the reporters that contain the quantities of interest
  for (unsigned int rep_i = 0; rep_i < _response_names.size(); ++rep_i)
  {
    const std::vector<Real> & data =
        getReporterValueByName<std::vector<Real>>(_response_names[rep_i]);
    for (unsigned int step_i = 0; step_i < data.size() - 1; ++step_i)
    {
      if (step_i == 0)
        _flattened_data[no_cols * step_i + rep_i] = data[step_i];
      else
        _flattened_data[no_cols * step_i + rep_i] = data[step_i - 1];
      _flattened_data[no_cols * step_i + rep_i + no_responses] = data[step_i];

      _flattened_response[no_responses * step_i + rep_i] = data[step_i + 1];
    }
  }

  // Fill the 1D containers with the reporters that contain the control values
  for (unsigned int rep_i = 0; rep_i < _control_names.size(); ++rep_i)
  {
    const std::vector<Real> & data =
        getReporterValueByName<std::vector<Real>>(_control_names[rep_i]);
    for (unsigned int step_i = 0; step_i < data.size() - 1; ++step_i)
      _flattened_data[no_cols * step_i + rep_i + 2 * no_responses] = data[step_i];
  }

  _communicator.allgather(_flattened_data);
  _communicator.allgather(_flattened_response);

  // We train the emulator to be able to emulate the system response
  trainEmulator();

  // We train the controller using the emulator to get a good control strategy
  trainController();
}

void
LibtorchSimpleNNControlTrainer::trainEmulator()
{

#ifdef LIBTORCH_ENABLED

  // This gets us the number of timesteps
  auto no_steps = getReporterValueByName<std::vector<Real>>(_response_names[0]).size();

  unsigned int n_rows = no_steps - 1;
  unsigned int n_cols = 2 * _response_names.size() + _control_names.size();
  unsigned int n_responses = _response_names.size();

  // The default data type in pytorch is float, while we use double in MOOSE.
  // Therefore, in some cases we have to convert Tensors to double.
  auto options = torch::TensorOptions().dtype(at::kDouble);
  torch::Tensor data_tensor =
      torch::from_blob(_flattened_data.data(), {n_rows, n_cols}, options).to(at::kDouble);
  torch::Tensor response_tensor =
      torch::from_blob(_flattened_response.data(), {n_rows, n_responses}, options).to(at::kDouble);

  // We create a custom data loader which can be used to select samples for the in
  // the training process. See the header file for the definition of this structure.
  MyData my_data(data_tensor, response_tensor);

  // We initialize a data_loader for the training part.
  unsigned int sample_per_batch = n_rows > _no_emulator_batches ? n_rows / _no_emulator_batches : 1;

  auto data_set = my_data.map(torch::data::transforms::Stack<>());
  auto data_loader = torch::data::make_data_loader<torch::data::samplers::SequentialSampler>(
      std::move(data_set), sample_per_batch);

  // We create a neural net
  _emulator_nn =
      std::make_shared<Moose::LibtorchArtificialNeuralNet>(_filename,
                                                           n_cols,
                                                           n_responses,
                                                           _no_emulator_neurons_per_layer,
                                                           _emulator_activation_functions);

  // Initialize the optimizer
  torch::optim::Adam optimizer(_emulator_nn->parameters(),
                               torch::optim::AdamOptions(_emulator_learning_rate));

  // Begin training loop
  for (size_t epoch = 1; epoch <= _no_emulator_epocs; ++epoch)
  {
    Real epoch_error = 0.0;
    for (auto & batch : *data_loader)
    {
      // Reset gradients
      optimizer.zero_grad();

      // Compute prediction
      torch::Tensor prediction = _emulator_nn->forward(batch.data);

      // Compute loss values using a MSE ( mean squared error)
      torch::Tensor loss = torch::mse_loss(prediction, batch.target);

      // Propagate error back
      loss.backward();

      // Use new gradients to update the parameters
      optimizer.step();

      epoch_error += loss.item<double>();
    }

    epoch_error = epoch_error / _no_emulator_batches;

    if (epoch % 10 == 0)
      _console << "Epoch: " << epoch << " | Loss: " << COLOR_GREEN << epoch_error << COLOR_DEFAULT
               << std::endl;
  }

#endif
}

void
LibtorchSimpleNNControlTrainer::trainController()
{

#ifdef LIBTORCH_ENABLED

  // We get the number of time steps
  auto no_steps = getReporterValueByName<std::vector<Real>>(_response_names[0]).size();

  unsigned int n_cols = 2 * _response_names.size();
  unsigned int n_responses = _response_names.size();
  unsigned int n_controls = _control_names.size();

  // We create the controller neural net
  _control_nn = std::make_shared<Moose::LibtorchArtificialNeuralNet>(
      _filename, n_cols, n_controls, _no_control_neurons_per_layer, _control_activation_functions);

  // Initialize the optimizer
  torch::optim::Adam optimizer(_control_nn->parameters(),
                               torch::optim::AdamOptions(_control_learning_rate));

  // We initialize the starting vector for the sweeps
  std::vector<Real> start_vector(&_flattened_data[0], &_flattened_data[n_cols]);
  auto options = torch::TensorOptions().dtype(at::kDouble);
  torch::Tensor input = torch::from_blob(start_vector.data(), {1, n_cols}, options).to(at::kDouble);

  // Begin controller training loop
  for (unsigned int loop_i = 0; loop_i < _no_control_loops; ++loop_i)
  {
    // We loop through the timesteps simulating the system using the emulator
    // neural net
    for (unsigned int step_i = 1; step_i <= no_steps - 1; ++step_i)
    {
      Real epoch_error = 100.0;
      unsigned int epoch_counter = 0;
      std::vector<Real> converted_prediction;

      // We iterate in each timestep until we get a good neural net control
      while (epoch_counter < _no_control_epocs && epoch_error > 1e-6)
      {
        epoch_counter += 1;

        optimizer.zero_grad();

        // We predict the controller values based on the previous two timesteps
        torch::Tensor control_prediction = _control_nn->forward(input);

        // We add the controller values to the input vectors
        torch::Tensor extended_input = torch::cat({input, control_prediction}, -1).to(at::kDouble);

        // We use the emulator to predict the values of the quantities interest in the
        // next timestep
        torch::Tensor value_prediction = _emulator_nn->forward(extended_input);

        converted_prediction =
            std::vector<Real>({value_prediction.data_ptr<Real>(),
                               value_prediction.data_ptr<Real>() + value_prediction.size(1)});

        // Evaluate the contraints for each response
        std::vector<Real> constraints;
        Point dummy;
        for (unsigned int resp_i = 0; resp_i < n_responses; ++resp_i)
        {
          const Function & constraint_function(getFunctionByName(_response_constraints[resp_i]));
          auto constraint = constraint_function.value(converted_prediction[resp_i], dummy);
          constraints.push_back(constraint);
        }

        // Build the constraint tensor for the loss calculation
        torch::Tensor constraint_tensor =
            torch::from_blob(constraints.data(), {1, n_responses}, options).to(at::kDouble);

        // Compute loss values using a MSE ( mean squared error)
        torch::Tensor loss = torch::mse_loss(value_prediction, constraint_tensor);

        // Propagate error back
        loss.backward();

        // Use new gradients to update the parameters
        optimizer.step();

        epoch_error = loss.item<double>();
      }

      _console << "Controller training step: " << step_i << " (" << epoch_counter
               << ") | Loss: " << COLOR_GREEN << epoch_error << COLOR_DEFAULT << std::endl;

      // Build a new input tensor using the outputs
      for (unsigned int resp_i = 0; resp_i < n_responses; ++resp_i)
      {
        start_vector[resp_i] = start_vector[resp_i + n_responses];
        start_vector[resp_i + n_responses] = converted_prediction[resp_i];
      }
      input = torch::from_blob(start_vector.data(), {1, n_cols}, options).to(at::kDouble);
    }
  }

  // Save the controller neural net so our controller can read it
  torch::save(_control_nn, _control_nn->name());

#endif
}
