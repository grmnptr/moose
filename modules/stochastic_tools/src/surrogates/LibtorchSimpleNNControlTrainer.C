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

registerMooseObject("StochasticToolsApp", LibtorchSimpleNNControlTrainer);

InputParameters
LibtorchSimpleNNControlTrainer::validParams()
{
  InputParameters params = SurrogateTrainer::validParams();

  params.addClassDescription("Trains a simple neural network using libtorch.");

  params.addRequiredParam<std::vector<ReporterName>>("response_reporter", "Response reporters.");
  params.addRequiredParam<std::vector<ReporterName>>("control_reporter", "Control reporters.");

  params.addParam<unsigned int>("no_batches", 1, "Number of batches.");
  params.addParam<unsigned int>("no_epochs", 1, "Number of epochs.");
  params.addParam<unsigned int>("no_hidden_layers", 0, "Number of hidden layers.");
  params.addParam<std::vector<unsigned int>>(
      "no_neurons_per_layer", std::vector<unsigned int>(), "Number of neurons per layer.");
  params.addParam<std::string>(
      "filename", "net.pt", "Filename used to output the neural net parameters.");
  params.addParam<bool>("read_from_file",
                        false,
                        "Switch to allow reading old trained neural nets for further training.");
  params.addParam<Real>("learning_rate", 0.001, "Learning rate (relaxation).");
  params.addParam<unsigned int>(
      "seed", 11, "Random number generator seed for stochastic optimizers.");

  return params;
}

LibtorchSimpleNNControlTrainer::LibtorchSimpleNNControlTrainer(const InputParameters & parameters)
  : SurrogateTrainer(parameters),
    _response_names(getParam<std::vector<ReporterName>>("response_reporter")),
    _control_names(getParam<std::vector<ReporterName>>("control_reporter")),
    _no_batches(getParam<unsigned int>("no_batches")),
    _no_epocs(getParam<unsigned int>("no_epochs")),
    _no_hidden_layers(declareModelData<unsigned int>("no_hidden_layers")),
    _no_neurons_per_layer(declareModelData<std::vector<unsigned int>>("no_neurons_per_layer")),
    _filename(getParam<std::string>("filename")),
    _read_from_file(getParam<bool>("read_from_file")),
    _learning_rate(getParam<Real>("learning_rate"))
#ifdef TORCH_ENABLED
    ,
    _control_nn(
        declareModelData<std::shared_ptr<StochasticTools::LibtorchSimpleNeuralNet>>("control_nn")),
    _emulator_nn(
        declareModelData<std::shared_ptr<StochasticTools::LibtorchSimpleNeuralNet>>("emulator_nn"))
#endif
{
  if (_no_hidden_layers != _no_neurons_per_layer.size())
    mooseError("The number of layers are not the same!");

  if (_response_names.size() == 0)
    mooseError("The number of reponses reporters should be more than 0!");

  if (_control_names.size() == 0)
    mooseError("The number of control reporters should be more than 0!");

  _no_hidden_layers = getParam<unsigned int>("no_hidden_layers");
  _no_neurons_per_layer = getParam<std::vector<unsigned int>>("no_neurons_per_layer");
  _filename = getParam<std::string>("filename");

#ifdef TORCH_ENABLED
  // Fixing the RNG seed to make sure every experiment is the same.
  // Otherwise sampling / stochastic gradient descent would be different.
  torch::manual_seed(getParam<unsigned int>("seed"));
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
  auto no_steps = getReporterValueByName<std::vector<Real>>(_response_names[0]).size();

  unsigned int no_rows = no_steps - 1;
  unsigned int no_cols = 2 * _response_names.size() + _control_names.size();
  _flattened_data.resize(no_rows * no_cols);
  _flattened_response.resize(no_rows * _response_names.size());

  for (unsigned int rep_i = 0; rep_i < _response_names.size(); ++rep_i)
  {
    const std::vector<Real> & data =
        getReporterValueByName<std::vector<Real>>(_response_names[rep_i]);
    for (unsigned int step_i = 0; step_i < data.size() - 1; ++step_i)
    {
      if (step_i == 0)
        _flattened_data[no_cols * step_i + 2 * rep_i] = data[step_i];
      else
        _flattened_data[no_cols * step_i + 2 * rep_i] = data[step_i - 1];
      _flattened_data[no_cols * step_i + 2 * rep_i + 1] = data[step_i];

      _flattened_response[_response_names.size() * step_i + rep_i] = data[step_i + 1];
    }
  }

  for (unsigned int rep_i = 0; rep_i < _control_names.size(); ++rep_i)
  {
    const std::vector<Real> & data =
        getReporterValueByName<std::vector<Real>>(_control_names[rep_i]);
    for (unsigned int step_i = 0; step_i < data.size() - 1; ++step_i)
      _flattened_data[no_cols * step_i + rep_i + 2 * _response_names.size()] = data[step_i];
  }

  _communicator.allgather(_flattened_data);
  _communicator.allgather(_flattened_response);

  trainEmulator();
}

void
LibtorchSimpleNNControlTrainer::trainEmulator()
{

#ifdef TORCH_ENABLED

  // Then, we create and load our Tensors
  auto no_steps = getReporterValueByName<std::vector<Real>>(_response_names[0]).size();

  unsigned int n_rows = no_steps - 1;
  unsigned int n_cols = 2 * _response_names.size() + _control_names.size();
  unsigned int n_cols_response = _response_names.size();

  // The default data type in pytorch is float, while we use double in MOOSE.
  // Therefore, in some cases we have to convert Tensors to double.
  auto options = torch::TensorOptions().dtype(at::kDouble);
  torch::Tensor data_tensor =
      torch::from_blob(_flattened_data.data(), {n_rows, n_cols}, options).to(at::kDouble);
  torch::Tensor response_tensor =
      torch::from_blob(_flattened_response.data(), {n_rows, n_cols_response}, options)
          .to(at::kDouble);

  // We create a custom data loader which can be used to select samples for the in
  // the training process. See the header file for the definition of this structure.
  MyData my_data(data_tensor, response_tensor);

  // We initialize a data_loader for the training part.
  unsigned int sample_per_batch = n_rows > _no_batches ? n_rows / _no_batches : 1;

  auto data_set = my_data.map(torch::data::transforms::Stack<>());
  auto data_loader = torch::data::make_data_loader<torch::data::samplers::SequentialSampler>(
      std::move(data_set), sample_per_batch);

  // We create a neural net (for the definition of the net see the header file)
  _emulator_nn = std::make_shared<StochasticTools::LibtorchSimpleNeuralNet>(
      _filename, n_cols, _no_hidden_layers, _no_neurons_per_layer, n_cols_response);

  // Initialize the optimizer
  torch::optim::Adam optimizer(_emulator_nn->parameters(),
                               torch::optim::AdamOptions(_learning_rate));

  if (_read_from_file)
    try
    {
      torch::load(_emulator_nn, _filename);
      _console << "Loaded requested .pt file." << std::endl;
    }
    catch (...)
    {
      mooseError("The requested pytorch file could not be loaded.");
    }

  // Begin training loop
  for (size_t epoch = 1; epoch <= _no_epocs; ++epoch)
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

    epoch_error = epoch_error / _no_batches;

    if (epoch % 10 == 0)
      _console << "Epoch: " << epoch << " | Loss: " << COLOR_GREEN << epoch_error << COLOR_DEFAULT
               << std::endl;
  }

#endif
}
