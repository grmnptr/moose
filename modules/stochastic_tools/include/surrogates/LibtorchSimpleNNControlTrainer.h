//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#ifdef TORCH_ENABLED
#include <torch/torch.h>
#include "LibtorchSimpleNeuralNet.h"
#endif

#include "libmesh/utility.h"
#include "SurrogateTrainer.h"

class LibtorchSimpleNNControlTrainer : public SurrogateTrainer
{
public:
  static InputParameters validParams();

  LibtorchSimpleNNControlTrainer(const InputParameters & parameters);

  virtual void preTrain() override;

  virtual void train() override;

  virtual void postTrain() override;

  void trainEmulator();

  void trainController();

#ifdef TORCH_ENABLED
  const std::shared_ptr<StochasticTools::LibtorchSimpleNeuralNet> & controlNeuralNet() const
  {
    return _control_nn;
  }
#endif

protected:
#ifdef TORCH_ENABLED

  // A custom strcuture which is used to organize data foor the training of
  // torch-based neural nets.
  struct MyData : torch::data::datasets::Dataset<MyData>
  {
  public:
    MyData(torch::Tensor dt, torch::Tensor rt) : _data_tensor(dt), _response_tensor(rt) {}
    torch::data::Example<> get(size_t index) override
    {
      return {_data_tensor[index], _response_tensor[index]};
    }

    torch::optional<size_t> size() const override { return _response_tensor.sizes()[0]; }

  private:
    torch::Tensor _data_tensor;
    torch::Tensor _response_tensor;
  };

#endif

private:
  /// Response reporter names
  std::vector<ReporterName> _response_names;

  /// Names of the functions describing the response constraints
  std::vector<FunctionName> _response_constraints;

  /// Control reporter names
  std::vector<ReporterName> _control_names;

  /// The gathered data in a flattened form to be able to convert easily to torch::Tensor.
  std::vector<Real> _flattened_data;

  /// The gathered response in a flattened form to be able to convert easily to torch::Tensor.
  std::vector<Real> _flattened_response;

  /// Number of batches we want to prepare for the emulator
  unsigned int _no_emulator_batches;

  /// Number of epochs for the training of the emulator
  unsigned int _no_emulator_epocs;

  /// Number of hidden layers in the emulator neural net
  unsigned int & _no_emulator_hidden_layers;

  /// Number of neurons within the hidden layers i nthe emulator neural net
  std::vector<unsigned int> & _no_emulator_neurons_per_layer;

  /// The learning rate for the optimization algorithm in the meulator
  Real _emulator_learning_rate;

  /// Number of control epochs for the training
  unsigned int _no_control_epocs;

  /// Number of control loops for the training
  unsigned int _no_control_loops;

  /// Number of hidden layers in the neural net
  const unsigned int & _no_control_hidden_layers;

  /// Number of neurons within the hidden layers in the control neural net
  const std::vector<unsigned int> & _no_control_neurons_per_layer;

  /// The control learning rate for the optimization algorithm
  Real _control_learning_rate;

  /// Name of the pytorch output file. This is used for loading and storing
  /// already existing data.
  std::string _filename;

  /// Switch indicating if an already existing neural net should be read from a
  /// file or not. This can be used to load existing torch files (from previous
  /// MOOSE or python runs for retraining and further manipulation)
  bool _read_from_file;

#ifdef TORCH_ENABLED
  /// Pointer to the emulator neural net object (initialized as null)
  std::shared_ptr<StochasticTools::LibtorchSimpleNeuralNet> & _control_nn;

  /// Pointer to the controller neural net object (initialized as null)
  std::shared_ptr<StochasticTools::LibtorchSimpleNeuralNet> & _emulator_nn;
#endif
};
