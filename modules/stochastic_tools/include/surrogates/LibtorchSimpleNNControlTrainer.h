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

  /// Control reporter names
  std::vector<ReporterName> _control_names;

  /// The gathered data in a flattened form to be able to convert easily to torch::Tensor.
  std::vector<Real> _flattened_data;

  /// The gathered response in a flattened form to be able to convert easily to torch::Tensor.
  std::vector<Real> _flattened_response;

  /// Number of batches we want to prepare
  unsigned int _no_batches;

  /// Number of epochs for the training
  unsigned int _no_epocs;

  /// Number of hidden layers in the neural net
  unsigned int & _no_hidden_layers;

  /// Number of neurons within the hidden layers (the length of this vector
  /// should be the same as _no_hidden_layers)
  std::vector<unsigned int> & _no_neurons_per_layer;

  /// Name of the pytorch output file. This is used for loading and storing
  /// already existing data.
  std::string _filename;

  /// Switch indicating if an already existing neural net should be read from a
  /// file or not. This can be used to load existing torch files (from previous
  /// MOOSE or python runs for retraining and further manipulation)
  bool _read_from_file;

  /// The learning rate for the optimization algorithm
  Real _learning_rate;

#ifdef TORCH_ENABLED
  /// Pointer to the neural net object (initialized as null)
  std::shared_ptr<StochasticTools::LibtorchSimpleNeuralNet> & _control_nn;

  /// Pointer to the neural net object (initialized as null)
  std::shared_ptr<StochasticTools::LibtorchSimpleNeuralNet> & _emulator_nn;
#endif
};
