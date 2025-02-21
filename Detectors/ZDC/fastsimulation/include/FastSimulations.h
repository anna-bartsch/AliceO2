// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

///
/// @file   FastSimulations.h
/// @author SwirtaB
///

#ifndef O2_ZDC_FAST_SIMULATIONS_H
#define O2_ZDC_FAST_SIMULATIONS_H

#include <onnxruntime/core/session/onnxruntime_cxx_api.h>
#include <optional>

namespace o2::zdc::fastsim
{
/**
 * @brief Abstract class providing interface for various specialized implementations.
 *
 */
class NeuralFastSimulation
{
 public:
  NeuralFastSimulation(const std::string& modelPath,
                       Ort::SessionOptions sessionOptions,
                       OrtAllocatorType allocatorType,
                       OrtMemType memoryType);
  virtual ~NeuralFastSimulation() = default;

  /**
   * @brief Wrapper for converting raw input to Ort::Value.
   *
   * @param input flattened input data
   * @return true on success
   * @return false on failure
   */
  virtual bool setInput(std::vector<std::vector<float>>& input) = 0;
  /**
   * @brief Wraps Session.Run()
   *        Result should be stored as private member.
   */
  virtual void run() = 0;

  /// returns model output as const &.
  virtual const std::vector<Ort::Value>& getResult() = 0;

 protected:
  /// Sets models metadata (input/output layers names, inputs shape) in onnx session
  void setInputOutputData();
  /// Converts flattend input data to Ort::Value. Tensor shapes are taken from loaded model metadata.
  void setTensors(std::vector<std::vector<float>>& input);

  /// ONNX specific attributes
  /// User shoudn't has direct access to those in derived classes
  Ort::Env mEnv;
  Ort::Session mSession;
  Ort::AllocatorWithDefaultOptions mAllocator;
  Ort::MemoryInfo mMemoryInfo;

  /// Input/Output names and input shape
  std::vector<char*> mInputNames;
  std::vector<char*> mOutputNames;
  std::vector<std::vector<int64_t>> mInputShapes;

  /// Container for input tensors
  std::vector<Ort::Value> mInputTensors;
};

/**
 * @brief Derived class implementing interface for specific types of models.
 *
 */
class ConditionalModelSimulation : public NeuralFastSimulation
{
 public:
  ConditionalModelSimulation(const std::string& modelPath);
  ~ConditionalModelSimulation() override = default;

  /**
   * @brief Implements setInput
   *
   * @param input flattend input
   * @return true on success
   * @return false on failure
   */
  bool setInput(std::vector<std::vector<float>>& input) override;
  /**
   * @brief Implements run().
   *
   */
  void run() override;
  /**
   * @brief Returns single model output as const&.
   *        Returned vector is of size 1.
   *
   * @return std::vector<Ort::Value> model output
   */
  const std::vector<Ort::Value>& getResult() override;

 private:
  std::vector<Ort::Value> mModelOutput;
};

std::optional<std::pair<std::vector<float>, std::vector<float>>> loadScales(const std::string& path);

} // namespace o2::zdc::fastsim
#endif // ONNX_API_FAST_SIMULATIONS_H