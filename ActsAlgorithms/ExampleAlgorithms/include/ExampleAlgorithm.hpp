#pragma once

#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"

#include <memory>
#include <string>
#include <utility>
#include <vector>

using namespace ActsExamples;

namespace AliceActsTrk
{

class ExampleAlgorithm final : public IAlgorithm
{

 public:
  struct Config {
    std::string testString;
  };

  ExampleAlgorithm(Config cfg, std::unique_ptr<const Acts::Logger> logger = nullptr);

  ActsExamples::ProcessCode execute(const AlgorithmContext& ctx) const final;

  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;

}; // AliceActsTrk::ExampleAlgorithm

} // namespace AliceActsTrk
