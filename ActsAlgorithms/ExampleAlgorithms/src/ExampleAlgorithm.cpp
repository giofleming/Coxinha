#include "ExampleAlgorithm.hpp"

using namespace ActsExamples;

namespace AliceActsTrk
{

ExampleAlgorithm::ExampleAlgorithm(ExampleAlgorithm::Config cfg,
                                   std::unique_ptr<const Acts::Logger> logger)
  : IAlgorithm("ExampleAlgorithm", std::move(logger)), m_cfg(std::move(cfg))
{

  std::cout << "PF:: Constructed ExampleAlgorithm" << std::endl;
}

ProcessCode ExampleAlgorithm::execute(const AlgorithmContext& ctx) const
{

  std::cout << "PF:: Executing ExampleAlgorithm" << std::endl;

  return ProcessCode::SUCCESS;
}

} // namespace AliceActsTrk
