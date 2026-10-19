// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <Eigen/Core>
#include <TH2F.h>

namespace AliceActsTrk
{

class ResidualAccumulator
{
 public:
  static constexpr int kSize = 5;

  using Vector = Eigen::Matrix<double, kSize, 1>;
  using Matrix = Eigen::Matrix<double, kSize, kSize>;

  struct GaussianCoreResult {
    Vector mean = Vector::Zero();
    Vector sigma = Vector::Zero();
    Matrix covariance = Matrix::Zero();

    std::size_t coreEntries = 0; ///< Number of samples in the common fitted core.
    bool valid = false;          ///< True when the fitted covariance is finite and usable.
  };

  explicit ResidualAccumulator() = default;

  /// Add one residual sample and update the running ordinary covariance.
  /// Non-finite samples are rejected and do not affect the accumulator.
  bool add(const Vector& residual)
  {
    if (!residual.allFinite()) {
      ++mRejected;
      return false;
    }

    // Keep the individual samples because they are needed
    // for the Gaussian-core determination.
    mSamples.push_back(residual);

    // Multivariate Welford update for the ordinary covariance.
    ++mEntries;

    const Vector delta = residual - mMean;
    mMean += delta / static_cast<double>(mEntries);

    const Vector deltaAfterUpdate = residual - mMean;
    mM2 += delta * deltaAfterUpdate.transpose();

    mLastResidual = residual;

    return true;
  }

  [[nodiscard]]
  std::size_t entries() const
  {
    return mEntries;
  }

  [[nodiscard]]
  std::size_t rejected() const
  {
    return mRejected;
  }

  [[nodiscard]]
  const Vector& mean() const
  {
    return mMean;
  }

  [[nodiscard]]
  Matrix covariance() const
  {
    if (mEntries < 2) {
      return Matrix::Zero();
    }

    Matrix cov = mM2 / static_cast<double>(mEntries - 1);

    // Remove tiny numerical asymmetries.
    cov = 0.5 * (cov + cov.transpose());

    return cov;
  }

  /// Determine the covariance matrix from the Gaussian core
  /// of the accumulated residual distribution.
  ///
  /// For every parameter:
  ///   - fit the 1D residual distribution with a Gaussian;
  ///   - extract the core mean and sigma.
  ///
  /// A common multidimensional core sample is then selected and used
  /// to determine the correlations.
  ///
  /// The final covariance is
  ///
  ///   C_ij = rho_ij(core) * sigma_i(fit) * sigma_j(fit)
  ///
  /// and therefore
  ///
  ///   C_ii = sigma_i(fit)^2.
  ///
  /// The result is invalid when there are too few samples, a fit fails,
  /// or the common core contains fewer than two samples.
  [[nodiscard]]
  GaussianCoreResult gaussianCore(double nSigma = 2.0, std::size_t minEntries = 20) const;

  /// Convenience function.
  ///
  /// Returns the Gaussian-core covariance when possible,
  /// otherwise falls back to the ordinary full covariance.
  [[nodiscard]]
  Matrix coreCovariance(double nSigma = 2.0, std::size_t minEntries = 20) const
  {
    const auto core = gaussianCore(nSigma, minEntries);

    if (core.valid) {
      return core.covariance;
    }

    return covariance();
  }

  [[nodiscard]]
  const Vector& lastResidual() const
  {
    return mLastResidual;
  }

  void computePulls(std::vector<Vector>& pulls) const
  {
    pulls.clear();
    pulls.reserve(mSamples.size());

    const auto core = gaussianCore();

    if (!core.valid) {
      return;
    }

    for (const auto& sample : mSamples) {
      Vector pull = Vector::Zero();
      for (int i = 0; i < kSize; ++i) {
        pull(i) = (sample(i) - core.mean(i)) / core.sigma(i);
      }
      pulls.push_back(pull);
    }
  }

  void computePulls(TH2F* pull, float pt, int i) const
  {
    std::vector<Vector> pullVectors;
    computePulls(pullVectors);
    for (const auto& p : pullVectors) {
      pull->Fill(pt, p(i));
    }
  }

 private:
  std::size_t mEntries = 0;
  std::size_t mRejected = 0;

  Vector mMean = Vector::Zero();
  Matrix mM2 = Matrix::Zero();

  Vector mLastResidual = Vector::Zero();

  std::vector<Vector> mSamples;
};

} // namespace AliceActsTrk