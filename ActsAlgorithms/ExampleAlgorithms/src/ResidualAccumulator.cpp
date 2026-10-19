// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <TF1.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <cstdlib>
#include "ResidualAccumulator.hpp"

using namespace AliceActsTrk;

ResidualAccumulator::GaussianCoreResult ResidualAccumulator::gaussianCore(double nSigma, std::size_t minEntries) const
{
  GaussianCoreResult result;

  if (mSamples.size() < minEntries) {
    return result;
  }

  // ------------------------------------------------------------
  // First determine the Gaussian core independently for each
  // of the five residual parameters.
  // ------------------------------------------------------------

  for (int ipar = 0; ipar < kSize; ++ipar) {

    std::vector<double> values;
    values.reserve(mSamples.size());

    for (const auto& sample : mSamples) {
      values.push_back(sample(ipar));
    }

    std::sort(values.begin(), values.end());

    const auto quantile = [&values](double q) -> double {
      if (values.empty()) {
        return 0.;
      }

      const double pos = q * static_cast<double>(values.size() - 1);
      const auto low = static_cast<std::size_t>(std::floor(pos));
      const auto high = static_cast<std::size_t>(std::ceil(pos));

      if (low == high) {
        return values[low];
      }

      const double fraction = pos - static_cast<double>(low);
      return values[low] * (1. - fraction) + values[high] * fraction;
    };

    // Robust starting values.
    double mean = quantile(0.5);

    double sigma = 0.5 * (quantile(0.8413447460685429) - quantile(0.1586552539314571));

    if (!std::isfinite(sigma) || sigma <= 0.) {

      // Fall back to the ordinary standard deviation
      // only for initialization.
      double sum = 0.;
      double sum2 = 0.;

      for (const double value : values) {
        sum += value;
        sum2 += value * value;
      }

      const double n = static_cast<double>(values.size());

      mean = sum / n;

      sigma = std::sqrt(std::max(0., sum2 / n - mean * mean));
    }

    if (!std::isfinite(mean) || !std::isfinite(sigma) || sigma <= 0.) {
      return result;
    }

    // Do not allow extreme tails to define the histogram range.
    const double histMin = mean - 6. * sigma;

    const double histMax = mean + 6. * sigma;

    const int nBins = std::clamp(static_cast<int>(4. * std::sqrt(static_cast<double>(values.size()))), 50, 200);

    TH1D hist(Form("hResidualCore_%p_%d", static_cast<const void*>(this), ipar), "", nBins, histMin, histMax);
    hist.SetDirectory(nullptr);

    for (const double value : values) {
      hist.Fill(value);
    }

    bool fitSucceeded = false;

    // Iterative Gaussian-core fit.
    //
    // First iteration is slightly wider, then converge
    // using the requested nSigma range.
    for (int iteration = 0; iteration < 3; ++iteration) {

      const double fitNSigma = (iteration == 0) ? 2.5 : nSigma;
      const double fitMin = mean - fitNSigma * sigma;
      const double fitMax = mean + fitNSigma * sigma;
      TF1 gaussian(Form("fResidualCore_%p_%d_%d", static_cast<const void*>(this), ipar, iteration), "gaus", fitMin, fitMax);
      gaussian.SetParameters(hist.GetMaximum(), mean, sigma);
      const auto fitResult = hist.Fit(&gaussian, "QNR0S");
      if (static_cast<int>(fitResult) != 0) {

        if (!fitSucceeded) {
          return result;
        }

        // A later refinement failed:
        // retain the previous successful iteration.
        break;
      }

      const double fittedMean = gaussian.GetParameter(1);
      const double fittedSigma = std::abs(gaussian.GetParameter(2));
      if (!std::isfinite(fittedMean) || !std::isfinite(fittedSigma) || fittedSigma <= 0.) {
        if (!fitSucceeded) {
          return result;
        }
        break;
      }

      mean = fittedMean;
      sigma = fittedSigma;

      fitSucceeded = true;
    }

    if (!fitSucceeded) {
      return result;
    }

    result.mean(ipar) = mean;
    result.sigma(ipar) = sigma;
  }

  // ------------------------------------------------------------
  // Select ONE common multidimensional Gaussian core.
  //
  // This is preferable to selecting a different population for
  // every covariance element because the resulting correlation
  // matrix remains mathematically consistent.
  // ------------------------------------------------------------

  std::vector<Vector> coreSamples;
  coreSamples.reserve(mSamples.size());

  for (const auto& sample : mSamples) {

    bool insideCore = true;

    for (int ipar = 0; ipar < kSize; ++ipar) {
      const double sigma = result.sigma(ipar);
      if (sigma <= 0.) {
        insideCore = false;
        break;
      }
      const double pull = (sample(ipar) - result.mean(ipar)) / sigma;
      if (!std::isfinite(pull) || std::abs(pull) > nSigma) {
        insideCore = false;
        break;
      }
    }

    if (insideCore) {
      coreSamples.push_back(sample);
    }
  }

  result.coreEntries = coreSamples.size();

  if (coreSamples.size() < 2) {
    return result;
  }

  // ------------------------------------------------------------
  // Determine the correlation matrix of the common core.
  // ------------------------------------------------------------

  Vector coreMean = Vector::Zero();
  for (const auto& sample : coreSamples) {
    coreMean += sample;
  }
  coreMean /= static_cast<double>(coreSamples.size());
  Matrix coreCov = Matrix::Zero();
  for (const auto& sample : coreSamples) {
    const Vector delta = sample - coreMean;
    coreCov += delta * delta.transpose();
  }
  coreCov /= static_cast<double>(coreSamples.size() - 1);
  coreCov = 0.5 * (coreCov + coreCov.transpose());
  Matrix correlation = Matrix::Identity();
  for (int i = 0; i < kSize; ++i) {
    for (int j = 0; j < i; ++j) {
      const double varianceI = coreCov(i, i);
      const double varianceJ = coreCov(j, j);
      double rho = 0.;
      if (varianceI > 0. && varianceJ > 0.) {
        rho = coreCov(i, j) / std::sqrt(varianceI * varianceJ);
        rho = std::clamp(rho, -1., 1.);
      }
      correlation(i, j) = rho;
      correlation(j, i) = rho;
    }
  }

  // ------------------------------------------------------------
  // Build final LUT covariance.
  //
  // The diagonal comes DIRECTLY from the Gaussian fits:
  //
  //   C_ii = sigma_i^2
  //
  // while off-diagonal entries retain the correlations measured
  // in the common core:
  //
  //   C_ij = rho_ij sigma_i sigma_j.
  // ------------------------------------------------------------

  result.covariance = result.sigma.asDiagonal() * correlation * result.sigma.asDiagonal();
  result.covariance = 0.5 * (result.covariance + result.covariance.transpose());
  result.valid = result.covariance.allFinite();
  return result;
}