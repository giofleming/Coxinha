// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "FlatLutEntry.h"

#include <array>
#include <cstddef>
#include <string>
#include <vector>

namespace AliceActsTrk
{

struct LutRootRow {
  int entryIndex = -1;
  int nchBin = -1;
  int radBin = -1;
  int etaBin = -1;
  int ptBin = -1;

  float nch = 0.f;
  float eta = 0.f;
  float pt = 0.f;

  bool valid = false;
  float eff = 0.f;
  float eff2 = 0.f;
  float itof = 0.f;
  float otof = 0.f;

  std::array<float, 15> covm{};
  std::array<float, 5> eigval{};
  std::array<float, 25> eigvec{};
  std::array<float, 25> eiginv{};
};

class LutRootData
{
 public:
  LutRootData(const o2::delphes::lutHeader_t& header, bool covarianceStoresSigma);

  void reserve(std::size_t n);
  void addEntry(int entryIndex,
                int nchBin,
                int radBin,
                int etaBin,
                int ptBin,
                const o2::delphes::lutEntry_t& entry);

  bool writeToFile(const std::string& rootFileName) const;

 private:
  o2::delphes::lutHeader_t mHeader;
  bool mCovarianceStoresSigma = false;
  std::vector<LutRootRow> mRows;
};

} // namespace AliceActsTrk
