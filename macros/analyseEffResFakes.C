// analyseEffResFakes.C
// =============================================================================
// Tracking-QA analysis of ACTS output (merged-over-iterations trees).
//
// It produces, for a chosen particle species, four families of plots vs pT, eta
// and phi (binning/ranges configurable; phi wrapped to [0, 2pi)):
//   (A) fake-hit rate per track  : fraction of tracks with EXACTLY k / AT LEAST k
//                                  fake hits (k = 1..maxFake);
//   (B) per-layer plots          : wrong-hit rate and hole rate per detector
//                                  layer (merged-layer aware);
//   (C) tracking efficiency      : "reconstructable" (denominator = findable
//                                  particles, filled markers) and "physics"
//                                  (denominator = all charged generator
//                                  primaries, open markers), split by the
//                                  matching track's fake count (all / exactly-k
//                                  / at-most-k);
//   (D) resolution & bias        : pT (relative), DCAxy, DCAz residual width
//                                  (Gaussian sigma in +-3 RMS) and mean, vs
//                                  pT/eta/phi, with QA residual distributions.
//
// ------------------------------- KEY DEFINITIONS ------------------------------
//   * A track state is a "measurement" if stateType == 0 (ACTS
//     RootTrackStatesWriter::StateType: 0=measurement,1=outlier,2=hole,3=material).
//   * Every particle is identified by a 5-component ACTS barcode
//       (vertex_primary, vertex_secondary, particle, generation, sub_particle).
//     The per-measurement contributors (trackstates particle_ids_*) and the
//     track's majority particle (tracksummary majorityParticleId_*) use it.
//   * FAKE hit  : a measurement whose contributor list does NOT contain the
//                 track's majority particle (so a hit with the majority among
//                 several contributors is NOT fake).
//   * FAKE track: a track with >= 1 fake hit.
//   * "matched" : the track has a real majority particle (nMajorityHits != UINT_MAX,
//                 or trackClassification, see useTrackClassification).
//
// --------------------------------- INPUTS ------------------------------------
//   tracksummary_..._merged.root  (tree "tracksummary") : 1 entry per EVENT,
//       branches are vectors over tracks (majority barcode, truth pT/eta/phi,
//       nMeasurements, fit residuals res_eLOC0/1_fit, eQOP_fit, eTHETA_fit, ...).
//   trackstates_..._merged.root   (tree "trackstates")  : 1 entry per TRACK,
//       branches are vectors over states (stateType, volume_id, layer_id,
//       particle_ids_*).
//   particles_*.root              (tree "particles")    : truth-PDG lookup,
//       efficiency denominators (selected + generator), per-particle kinematics.
//   hits.root                     (tree "hits")         : per-hit (volume,layer)
//       + barcode (masked-excluded counts and the required-layers cut).
//   digi(...).json (optional)     : volume_id -> volume name (labels, masking,
//       merging, required layers).
// The trees are joined by (event_nr / event_id, track_nr, barcode).
//
// --------------------------------- PIPELINE ----------------------------------
//   0) optional: read truth PDG map (particles) and per-hit info (hits.root);
//   1) read tracksummary -> trackInfo[(event,track)] = {majority, pt/eta/phi,
//      nMeas, matched, residuals};
//   2) book all histograms;
//   3) loop trackstates (one entry per track): apply the selection, count fake
//      hits, fill fake-rate / per-layer / resolution histos, and record the
//      track in recoMap (for the efficiency numerator);
//   4) turn counts into rates; build per-layer histos; efficiency (4b),
//      physics efficiency (4c), resolution/bias fits (4d);
//   5) draw everything and write a PDF + a ROOT file.
//
// Usage (from the directory that contains the .root files):
//   root -l -b -q analyseEffResFakes.C                 // all defaults
//   root -l -b -q 'analyseEffResFakes.C("." , 0.0,5.0,25, -2.5,2.5,32)'   // set ranges
// All knobs are function arguments with defaults (see the signature below).
// =============================================================================

#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TString.h>
#include <TLatex.h>
#include <TH2.h>
#include <TF1.h>

#include <array>
#include <algorithm>
#include <cmath>
#include <limits>
#include <cstdint>
#include <cstdlib>
#include <map>
#include <set>
#include <utility>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

// One truth barcode, kept as its 5 components so we can compare component-wise.
struct Barcode {
  unsigned int vp{0}, vs{0}, p{0}, g{0}, sp{0};
  bool operator==(const Barcode& o) const {
    return vp == o.vp && vs == o.vs && p == o.p && g == o.g && sp == o.sp;
  }
};

// Key into the truth-PDG map: (event, vertex_primary, vertex_secondary,
// particle, generation, sub_particle).
typedef std::array<unsigned int, 6> PKey;

// Info we need about each track, taken from the track summary tree.
struct TrackInfo {
  Barcode majority;
  bool matched{false};
  float pt{0.f};
  float eta{0.f};
  float phi{0.f};
  unsigned int nMeas{0};
  int majPdg{0};  // PDG of the majority particle (0 = unknown)
  float vz{0.f};  // production-vertex z of the majority particle [mm]
  // fit residuals (NaN if no fitted params): DCA in micrometers, pT relative
  float resD0Um{0.f};
  float resZ0Um{0.f};
  float resPtRel{0.f};
  bool hasRes{false};
};

// One fake measurement, kept for the optional printout.
struct FakeHit {
  int vol;
  int lay;
  std::vector<int> pdgs;  // PDG of EVERY particle contributing to the measurement
};

// Helper: make a ratio histogram (num/den) with binomial errors.
static TH1D* makeRatio(const TH1D* num, const TH1D* den, const char* name,
                       const char* title) {
  TH1D* r = (TH1D*)num->Clone(name);
  r->SetTitle(title);
  r->Divide(num, den, 1.0, 1.0, "B");  // "B" = binomial errors
  r->SetDirectory(nullptr);
  return r;
}

// Read volume_id -> volume name from an ACTS digitization-config JSON.
// Each entry contains a "NAME" followed by a "volume" id, so we simply pair
// every "volume" with the most recently seen "NAME". This minimal text scan
// avoids needing a JSON library inside the macro.
static bool loadVolumeNames(const std::string& path,
                            std::map<int, std::string>& names) {
  std::ifstream in(path.c_str());
  if (!in.good()) return false;
  std::stringstream buf;
  buf << in.rdbuf();
  const std::string s = buf.str();

  std::string lastName;
  std::size_t pos = 0;
  while (pos < s.size()) {
    std::size_t kn = s.find("\"NAME\"", pos);
    std::size_t kv = s.find("\"volume\"", pos);
    if (kn == std::string::npos && kv == std::string::npos) break;

    if (kn != std::string::npos && (kv == std::string::npos || kn < kv)) {
      // read the NAME string value
      std::size_t colon = s.find(':', kn);
      std::size_t q1 = s.find('"', colon + 1);
      std::size_t q2 = s.find('"', q1 + 1);
      if (q1 == std::string::npos || q2 == std::string::npos) break;
      lastName = s.substr(q1 + 1, q2 - q1 - 1);
      pos = q2 + 1;
    } else {
      // read the volume integer value
      std::size_t colon = s.find(':', kv);
      int vol = std::atoi(s.c_str() + colon + 1);
      if (!lastName.empty()) names[vol] = lastName;
      pos = colon + 1;
    }
  }
  return !names.empty();
}

// Trim leading/trailing whitespace.
static std::string trimStr(const std::string& s) {
  std::size_t a = s.find_first_not_of(" \t\n\r");
  if (a == std::string::npos) return "";
  std::size_t b = s.find_last_not_of(" \t\n\r");
  return s.substr(a, b - a + 1);
}

// Split a string on 'sep', trimming tokens and dropping empty ones.
static std::vector<std::string> splitStr(const std::string& s, char sep) {
  std::vector<std::string> out;
  std::stringstream ss(s);
  std::string tok;
  while (std::getline(ss, tok, sep)) {
    std::string t = trimStr(tok);
    if (!t.empty()) out.push_back(t);
  }
  return out;
}

// Pretty per-layer x-axis labels, keyed by (volume name, RAW layer_id), in an
// ORDERED list. The list order is also the x-axis bin order of the per-layer
// plots (bins not in this list are appended afterwards in natural order).
// layer_id is the per-layer-plot bin's (merged) representative raw id, i.e. the
// "keep" layer of a merge triplet (merge 8,10 -> bin id 8). Keying on the raw id
// makes the mapping robust to masking/merging: an id that is not present simply
// falls back to the generic label, so a mismatch never produces a *wrong* name.
// NOTE: the raw ids are derived from the default mask/merge config -- the
// per-layer terminal table prints the raw (volume, layer) of every bin, so any
// label that fell back to the generic form is easy to spot and fix.
struct PrettyLayerDef {
  const char* vol;
  int layer;
  const char* label;
};
static const std::vector<PrettyLayerDef> kPrettyLayerOrder = {
    // barrel
    {"InnerPixels::Barrel", 2, "VD0"},
    {"InnerPixels::Barrel", 4, "VD1"},
    {"InnerPixels::Barrel", 6, "VD2"},
    {"OuterPixels::Barrel", 2, "ML0"},
    {"OuterPixels::Barrel", 4, "ML1"},
    {"OuterPixels::Barrel", 6, "ML2"},
    {"OuterPixels::Barrel", 8, "ML3"},    // merged 8,10 (double-sided)
    {"OuterPixels::Barrel", 12, "iTOF"},
    {"OuterPixels::Barrel", 14, "ML4"},   // merged 14,16 (double-sided); iTOF (12) masked
    // OuterPixels::Barrel iTOF/ML5: layer 12 is masked by default; add an entry
    // here if you keep it, e.g. {"OuterPixels::Barrel", 12, "iTOF"}.
    {"OuterMostBarrel::Barrel", 2, "OT0"},
    {"OuterMostBarrel::Barrel", 4, "OT1"},
    {"OuterMostBarrel::Barrel", 6, "OT2"},
    // forward disks (each endcap: merged 2,4 -> 2 ; 6,8 -> 6 ; 10,12 -> 10)
    {"OuterPixels::PositiveEndcap", 2, "FD0"},
    {"OuterPixels::PositiveEndcap", 6, "FD1"},
    {"OuterPixels::PositiveEndcap", 10, "FD2"},
    {"OuterMostEndCap::PositiveEndcap", 2, "FD3"},
    {"OuterMostEndCap::PositiveEndcap", 6, "FD4"},
    {"OuterMostEndCap::PositiveEndcap", 10, "FD5"},
    // backward disks
    {"OuterPixels::NegativeEndcap", 2, "BD0"},
    {"OuterPixels::NegativeEndcap", 6, "BD1"},
    {"OuterPixels::NegativeEndcap", 10, "BD2"},
    {"OuterMostEndCap::NegativeEndcap", 2, "BD3"},
    {"OuterMostEndCap::NegativeEndcap", 6, "BD4"},
    {"OuterMostEndCap::NegativeEndcap", 10, "BD5"},

    {"oTOF::Barrel", 2, "oTOF"},
};

void analyseEffResFakes(const char* dirRec = ".",
                 const char* dirSim = ".",  // folder with the sim .root files (particles*, hits)
                 double ptMin = 0.05, double ptMax = 2.0, int ptNbins = 25,
                 double etaMin = -2.5, double etaMax = 2.5, int etaNbins = 32,
                 // phi axis + acceptance. phi is wrapped to [0, 2pi).
                 const char* pidList = "999",  // comma/space separated PDG codes, e.g. "211,321,2212,11,13"; keep a track if |majority PDG| matches ANY (antiparticles always included). "999" or empty = all species.
                 double phiMin = 0.0, double phiMax = 6.283185307179586,
                 int phiNbins = 32,
                 // Collision-vertex z acceptance [cm] applied to ALL plots: a
                 // track/particle is kept only if its production-vertex z (truth
                 // t_vz / vz) is in [z_vtx_min, z_vtx_max]. (converted to mm below)
                 double z_vtx_min = -10.0, double z_vtx_max = 10.0,
                 int maxFake = 3,
                 unsigned int nMeasMin = 7,
                 const char* digiConfig = "/Users/igor/ALICE3/actso2/geometries/geometry_2026_08_06_realistic_OT_with_eosCards_4percBarrel_4percDiscs_Gen3/digiConfigurations/digi-segm_config_pitch20um_with_barrel_TOFs-gen3.json",  // path to digi(...).json; "" = use volume numbers
                 bool logPt = true,            // logarithmic pT axis (log-spaced bins)
                 double ptMinLog = 0.05,       // lowest pT bin edge when logPt (must be > 0)
                 int minTracksPerLayer = 10,   // per-layer bins with fewer crossing tracks are set to 0
                 // Layers to EXCLUDE from the fake-hit calculation (everywhere).
                 //   maskById   : "vol,layer; vol,layer; ..."     e.g. "26,2; 25,4"
                 //   maskByName : "VolName,layer; ..." (needs digiConfig) e.g. "RICH::Barrel,2"
                 const char* maskById = "",
                 const char* maskByName = "", // "RICH::Barrel,2; RICH::Barrel,4", //; oTOF::Barrel,2; OuterPixels::Barrel,12",
                 // Layers to MERGE into one per-layer bin (merges layerB into layerA):
                 //   mergeTriplets       : "vol,layerA,layerB; ..."          e.g. "19,2,4"
                 //   mergeTripletsByName : "VolName,layerA,layerB; ..." (needs digiConfig)
                 //                                                   e.g. "OuterPixels::Barrel,2,4"
                 const char* mergeTriplets = "",
                 const char* mergeTripletsByName = "", // "OuterPixels::Barrel,8,10;OuterPixels::Barrel,14,16; OuterMostEndCap::NegativeEndcap,2,4;OuterMostEndCap::NegativeEndcap,6,8;OuterMostEndCap::NegativeEndcap,10,12; OuterMostEndCap::PositiveEndcap,2,4;OuterMostEndCap::PositiveEndcap,6,8;OuterMostEndCap::PositiveEndcap,10,12;   OuterPixels::NegativeEndcap,2,4;OuterPixels::NegativeEndcap,6,8;OuterPixels::NegativeEndcap,10,12; OuterPixels::PositiveEndcap,2,4;OuterPixels::PositiveEndcap,6,8;OuterPixels::PositiveEndcap,10,12",
                 bool printFakes = true,    // print details of tracks with >=1 fake hit
                 int nFakesToPrint = 10,    // how many such tracks to print
                 const char* particlesFile = "particles_digitized_selected.root",
                 const char* particlesTree = "particles",
                 // Efficiency: denominator = selectable truth particles.
                 const char* particlesEffFile = "particles_digitized_selected.root",
                 const char* particlesEffTree = "particles",
                 int nHitsMin = 5, //-1,  // denominator threshold; -1 => use nMeasMin. Interpreted as min #distinct layers when usePerLayerCountsForMeasAndHits.
                 const char* summaryFile = "tracksummary_ambi-tracks-merged.root",
                 const char* statesFile = "trackstates_ambi-tracks-merged.root",
                 const char* summaryTree = "tracksummary",
                 const char* statesTree = "trackstates",
                 const char* outFile = "EffResFakes_out.root",
                 bool printPerLayer = true,      // print the per-layer cross/wrong/rate table
                 bool dropLowStatLayers = true,  // drop bins with cross < minTracksPerLayer
                 // If true, nMeasMin and nHitsMin are interpreted as a MINIMUM
                 // NUMBER OF DISTINCT LAYERS carrying a measurement (reco/numerator)
                 // or a sim hit (truth/denominator), instead of the number of
                 // measurements / hits. Masked layers (maskById/maskByName) are
                 // always excluded and layer merging (mergeTriplets*) is applied in
                 // both modes. Requires hitsFile for the denominator layer count.
                 bool usePerLayerCountsForMeasAndHits = true,
                 bool strictCountCheck = true,   // mismatch vs stored count: true=abort, false=warn
                 const char* hitsFile = "hits.root",
                 const char* hitsTree = "hits",
                 // Require the (truth) particle to have sim hits on ALL listed layers,
                 // else the track/particle is excluded from ALL plots. Empty = no requirement.
                 //   requiredById   : "vol,layer; vol,layer; ..."
                 //   requiredByName : "VolName,layer; ..." (needs digiConfig)
                 // Matched against raw sim-hit layers (no merging). Uses hits.root.
                 const char* requiredById = "",
                 const char* requiredByName = "",
                 // Keep only primary particles (generation==0 && vertex_secondary==0),
                 // applied to the efficiency denominator and to each track's majority particle.
                 bool primaryOnly = true,
                 // How to decide a track is "matched" to a truth particle:
                 //   false (default): nMajorityHits != UINT_MAX (= Matched or Duplicate);
                 //   true           : read trackClassification == Matched (and Duplicate if
                 //                    includeDuplicates).
                 bool useTrackClassification = false,
                 bool includeDuplicates = true,
                 // "Physics" efficiency denominator: all charged primaries from the
                 // generator file (already-primary by construction). Drawn as open
                 // markers on the efficiency canvases. Empty file => skipped.
                 const char* particlesGenFile = "particles.root",
                 const char* particlesGenTree = "particles",
                 // ---- resolution block (pT, DCAxy, DCAz vs pT/eta/phi) ----
                 // Residual axis half-ranges: pT relative; DCAxy/DCAz in micrometers.
                 double resPtRelMax = 0.5,
                 double resDcaXyUmMax = 250.0,
                 double resDcaZUmMax = 800.0,
                 int resNbins = 400,
                 bool drawResiduals = false,  // QA residual-distribution canvases on/off
                 // QA (pT, eta) bins, "min:max,min:max,..." (':' separator allows negative eta)
                 const char* qaPtBins = "0.1:0.2,0.2:0.3,0.5:0.6,0.9:1.1,1.9:2.1,3.9:4.1",
                 const char* qaEtaBins = "-0.1:0.1,0.1:0.3,0.9:1.1,1.4:1.6,1.9:2.1,2.3:2.5",
                 bool drawBias = false,    // draw the bias (mean residual) canvases
                 // Use the friendly per-layer x-axis labels (VD0, ML0, BD0, FD0, ...)
                 // from kPrettyLayer (keyed on volume name + RAW layer id) instead of
                 // "<volume name> L<counter>". Needs digiConfig. Unmatched ids fall back.
                 bool prettyLayerNames = true,
                 // ---- fake-hit ORIGIN breakdown ----------------------------------
                 // Split fake hits by the origin of their contributing particle:
                 //   primary  : contributor is a generator primary (generation==0 &&
                 //              vertex_secondary==0), from ANY collision (signal/pileup);
                 //   material : contributor is a simulation secondary produced by a
                 //              material interaction (photon conversion / bremsstrahlung
                 //              / nuclear interaction);
                 //   decay    : contributor is a simulation secondary produced by a decay.
                 // The production process is read from originParticlesFile (must contain
                 // secondaries WITH their "process" branch -- the full simulated set,
                 // NOT the primaries-only digitized selection). Secondaries absent from
                 // that file (e.g. very low-hit delta rays) fall into neither series.
                 bool drawFakeOrigins = true,
                 const char* originParticlesFile = "particles_simulation.root",
                 // Per-layer mean number of sim particles contributing to each
                 // on-track hit (cluster truth multiplicity), averaged over the
                 // measurements of the selected tracks. Useful e.g. to see how
                 // occupancy/merging depends on the pixel size of each layer.
                 bool drawContribMult = true,
                 // Draw the injected-noise origin series on the fake-origin
                 // plots (per-layer + vs pt/eta/phi). Set false to hide it.
                 bool drawNoise = true) {
  gStyle->SetOptStat(0);

  // Collision-vertex z acceptance, converted from cm (params) to mm (the unit of
  // the t_vz / vz branches). Applied to every plot via the track + denominator
  // selections below.
  const double zVtxMinMm = z_vtx_min * 10.0;
  const double zVtxMaxMm = z_vtx_max * 10.0;

  // phi is wrapped into [0, 2pi) so the default [0, 2pi] range covers everything.
  auto wrapPhi = [](double p) {
    const double twoPi = 6.283185307179586;
    while (p < 0) p += twoPi;
    while (p >= twoPi) p -= twoPi;
    return p;
  };

  // Optional: volume_id -> volume name, read from the digitization-config JSON.
  // digiConfig is a FULL PATH to the JSON (it lives in its own directory).
  // If no file is given (or it cannot be read) we fall back to volume numbers.
  std::map<int, std::string> volNames;
  if (digiConfig && std::string(digiConfig).size() > 0) {
    loadVolumeNames(digiConfig, volNames);
    if (!volNames.empty()) {
      std::cout << "Loaded " << volNames.size() << " volume names from "
                << digiConfig << std::endl;
    } else {
      std::cout << "WARNING: could not read volume names from '" << digiConfig
                << "'; using volume numbers." << std::endl;
    }
  }
  // reverse map name -> volume id, used for name-based masking
  std::map<std::string, int> name2vol;
  for (const auto& kv : volNames) name2vol[kv.second] = kv.first;

  // small helper to print "volume X (Name)" or just "volume X"
  auto volStr = [&](int v) {
    auto it = volNames.find(v);
    return it != volNames.end() ? Form("volume %d (%s)", v, it->second.c_str())
                                : Form("volume %d", v);
  };

  // ---- masked (volume, layer) pairs: excluded from the fake-hit calc -------
  std::set<std::pair<int, int>> masked;
  for (const auto& e : splitStr(maskById, ';')) {
    auto f = splitStr(e, ',');
    if (f.size() == 2) masked.insert({std::atoi(f[0].c_str()), std::atoi(f[1].c_str())});
  }
  for (const auto& e : splitStr(maskByName, ';')) {
    auto f = splitStr(e, ',');
    if (f.size() != 2) continue;
    auto it = name2vol.find(f[0]);
    if (it == name2vol.end()) {
      std::cout << "Warning: mask volume name '" << f[0]
                << "' not found in digi config; entry ignored." << std::endl;
      continue;
    }
    masked.insert({it->second, std::atoi(f[1].c_str())});
  }
  if (!masked.empty()) {
    std::cout << "Warning: the following (volume, layer) are EXCLUDED from the "
                 "fake-hit calculation:" << std::endl;
    for (const auto& m : masked)
      std::cout << "  " << volStr(m.first) << ", layer " << m.second << std::endl;
  }

  // ---- layer merging: child (vol,layer) -> parent (vol,layer) --------------
  std::map<std::pair<int, int>, std::pair<int, int>> mergeParent;
  auto addMerge = [&](int v, int a, int b) {
    mergeParent[{v, b}] = {v, a};  // fold layer b into layer a
    std::cout << "Warning: layers " << a << " and " << b << " (" << volStr(v)
              << ") are merged into layer " << a << std::endl;
  };
  for (const auto& e : splitStr(mergeTriplets, ';')) {
    auto f = splitStr(e, ',');
    if (f.size() == 3) {
      addMerge(std::atoi(f[0].c_str()), std::atoi(f[1].c_str()), std::atoi(f[2].c_str()));
    }
  }
  for (const auto& e : splitStr(mergeTripletsByName, ';')) {
    auto f = splitStr(e, ',');
    if (f.size() != 3) continue;
    auto it = name2vol.find(f[0]);
    if (it == name2vol.end()) {
      std::cout << "Warning: merge volume name '" << f[0]
                << "' not found in digi config; entry ignored." << std::endl;
      continue;
    }
    addMerge(it->second, std::atoi(f[1].c_str()), std::atoi(f[2].c_str()));
  }
  // follow the chain so e.g. 6->4->2 resolves to 2
  auto resolveLayer = [&](std::pair<int, int> key) {
    int guard = 0;
    while (mergeParent.count(key) && guard++ < 1000) key = mergeParent[key];
    return key;
  };

  // ---- required (volume, layer) pairs: a track's truth particle must have sim
  //      hits on ALL of these, otherwise it is excluded from every plot. ------
  std::set<std::pair<int, int>> requiredLayers;
  for (const auto& e : splitStr(requiredById, ';')) {
    auto f = splitStr(e, ',');
    if (f.size() == 2) requiredLayers.insert({std::atoi(f[0].c_str()), std::atoi(f[1].c_str())});
  }
  for (const auto& e : splitStr(requiredByName, ';')) {
    auto f = splitStr(e, ',');
    if (f.size() != 2) continue;
    auto it = name2vol.find(f[0]);
    if (it == name2vol.end()) {
      std::cout << "Warning: required volume name '" << f[0]
                << "' not found in digi config; entry ignored." << std::endl;
      continue;
    }
    requiredLayers.insert({it->second, std::atoi(f[1].c_str())});
  }
  if (!requiredLayers.empty()) {
    std::cout << "Warning: only particles with sim hits on ALL of these layers are kept:"
              << std::endl;
    for (const auto& r : requiredLayers)
      std::cout << "  " << volStr(r.first) << ", layer " << r.second << std::endl;
  }

  // ---- per-particle sim-hit info from hits.root (read once, up-front) --------
  // hitCounts: (event,barcode) -> (#all hits, #hits on non-masked layers)
  // hitLayers: (event,barcode) -> set of raw (volume,layer) the particle hit
  // Needed when masked counts are requested and/or a required-layer list is set.
  const bool needRequired = !requiredLayers.empty();
  // Per-particle layer sets (hitLayers) are needed for the required-layers cut and
  // for per-layer denominator counting; non-masked hit counts also need hits.root.
  const bool needLayerSets = needRequired || usePerLayerCountsForMeasAndHits;
  const bool needHits = needLayerSets || !masked.empty();
  // Announce how the nMeasMin / nHitsMin thresholds are interpreted this run.
  {
    const int effHitsMinMsg = (nHitsMin < 0) ? (int)nMeasMin : nHitsMin;
    if (usePerLayerCountsForMeasAndHits)
      std::cout << "cut mode: per-layer -- require >= " << nMeasMin
                << " distinct layers with a measurement (numerator) and >= "
                << effHitsMinMsg << " with a sim hit (denominator); "
                << "masked layers excluded, merging applied." << std::endl;
    else
      std::cout << "cut mode: measurements/hits -- require >= " << nMeasMin
                << " measurements (numerator) and >= " << effHitsMinMsg
                << " sim hits (denominator); masked layers excluded." << std::endl;
  }
  std::map<PKey, std::pair<long, long>> hitCounts;
  std::map<PKey, std::set<std::pair<int, int>>> hitLayers;
  if (needHits) {
    std::string hitsPath = std::string(dirSim) + "/" + hitsFile;
    TFile* fHit = TFile::Open(hitsPath.c_str(), "READ");
    if (fHit == nullptr || fHit->IsZombie()) {
      std::cout << "ERROR: cannot open " << hitsPath
                << " (needed for masking / required / per-layer counts). Aborting."
                << std::endl;
      return;
    }
    TTreeReader rHit(hitsTree, fHit);
    TTreeReaderValue<std::uint32_t> hEv(rHit, "event_id");
    TTreeReaderValue<std::uint32_t> hVol(rHit, "volume_id");
    TTreeReaderValue<std::uint32_t> hLay(rHit, "layer_id");
    TTreeReaderValue<std::uint32_t> hVP(rHit, "barcode_vertex_primary");
    TTreeReaderValue<std::uint32_t> hVS(rHit, "barcode_vertex_secondary");
    TTreeReaderValue<std::uint32_t> hPa(rHit, "barcode_particle");
    TTreeReaderValue<std::uint32_t> hGe(rHit, "barcode_generation");
    TTreeReaderValue<std::uint32_t> hSu(rHit, "barcode_sub_particle");
    while (rHit.Next()) {
      PKey key{*hEv, *hVP, *hVS, *hPa, *hGe, *hSu};
      std::pair<int, int> lay{(int)*hVol, (int)*hLay};
      auto& c = hitCounts[key];
      c.first += 1;                            // all layers
      if (!masked.count(lay)) c.second += 1;   // non-masked
      if (needLayerSets) hitLayers[key].insert(lay);  // raw layer (no merging)
    }
    std::cout << "Loaded sim hits for " << hitCounts.size() << " particles from "
              << hitsFile << std::endl;
  }
  // Does the particle (event,barcode) have sim hits on ALL required layers?
  auto passesRequired = [&](const PKey& key) -> bool {
    if (requiredLayers.empty()) return true;
    auto it = hitLayers.find(key);
    if (it == hitLayers.end()) return false;
    for (const auto& req : requiredLayers) {
      if (!it->second.count(req)) return false;
    }
    return true;
  };

  // ------------------------------------------------------------------
  // 0) Optional: truth PDG per particle, keyed by (event, barcode).
  //    Needed for the pid filter and for the fake-track printout.
  // ------------------------------------------------------------------
  // Parse the PID selection list: a comma/space/semicolon separated list of
  // integer PDG codes (e.g. "211,321,2212,11,13"). Antiparticles are always
  // included (the comparison uses |PDG|). "999" or an empty list => all species.
  std::vector<int> pidCodes;   // codes as given, for the file-name tag
  std::set<int> pidAbsSet;     // |PDG| values to keep
  {
    std::string tok;
    auto flush = [&]() {
      if (!tok.empty()) {
        const int v = std::atoi(tok.c_str());
        pidCodes.push_back(v);
        pidAbsSet.insert(std::abs(v));
        tok.clear();
      }
    };
    for (const char* p = (pidList ? pidList : ""); *p; ++p) {
      if (*p == ',' || *p == ';' || *p == ' ' || *p == '\t')
        flush();
      else
        tok.push_back(*p);
    }
    flush();
  }
  const bool allSpecies = pidCodes.empty() || pidAbsSet.count(999) > 0;
  // Keep a track/particle whose species is selected (all-species => always true).
  auto passesPid = [&](int pdg) -> bool {
    return allSpecies || pidAbsSet.count(std::abs(pdg)) > 0;
  };
  if (allSpecies) {
    std::cout << "species filter: all species" << std::endl;
  } else {
    std::cout << "species filter: |PDG| in {";
    bool firstPid = true;
    for (int v : pidAbsSet) {
      std::cout << (firstPid ? "" : ",") << v;
      firstPid = false;
    }
    std::cout << "} (antiparticles included)" << std::endl;
  }

  const bool needPdg = !allSpecies || printFakes;
  std::map<PKey, int> pdgMap;
  if (needPdg) {
    std::string parPath = std::string(dirRec) + "/" + particlesFile;
    TFile* fPar = TFile::Open(parPath.c_str(), "READ");
    if (fPar == nullptr || fPar->IsZombie()) {
      std::cout << "WARNING: cannot open " << parPath
                << " -> truth PDG not available." << std::endl;
    } else {
      TTreeReader rPar(particlesTree, fPar);
      TTreeReaderValue<std::uint32_t> pEv(rPar, "event_id");
      TTreeReaderValue<std::vector<std::int32_t>> pType(rPar, "particle_type");
      TTreeReaderValue<std::vector<std::uint32_t>> pVP(rPar, "vertex_primary");
      TTreeReaderValue<std::vector<std::uint32_t>> pVS(rPar, "vertex_secondary");
      TTreeReaderValue<std::vector<std::uint32_t>> pPa(rPar, "particle");
      TTreeReaderValue<std::vector<std::uint32_t>> pGe(rPar, "generation");
      TTreeReaderValue<std::vector<std::uint32_t>> pSu(rPar, "sub_particle");
      while (rPar.Next()) {
        for (std::size_t i = 0; i < pType->size(); ++i) {
          pdgMap[PKey{*pEv, (*pVP)[i], (*pVS)[i], (*pPa)[i], (*pGe)[i], (*pSu)[i]}] =
              (*pType)[i];
        }
      }
      std::cout << "Loaded PDG for " << pdgMap.size() << " particles from "
                << particlesFile << std::endl;
    }
    if (!allSpecies && pdgMap.empty()) {
      std::cout << "ERROR: pid filter requested (pidList=\"" << pidList
                << "\") but no PDG information is available." << std::endl;
      return;
    }
  }
  // (event, barcode) -> PDG; returns 0 if the particle is not in the map.
  auto pdgOf = [&](unsigned int ev, const Barcode& b) -> int {
    auto it = pdgMap.find(PKey{ev, b.vp, b.vs, b.p, b.g, b.sp});
    return it == pdgMap.end() ? 0 : it->second;
  };

  // ------------------------------------------------------------------
  // 0b) Optional: production process per particle, keyed by (event, barcode).
  //     Used to label a fake hit's secondary contributor as decay vs material.
  //     ActsFatras::GenerationProcess codes: 0=undefined, 1=decay,
  //     2=photon conversion, 3=bremsstrahlung, 4=nuclear interaction.
  //     Read from originParticlesFile (the FULL simulated set, which contains
  //     secondaries with their "process" branch).
  std::map<PKey, int> procMap;
  if (drawFakeOrigins) {
    std::string oPath = std::string(dirSim) + "/" + originParticlesFile;
    TFile* fOri = TFile::Open(oPath.c_str(), "READ");
    if (fOri == nullptr || fOri->IsZombie()) {
      std::cout << "WARNING: cannot open " << oPath
                << " -> decay/material split of fake hits not available." << std::endl;
    } else {
      TTreeReader rOri(particlesTree, fOri);
      TTreeReaderValue<std::uint32_t> oEv(rOri, "event_id");
      TTreeReaderValue<std::vector<std::uint32_t>> oVP(rOri, "vertex_primary");
      TTreeReaderValue<std::vector<std::uint32_t>> oVS(rOri, "vertex_secondary");
      TTreeReaderValue<std::vector<std::uint32_t>> oPa(rOri, "particle");
      TTreeReaderValue<std::vector<std::uint32_t>> oGe(rOri, "generation");
      TTreeReaderValue<std::vector<std::uint32_t>> oSu(rOri, "sub_particle");
      TTreeReaderValue<std::vector<std::uint32_t>> oPr(rOri, "process");
      while (rOri.Next()) {
        for (std::size_t i = 0; i < oPr->size(); ++i) {
          procMap[PKey{*oEv, (*oVP)[i], (*oVS)[i], (*oPa)[i], (*oGe)[i], (*oSu)[i]}] =
              (int)(*oPr)[i];
        }
      }
      std::cout << "Loaded production process for " << procMap.size()
                << " particles from " << originParticlesFile << std::endl;
    }
  }
  // (event, barcode) -> process code; -1 if the particle is not in the map.
  auto procOf = [&](unsigned int ev, const Barcode& b) -> int {
    auto it = procMap.find(PKey{ev, b.vp, b.vs, b.p, b.g, b.sp});
    return it == procMap.end() ? -1 : it->second;
  };

  // Classify a FAKE measurement by its contributor list into an origin series:
  //   0 = primary (>=1 contributor is a generator primary, any collision),
  //   1 = material (a secondary from conversion/brems/nuclear interaction),
  //   2 = decay    (a secondary from a decay),
  //   3 = noise    (NO contributor at all = injected/unmapped noise: a
  //                 pure-noise cluster from PerPixelNoiseAlgorithm has no truth
  //                 link, so its contributor list is empty),
  //   4 = unknown  (a secondary with undefined/unavailable process).
  //   Priority: primary > decay > material > noise > unknown.
  // When the origin is material (return value 1), matSub is set to the material
  // sub-process of the first material contributor (sub = process_code - 2):
  // 0 = photon conversion, 1 = bremsstrahlung, 2 = nuclear interaction,
  // 3 = delta-ray (ionisation), 4 = other. Otherwise matSub is -1.
  auto classifyOrigin = [&](unsigned int ev,
                            const std::vector<std::uint32_t>& vp,
                            const std::vector<std::uint32_t>& vs,
                            const std::vector<std::uint32_t>& pp,
                            const std::vector<std::uint32_t>& gg,
                            const std::vector<std::uint32_t>& ss, int& matSub) -> int {
    matSub = -1;
    bool anyDecay = false;
    int firstMat = -1;  // material sub-process of the first material contributor
    for (std::size_t c = 0; c < vp.size(); ++c) {
      Barcode b{vp[c], vs[c], pp[c], gg[c], ss[c]};
      if (b.g == 0 && b.vs == 0) return 0;  // primary contributor wins
      const int proc = procOf(ev, b);       // GenerationProcess of the secondary
      if (proc == 1) anyDecay = true;
      // proc 2..6 = conversion/brems/nuclear/delta-ray/other -> sub 0..4 (proc-2).
      // Bremsstrahlung (proc 3, sub 1) is folded into "other" (sub 4).
      else if (proc >= 2 && proc <= 6 && firstMat < 0) {
        int sub = proc - 2;
        if (sub == 1) sub = 4;  // brems -> other
        firstMat = sub;
      }
    }
    if (anyDecay) return 2;
    if (firstMat >= 0) { matSub = firstMat; return 1; }
    if (vp.empty()) return 3;  // no contributor at all = injected/unmapped noise
    return 4;                  // has contributor(s) but undefined process
  };

  // ------------------------------------------------------------------
  // 1) Read the track summary and remember, per (event, track), the
  //    majority particle and the truth pT / eta.
  // ------------------------------------------------------------------
  std::string summaryPath = std::string(dirRec) + "/" + summaryFile;
  TFile* fSum = TFile::Open(summaryPath.c_str(), "READ");
  if (fSum == nullptr || fSum->IsZombie()) {
    std::cerr << "ERROR: cannot open " << summaryPath << std::endl;
    return;
  }

  TTreeReader rSum(summaryTree, fSum);
  TTreeReaderValue<std::uint32_t> sEvent(rSum, "event_nr");
  TTreeReaderValue<std::vector<std::uint32_t>> sTrackNr(rSum, "track_nr");
  TTreeReaderValue<std::vector<unsigned int>> sNMeas(rSum, "nMeasurements");
  TTreeReaderValue<std::vector<unsigned int>> sNMaj(rSum, "nMajorityHits");
  TTreeReaderValue<std::vector<int>> sCls(rSum, "trackClassification");
  TTreeReaderValue<std::vector<std::uint32_t>> sMVP(rSum, "majorityParticleId_vertex_primary");
  TTreeReaderValue<std::vector<std::uint32_t>> sMVS(rSum, "majorityParticleId_vertex_secondary");
  TTreeReaderValue<std::vector<std::uint32_t>> sMP(rSum, "majorityParticleId_particle");
  TTreeReaderValue<std::vector<std::uint32_t>> sMG(rSum, "majorityParticleId_generation");
  TTreeReaderValue<std::vector<std::uint32_t>> sMSP(rSum, "majorityParticleId_sub_particle");
  TTreeReaderValue<std::vector<float>> sPt(rSum, "t_pT");
  TTreeReaderValue<std::vector<float>> sEta(rSum, "t_eta");
  TTreeReaderValue<std::vector<float>> sPhi(rSum, "t_phi");
  TTreeReaderValue<std::vector<float>> sVz(rSum, "t_vz");  // majority-particle prod. z [mm]
  TTreeReaderValue<std::vector<float>> sResD0(rSum, "res_eLOC0_fit");  // d0 residual [mm]
  TTreeReaderValue<std::vector<float>> sResZ0(rSum, "res_eLOC1_fit");  // z0 residual [mm]
  TTreeReaderValue<std::vector<float>> sQOPfit(rSum, "eQOP_fit");
  TTreeReaderValue<std::vector<float>> sThetaFit(rSum, "eTHETA_fit");

  std::map<std::pair<unsigned int, unsigned int>, TrackInfo> trackInfo;
  // Set of reconstructed (processed) event ids. The summary tree has one row per
  // processed event -- even events that yielded zero tracks -- so this is exactly
  // the set of reco events, used to restrict the efficiency denominators below.
  std::set<std::uint32_t> recoEvents;

  while (rSum.Next()) {
    recoEvents.insert(*sEvent);
    const unsigned int nTracks = sTrackNr->size();
    for (unsigned int i = 0; i < nTracks; ++i) {
      TrackInfo info;
      info.majority = Barcode{(*sMVP)[i], (*sMVS)[i], (*sMP)[i], (*sMG)[i], (*sMSP)[i]};
      // "matched" definition. ACTS marks unmatched/fake tracks with
      // nMajorityHits == UINT_MAX and an all-zero majority barcode.
      // TrackMatchClassification: 0=Unknown, 1=Matched, 2=Duplicate, 3=Fake.
      if (useTrackClassification) {
        const int cls = (*sCls)[i];
        info.matched = (cls == 1) || (includeDuplicates && cls == 2);
      } else {
        info.matched = (*sNMaj)[i] != std::numeric_limits<unsigned int>::max();
      }
      info.pt = (*sPt)[i];
      info.eta = (*sEta)[i];
      info.phi = wrapPhi((*sPhi)[i]);
      info.vz = (*sVz)[i];
      info.nMeas = (*sNMeas)[i];
      info.majPdg = needPdg ? pdgOf(*sEvent, info.majority) : 0;
      // fit residuals (NaN unless matched + fitted): DCA mm->um, pT relative
      info.resD0Um = (*sResD0)[i] * 1000.f;
      info.resZ0Um = (*sResZ0)[i] * 1000.f;
      const float qop = (*sQOPfit)[i];
      const float thf = (*sThetaFit)[i];
      if (std::isfinite(qop) && qop != 0.f && std::isfinite(thf) && info.pt > 0.f) {
        const float ptFit = std::sin(thf) / std::fabs(qop);
        info.resPtRel = (ptFit - info.pt) / info.pt;
      } else {
        info.resPtRel = std::numeric_limits<float>::quiet_NaN();
      }
      info.hasRes = std::isfinite(info.resD0Um) && std::isfinite(info.resZ0Um) &&
                    std::isfinite(info.resPtRel);
      trackInfo[{*sEvent, (*sTrackNr)[i]}] = info;
    }
  }
  std::cout << "Loaded " << trackInfo.size() << " tracks from the summary ("
            << recoEvents.size() << " reconstructed events)." << std::endl;

  // ------------------------------------------------------------------
  // 2) Book histograms.
  // ------------------------------------------------------------------
  // pT bin edges: equidistant, or log-spaced when logPt is requested.
  std::vector<double> ptEdges;
  {
    const double lo = (logPt && ptMin <= 0.0) ? ptMinLog : ptMin;
    if (logPt) {
      const double lmin = std::log10(lo), lmax = std::log10(ptMax);
      for (int i = 0; i <= ptNbins; ++i)
        ptEdges.push_back(std::pow(10.0, lmin + (lmax - lmin) * i / ptNbins));
    } else {
      for (int i = 0; i <= ptNbins; ++i)
        ptEdges.push_back(ptMin + (ptMax - ptMin) * i / ptNbins);
    }
  }
  auto makePtHist = [&](const char* name) {
    TH1D* h = new TH1D(name, ";p_{T} [GeV];tracks",
                       (int)ptEdges.size() - 1, ptEdges.data());
    h->SetDirectory(nullptr);
    return h;
  };
  // helper for the (linear) phi histograms, reused everywhere
  auto makePhiHist = [&](const char* name) {
    TH1D* h = new TH1D(name, ";#phi;entries", phiNbins, phiMin, phiMax);
    h->SetDirectory(nullptr);
    return h;
  };

  // Denominators: all SELECTED tracks.
  TH1D* hAllPt = makePtHist("hAllPt");
  TH1D* hAllEta = new TH1D("hAllEta", ";#eta;tracks", etaNbins, etaMin, etaMax);
  hAllEta->SetDirectory(nullptr);
  TH1D* hAllPhi = makePhiHist("hAllPhi");

  // ---- QA count distributions (filled for every track passing the kinematic +
  //      matching selection, BEFORE the nMeas/nLayers cut): number of
  //      measurements and number of distinct layers-with-measurement for the
  //      reco track; number of sim hits and distinct layers-with-hit for its
  //      majority (truth) particle. All masked-excluded and merge-aware, i.e.
  //      consistent with the actual cut quantities. ---------------------------
  auto makeCountHist = [&](const char* name, const char* xt) {
    TH1D* h = new TH1D(name, Form(";%s;tracks", xt), 41, -0.5, 40.5);
    h->SetDirectory(nullptr);
    return h;
  };
  TH1D* hQaNMeas    = makeCountHist("qa_nMeasurements",    "n measurements (reco track)");
  TH1D* hQaNLayMeas = makeCountHist("qa_nLayersWithMeas",  "n layers with measurement (reco track)");
  TH1D* hQaNHits    = makeCountHist("qa_nHits",            "n hits (majority particle)");
  TH1D* hQaNLayHits = makeCountHist("qa_nLayersWithHits",  "n layers with hit (majority particle)");
  // Overlap/merge multiplicity: extra measurements beyond one-per-layer (>=0).
  TH1D* hQaMeasMinusLay =
      makeCountHist("qa_nMeasMinusLayers", "n measurements - n layers (reco track)");

  // Numerators for "exactly k" and "at least k" fake hits (k = 1..maxFake).
  std::vector<TH1D*> hExactPt(maxFake + 1, nullptr), hExactEta(maxFake + 1, nullptr);
  std::vector<TH1D*> hAtLeastPt(maxFake + 1, nullptr), hAtLeastEta(maxFake + 1, nullptr);
  std::vector<TH1D*> hExactPhi(maxFake + 1, nullptr), hAtLeastPhi(maxFake + 1, nullptr);
  for (int k = 1; k <= maxFake; ++k) {
    hExactPt[k] = makePtHist(Form("hExactPt_%d", k));
    hExactEta[k] = new TH1D(Form("hExactEta_%d", k), ";#eta;tracks", etaNbins, etaMin, etaMax);
    hAtLeastPt[k] = makePtHist(Form("hAtLeastPt_%d", k));
    hAtLeastEta[k] = new TH1D(Form("hAtLeastEta_%d", k), ";#eta;tracks", etaNbins, etaMin, etaMax);
    hExactEta[k]->SetDirectory(nullptr);
    hAtLeastEta[k]->SetDirectory(nullptr);
    hExactPhi[k] = makePhiHist(Form("hExactPhi_%d", k));
    hAtLeastPhi[k] = makePhiHist(Form("hAtLeastPhi_%d", k));
  }

  // Fake-hit ORIGIN numerators: per-track "has >=1 fake hit of origin o", with
  // o = 0 primary, 1 material, 2 decay. Denominator = all selected tracks.
  // Plotted fake-origin series: primary, material, decay, noise (injected).
  const char* kOriginName[4] = {"primary", "material", "decay", "noise"};
  // Material sub-processes (the breakdown of the "material" origin), in the order
  // 0 = photon conversion, 1 = bremsstrahlung, 2 = nuclear interaction. Long
  // labels for printout/legends, short tags for histogram object names.
  // NB: bremsstrahlung is folded into "other" (see classifyOrigin), so the
  // "bremsstrahlung" bin stays empty and is not drawn / printed.
  const char* kMatSubName[5] = {"photon-conversion", "bremsstrahlung",
                                "nuclear-interaction", "delta-ray", "other"};
  const char* kMatSubTag[5] = {"conv", "brems", "nuclear", "delta", "other"};
  const int kNMatSub = 5;
  std::vector<TH1D*> hOrigPt(4, nullptr), hOrigEta(4, nullptr), hOrigPhi(4, nullptr);
  if (drawFakeOrigins) {
    for (int o = 0; o < 4; ++o) {
      hOrigPt[o] = makePtHist(Form("hOrigPt_%s", kOriginName[o]));
      hOrigEta[o] = new TH1D(Form("hOrigEta_%s", kOriginName[o]), ";#eta;tracks",
                             etaNbins, etaMin, etaMax);
      hOrigEta[o]->SetDirectory(nullptr);
      hOrigPhi[o] = makePhiHist(Form("hOrigPhi_%s", kOriginName[o]));
    }
  }

  // ---- resolution: residual TH2 (binning axis vs residual) ------------------
  // parameters: 0 = pT (relative), 1 = DCAxy [um], 2 = DCAz [um]
  // axes:       0 = pT, 1 = eta, 2 = phi
  // categories: 0 = all matched; 1+j = "<= j fake" (j = 0..maxFake)
  const char* resParName[3] = {"ptrel", "dcaxy", "dcaz"};
  const double resParMax[3] = {resPtRelMax, resDcaXyUmMax, resDcaZUmMax};
  const char* resSigTitle[3] = {"#sigma(p_{T})/p_{T}", "#sigma(DCA_{xy}) [#mum]",
                                "#sigma(DCA_{z}) [#mum]"};
  const char* resBiasTitle[3] = {"#LT#Deltap_{T}/p_{T}#GT", "#LTDCA_{xy}#GT [#mum]",
                                 "#LTDCA_{z}#GT [#mum]"};
  const char* resResTitle[3] = {"#Deltap_{T}/p_{T}", "DCA_{xy} residual [#mum]",
                                "DCA_{z} residual [#mum]"};
  const char* resAxName[3] = {"pt", "eta", "phi"};
  const char* resAxTitle[3] = {"p_{T} [GeV]", "#eta", "#phi"};
  const int nResCat = maxFake + 2;  // 0 = all matched; 1+j = "<= j fake"

  auto makeRes2 = [&](int p, int a, const char* name) -> TH2D* {
    const double ym = resParMax[p];
    TH2D* h = nullptr;
    if (a == 0) {
      h = new TH2D(name, "", (int)ptEdges.size() - 1, ptEdges.data(), resNbins, -ym, ym);
    } else if (a == 1) {
      h = new TH2D(name, "", etaNbins, etaMin, etaMax, resNbins, -ym, ym);
    } else {
      h = new TH2D(name, "", phiNbins, phiMin, phiMax, resNbins, -ym, ym);
    }
    h->SetDirectory(nullptr);
    return h;
  };
  std::vector<std::vector<std::vector<TH2D*>>> res2(
      3, std::vector<std::vector<TH2D*>>(3, std::vector<TH2D*>(nResCat, nullptr)));
  for (int p = 0; p < 3; ++p)
    for (int a = 0; a < 3; ++a)
      for (int c = 0; c < nResCat; ++c)
        res2[p][a][c] = makeRes2(p, a, Form("res2_%s_%s_c%d", resParName[p], resAxName[a], c));

  // QA residual histograms per (param, eta-bin, pt-bin); bins from "min:max,..."
  std::vector<std::pair<double, double>> qaPt, qaEta;
  for (const auto& e : splitStr(qaPtBins, ',')) {
    auto f = splitStr(e, ':');
    if (f.size() == 2) qaPt.push_back({std::atof(f[0].c_str()), std::atof(f[1].c_str())});
  }
  for (const auto& e : splitStr(qaEtaBins, ',')) {
    auto f = splitStr(e, ':');
    if (f.size() == 2) qaEta.push_back({std::atof(f[0].c_str()), std::atof(f[1].c_str())});
  }
  const int nQaPt = (int)qaPt.size();
  const int nQaEta = (int)qaEta.size();
  std::vector<std::vector<std::vector<TH1D*>>> qaRes(
      3, std::vector<std::vector<TH1D*>>(nQaEta, std::vector<TH1D*>(nQaPt, nullptr)));
  if (drawResiduals) {
    for (int p = 0; p < 3; ++p)
      for (int ie = 0; ie < nQaEta; ++ie)
        for (int ip = 0; ip < nQaPt; ++ip) {
          TH1D* h = new TH1D(Form("qa_%s_e%d_p%d", resParName[p], ie, ip), "", resNbins,
                             -resParMax[p], resParMax[p]);
          h->SetDirectory(nullptr);
          qaRes[p][ie][ip] = h;
        }
  }

  // Per-layer counters, keyed by the MERGED (volume_id, layer_id):
  //   layerStats: (#tracks with a wrong hit, #tracks crossing with a measurement)
  //   holeStats : (#tracks with a hole,      #tracks expected to cross the layer)
  // "expected to cross" = the track has a sensitive state (measurement/outlier/hole)
  // on that layer. Hole rate is then defined as holes / tracks expected to cross.
  std::map<std::pair<int, int>, std::pair<long, long>> layerStats;
  std::map<std::pair<int, int>, std::pair<long, long>> holeStats;
  // Per-layer wrong-hit counts split by origin (same crossing denominator as
  // layerStats): index 0 = primary, 1 = material, 2 = decay.
  std::vector<std::map<std::pair<int, int>, std::pair<long, long>>> layerStatsOrig(4);
  // Per-layer wrong-hit counts for the material sub-processes (0 = photon
  // conversion, 1 = bremsstrahlung, 2 = nuclear interaction). Same crossing
  // denominator; the numerators are subsets of the material numerator above.
  std::vector<std::map<std::pair<int, int>, std::pair<long, long>>> layerStatsMatSub(kNMatSub);
  // Per (merged) layer: accumulator for the mean contributor multiplicity of
  // on-track hits. Stored as {n hits, sum of multiplicities, sum of squares}
  // so we can form the mean and its error (RMS/sqrt(n)).
  std::map<std::pair<int, int>, std::array<double, 3>> contribStats;

  // For efficiency: per truth particle matched by a selected track, remember the
  // fake-hit count of its best (most-measurements) matching track.
  struct MatchInfo {
    int bestNMeas = -1;
    int nFake = 0;
  };
  std::map<PKey, MatchInfo> recoMap;

  // ------------------------------------------------------------------
  // 3) Loop over track states (one entry per track) and fill.
  // ------------------------------------------------------------------
  std::string statesPath = std::string(dirRec) + "/" + statesFile;
  TFile* fSta = TFile::Open(statesPath.c_str(), "READ");
  if (fSta == nullptr || fSta->IsZombie()) {
    std::cerr << "ERROR: cannot open " << statesPath << std::endl;
    return;
  }

  TTreeReader rSta(statesTree, fSta);
  TTreeReaderValue<std::uint32_t> tEvent(rSta, "event_nr");
  TTreeReaderValue<std::uint32_t> tTrackNr(rSta, "track_nr");
  TTreeReaderValue<std::vector<int>> tStateType(rSta, "stateType");
  TTreeReaderValue<std::vector<int>> tVol(rSta, "volume_id");
  TTreeReaderValue<std::vector<int>> tLay(rSta, "layer_id");
  TTreeReaderValue<std::vector<std::vector<std::uint32_t>>> tVP(rSta, "particle_ids_vertex_primary");
  TTreeReaderValue<std::vector<std::vector<std::uint32_t>>> tVS(rSta, "particle_ids_vertex_secondary");
  TTreeReaderValue<std::vector<std::vector<std::uint32_t>>> tP(rSta, "particle_ids_particle");
  TTreeReaderValue<std::vector<std::vector<std::uint32_t>>> tG(rSta, "particle_ids_generation");
  TTreeReaderValue<std::vector<std::vector<std::uint32_t>>> tSP(rSta, "particle_ids_sub_particle");

  long nSelected = 0, nFakeTracks = 0;
  // Among fake tracks, how many have >=1 fake hit of each origin (a track with
  // mixed-origin fakes is counted in several): 0 primary, 1 material, 2 decay,
  // 3 unknown. Used for the terminal fraction printout.
  long nFakeOrig[5] = {0, 0, 0, 0, 0};
  // Among fake tracks, how many have >=1 material fake hit of each sub-process
  // (0 conversion, 1 brems, 2 nuclear). Subsets of nFakeOrig[1] (material).
  long nFakeMatSub[5] = {0, 0, 0, 0, 0};
  long nMajZero = 0;  // selected tracks whose majority particle has PDG == 0 (unknown)
  long nPrinted = 0;
  if (printFakes) {
    std::cout << "--- Fake tracks (printing up to " << nFakesToPrint
              << "; PID=0 means the measurement has no associated truth particle, "
                 "e.g. a noise / unmapped hit) ---"
              << std::endl;
  }

  // One iteration = one reconstructed track. We look up its per-track summary
  // info (majority particle, truth kinematics, residuals) via (event, track),
  // apply the common track selection, walk the track's states to count fake
  // hits and per-layer / hole info, and then fill all the track-based histograms
  // and record the track in recoMap (used later as the efficiency numerator).
  while (rSta.Next()) {
    auto it = trackInfo.find({*tEvent, *tTrackNr});
    if (it == trackInfo.end()) {
      continue;  // no summary info for this track (should not happen)
    }
    const TrackInfo& info = it->second;

    // ---- common track selection (applied to ALL track-based plots) ----------
    if (!info.matched) continue;
    if (primaryOnly && (info.majority.g != 0 || info.majority.vs != 0))
      continue;  // primaries only (generation==0 && vertex_secondary==0)
    if (!passesPid(info.majPdg)) continue;  // species filter (any of pidList)
    if (info.pt < ptMin || info.pt > ptMax) continue;
    if (info.eta < etaMin || info.eta > etaMax) continue;
    if (info.phi < phiMin || info.phi > phiMax) continue;
    if (info.vz < zVtxMinMm || info.vz > zVtxMaxMm) continue;  // collision-vertex z
    // Required-layers cut: the majority particle must have sim hits on all of them.
    if (!passesRequired(PKey{*tEvent, info.majority.vp, info.majority.vs,
                             info.majority.p, info.majority.g, info.majority.sp}))
      continue;
    // nMeasMin is applied AFTER the state loop, so it can optionally use the
    // masked-excluded measurement count.

    // Walk the states of this track.
    int nFake = 0;
    int nMeasAll = 0;   // all measurement states (should equal summary nMeasurements)
    int nMeasUsed = 0;  // measurement states on non-masked layers
    // For the per-layer plot, one bool per (merged) layer this track touches:
    // "does the track have >=1 wrong hit on that layer?". Several measurements on
    // the same layer (overlapping sensors, or layers combined via merging) are
    // OR-ed into this single bool -> the track is counted ONCE per layer
    // (crossing) and flagged wrong if ANY of those measurements is fake.
    std::map<std::pair<int, int>, bool> layerWrong;
    // For the per-layer HOLE plot, per (merged) layer this track touches:
    //   layerExpected = the track has a sensitive state (measurement/outlier/hole) there;
    //   layerHasHole  = the track has >=1 hole there.
    std::set<std::pair<int, int>> layerExpected;
    std::set<std::pair<int, int>> layerHasHole;
    // For the optional printout: each fake measurement on this track, with the
    // PDGs of all particles that contributed to it.
    std::vector<FakeHit> fakeHits;
    // Fake-hit origin bookkeeping for this track: per-track "has >=1 fake hit of
    // origin o" (o = 0 primary, 1 material, 2 decay, 3 unknown), and the set of
    // (merged) layers on which the track has >=1 wrong hit of that origin (only
    // the 3 plotted origins; unknown is tracked for the terminal printout only).
    bool fakeFromOrigin[5] = {false, false, false, false, false};
    std::vector<std::set<std::pair<int, int>>> originWrong(4);
    // Same, restricted to material fakes split by sub-process (0 conversion,
    // 1 brems, 2 nuclear): per-track flag + the layers carrying such a fake hit.
    bool fakeFromMatSub[5] = {false, false, false, false, false};
    std::vector<std::set<std::pair<int, int>>> originMatSub(5);

    const unsigned int nStates = tStateType->size();
    for (unsigned int s = 0; s < nStates; ++s) {
      const int st = (*tStateType)[s];  // 0=measurement, 1=outlier, 2=hole, 3=material
      if (st == 0) ++nMeasAll;          // all measurement states (summary cross-check)

      // Only sensitive-surface states (measurement/outlier/hole) define a layer
      // crossing; material/unknown states are ignored everywhere.
      if (st != 0 && st != 1 && st != 2) continue;

      // Masked layers are excluded from every per-layer calculation.
      std::pair<int, int> rawLayer{(*tVol)[s], (*tLay)[s]};
      if (masked.count(rawLayer)) continue;
      const std::pair<int, int> mlayer = resolveLayer(rawLayer);  // merged key

      // Per-layer HOLE bookkeeping: any sensitive state => expected to cross.
      layerExpected.insert(mlayer);
      if (st == 2) layerHasHole.insert(mlayer);

      if (st != 0) continue;  // the rest (fakes / wrong-hit plot) is measurements only
      ++nMeasUsed;

      // Is the majority particle among the contributors of this measurement?
      bool correct = false;
      const auto& vp = (*tVP)[s];
      const auto& vs = (*tVS)[s];
      const auto& pp = (*tP)[s];
      const auto& gg = (*tG)[s];
      const auto& ss = (*tSP)[s];
      // Per-layer mean contributor multiplicity: number of sim particles
      // associated to THIS hit (= contributor list length), accumulated per
      // (merged) layer over all on-track measurements.
      if (drawContribMult) {
        const double m = (double)vp.size();
        auto& acc = contribStats[mlayer];
        acc[0] += 1.0;      // hits
        acc[1] += m;        // sum of multiplicities
        acc[2] += m * m;    // sum of squares
      }
      for (unsigned int c = 0; c < vp.size(); ++c) {
        Barcode b{vp[c], vs[c], pp[c], gg[c], ss[c]};
        if (b == info.majority) {
          correct = true;
          break;
        }
      }
      const bool wrong = !correct;  // empty contributors also count as wrong
      if (wrong) {
        ++nFake;
        // Origin of this fake hit (0 primary, 1 material, 2 decay, 3 unknown).
        if (drawFakeOrigins) {
          int matSub = -1;
          const int o = classifyOrigin(*tEvent, vp, vs, pp, gg, ss, matSub);
          if (o >= 0 && o < 5) fakeFromOrigin[o] = true;
          if (o >= 0 && o < 4) originWrong[o].insert(mlayer);  // 4 plotted origins
          if (o == 1 && matSub >= 0 && matSub < 5) {           // material sub-process
            fakeFromMatSub[matSub] = true;
            originMatSub[matSub].insert(mlayer);
          }
        }
        if (printFakes) {
          std::vector<int> pdgs;
          for (unsigned int c = 0; c < vp.size(); ++c) {
            pdgs.push_back(pdgOf(*tEvent, Barcode{vp[c], vs[c], pp[c], gg[c], ss[c]}));
          }
          fakeHits.push_back({rawLayer.first, rawLayer.second, std::move(pdgs)});
        }
      }

      // Group by the MERGED layer key (mlayer, computed above): merged layers
      // share one entry, so the OR spans all constituent layers of the group.
      auto lit = layerWrong.find(mlayer);
      if (lit == layerWrong.end()) {
        layerWrong[mlayer] = wrong;
      } else {
        lit->second = lit->second || wrong;
      }
    }

    // Consistency check: all measurement states must match the summary count.
    if (nMeasAll != (int)info.nMeas) {
      std::cout << (strictCountCheck ? "ERROR" : "Warning")
                << ": event " << *tEvent << " track " << *tTrackNr
                << " nMeasurements mismatch: trackstates=" << nMeasAll
                << " summary=" << info.nMeas << std::endl;
      if (strictCountCheck) {
        std::cout << "Aborting (strictCountCheck=true)." << std::endl;
        return;
      }
    }
    // ---- QA distributions (all selected tracks, BEFORE the nMeas/nLayers cut) --
    hQaNMeas->Fill(nMeasUsed);                    // measurements (masked-excluded)
    hQaNLayMeas->Fill((int)layerWrong.size());    // distinct merged non-masked layers
    hQaMeasMinusLay->Fill(nMeasUsed - (int)layerWrong.size());  // overlap/merge extra
    {
      const PKey mkey{*tEvent, info.majority.vp, info.majority.vs,
                      info.majority.p, info.majority.g, info.majority.sp};
      auto hcQA = hitCounts.find(mkey);
      if (hcQA != hitCounts.end())
        hQaNHits->Fill((int)hcQA->second.second);  // non-masked sim hits
      auto hlQA = hitLayers.find(mkey);
      if (hlQA != hitLayers.end()) {
        std::set<std::pair<int, int>> lays;
        for (const auto& raw : hlQA->second)
          if (!masked.count(raw)) lays.insert(resolveLayer(raw));
        hQaNLayHits->Fill((int)lays.size());       // distinct merged non-masked layers
      }
    }

    // nMeasMin cut. With usePerLayerCountsForMeasAndHits, count the distinct
    // (merged, non-masked) layers carrying a measurement -- exactly the set of
    // layerWrong keys built above; otherwise use the (masked-excluded) measurement
    // count. Masked layers are excluded in both cases (they are skipped earlier).
    const int nMeasForCut = usePerLayerCountsForMeasAndHits
                                ? (int)layerWrong.size()
                                : nMeasUsed;
    if (nMeasForCut < (int)nMeasMin) continue;

    // Fill denominators (all selected tracks).
    ++nSelected;
    if (info.majPdg == 0) ++nMajZero;
    hAllPt->Fill(info.pt);
    hAllEta->Fill(info.eta);
    hAllPhi->Fill(info.phi);

    // Remember this track's majority particle as reconstructed, keeping the best
    // (most-measurements) matching track to represent it in the efficiency.
    {
      PKey mk{*tEvent, info.majority.vp, info.majority.vs, info.majority.p,
              info.majority.g, info.majority.sp};
      MatchInfo& mm = recoMap[mk];
      if (nMeasForCut > mm.bestNMeas) {
        mm.bestNMeas = nMeasForCut;
        mm.nFake = nFake;
      }
    }

    if (nFake >= 1) ++nFakeTracks;

    // "exactly k"
    if (nFake >= 1 && nFake <= maxFake) {
      hExactPt[nFake]->Fill(info.pt);
      hExactEta[nFake]->Fill(info.eta);
      hExactPhi[nFake]->Fill(info.phi);
    }
    // "at least k"
    for (int k = 1; k <= maxFake && k <= nFake; ++k) {
      hAtLeastPt[k]->Fill(info.pt);
      hAtLeastEta[k]->Fill(info.eta);
      hAtLeastPhi[k]->Fill(info.phi);
    }
    // fake-origin: a track contributes to series o if it has >=1 fake hit of
    // that origin (a track with mixed-origin fakes can contribute to several).
    if (drawFakeOrigins) {
      for (int o = 0; o < 4; ++o)
        if (fakeFromOrigin[o]) {
          hOrigPt[o]->Fill(info.pt);
          hOrigEta[o]->Fill(info.eta);
          hOrigPhi[o]->Fill(info.phi);
        }
      for (int o = 0; o < 5; ++o)
        if (fakeFromOrigin[o]) ++nFakeOrig[o];  // for the terminal fraction printout
      for (int s = 0; s < kNMatSub; ++s)
        if (fakeFromMatSub[s]) ++nFakeMatSub[s];  // material sub-process breakdown
    }

    // Resolution: residual vs (pT/eta/phi) for all matched + "<=k" categories,
    // and the QA residual distributions per (eta-bin, pt-bin).
    if (info.hasRes) {
      const double axval[3] = {info.pt, info.eta, info.phi};
      const double resval[3] = {info.resPtRel, info.resD0Um, info.resZ0Um};
      for (int p = 0; p < 3; ++p) {
        for (int a = 0; a < 3; ++a) {
          res2[p][a][0]->Fill(axval[a], resval[p]);            // all matched
          for (int j = nFake; j <= maxFake; ++j)
            res2[p][a][1 + j]->Fill(axval[a], resval[p]);      // "<= j fake"
        }
      }
      if (drawResiduals) {
        int ie = -1, ip = -1;
        for (int b = 0; b < nQaEta; ++b)
          if (info.eta >= qaEta[b].first && info.eta < qaEta[b].second) { ie = b; break; }
        for (int b = 0; b < nQaPt; ++b)
          if (info.pt >= qaPt[b].first && info.pt < qaPt[b].second) { ip = b; break; }
        if (ie >= 0 && ip >= 0)
          for (int p = 0; p < 3; ++p) qaRes[p][ie][ip]->Fill(resval[p]);
      }
    }

    // Per-layer accumulation: one count per (merged) layer this track crossed,
    // regardless of how many measurements it had there. -> per-layer rate is
    // "fraction of crossing tracks with >=1 wrong hit on the layer".
    for (const auto& lw : layerWrong) {
      layerStats[lw.first].second += 1;             // crossing (denominator)
      if (lw.second) layerStats[lw.first].first += 1;  // wrong (numerator)
      // Same crossing denominator for the per-origin breakdown; numerator counts
      // tracks with >=1 wrong hit of that origin on the layer.
      if (drawFakeOrigins) {
        for (int o = 0; o < 4; ++o) {
          layerStatsOrig[o][lw.first].second += 1;
          if (originWrong[o].count(lw.first)) layerStatsOrig[o][lw.first].first += 1;
        }
        for (int s = 0; s < kNMatSub; ++s) {  // material sub-processes (same denominator)
          layerStatsMatSub[s][lw.first].second += 1;
          if (originMatSub[s].count(lw.first)) layerStatsMatSub[s][lw.first].first += 1;
        }
      }
    }

    // Per-layer hole accumulation: one count per (merged) layer the track was
    // expected to cross; numerator counts those where it had a hole.
    for (const auto& lay : layerExpected) {
      holeStats[lay].second += 1;                          // expected (denominator)
      if (layerHasHole.count(lay)) holeStats[lay].first += 1;  // holes (numerator)
    }

    // Optional printout of fake tracks.
    if (printFakes && nFake >= 1 && nPrinted < nFakesToPrint) {
      ++nPrinted;
      std::cout << "[fake " << nPrinted << "] event=" << *tEvent
                << " track=" << *tTrackNr << " nMeas=" << nMeasForCut
                << " majorityPID=" << info.majPdg
                << " pt_truth=" << info.pt << " eta_truth=" << info.eta
                << " nFakeHits=" << nFake << std::endl;
      for (const auto& fh : fakeHits) {
        auto vn = volNames.find(fh.vol);
        TString vstr =
            (vn != volNames.end()) ? vn->second.c_str() : Form("%d", fh.vol);
        std::cout << "        fake hit: volume=" << vstr << " layer=" << fh.lay
                  << "  contributor PIDs: ";
        if (fh.pdgs.empty()) {
          std::cout << "(none)";
        } else {
          for (std::size_t i = 0; i < fh.pdgs.size(); ++i) {
            std::cout << (i ? ", " : "") << fh.pdgs[i];
          }
        }
        std::cout << std::endl;
      }
    }
  }

  std::cout << "Selected tracks: " << nSelected
            << ", fake tracks (>=1 fake hit): " << nFakeTracks;
  if (nSelected > 0) {
    std::cout << "  (" << 100.0 * nFakeTracks / nSelected << " %)";
  }
  std::cout << std::endl;

  // Fraction of fake tracks that have >=1 fake hit of each origin (the same
  // track can appear in several categories if its fakes have mixed origins).
  // Categories can overlap when a track's fakes have mixed origins, so the
  // percentages needn't sum to 100.
  if (drawFakeOrigins && nFakeTracks > 0) {
    const char* lbl[5] = {"primary", "decay", "material", "noise", "unknown"};
    const int idx[5] = {0, 2, 1, 3, 4};  // print order: primary, decay, material, noise, unknown
    std::cout << "Fake-track origin breakdown (fraction of "
              << nFakeTracks << " fake tracks with >=1 fake hit from ...; "
                 "categories can overlap, so they need not sum to 100%):" << std::endl;
    for (int i = 0; i < 5; ++i) {
      const long n = nFakeOrig[idx[i]];
      std::cout << "    " << lbl[i] << ": " << n << "  ("
                << 100.0 * n / nFakeTracks << " %)" << std::endl;
      // After the material line, break it down into its sub-processes (extra
      // indentation). These are subsets of "material" and can overlap.
      if (idx[i] == 1) {
        for (int s = 0; s < kNMatSub; ++s) {
          if (s == 1) continue;  // brems folded into "other"
          std::cout << "        " << kMatSubName[s] << ": " << nFakeMatSub[s] << "  ("
                    << 100.0 * nFakeMatSub[s] / nFakeTracks << " %)" << std::endl;
        }
      }
    }
  }

  std::cout << "Tracks with majorityPID==0 (truth PDG unavailable): " << nMajZero
            << " / " << nSelected;
  if (nSelected > 0) {
    std::cout << "  (" << 100.0 * nMajZero / nSelected << " %)";
  }
  std::cout << std::endl;
  if (!needPdg) {
    std::cout << "   note: truth PDG was not loaded (all species selected and printFakes==false),"
                 " so majorityPID is 0 for all tracks by construction." << std::endl;
  }

  // ------------------------------------------------------------------
  // 4) Turn the count histograms into rates (divide by all selected tracks).
  // ------------------------------------------------------------------
  std::vector<TH1D*> rExactPt(maxFake + 1, nullptr), rExactEta(maxFake + 1, nullptr);
  std::vector<TH1D*> rAtLeastPt(maxFake + 1, nullptr), rAtLeastEta(maxFake + 1, nullptr);
  std::vector<TH1D*> rExactPhi(maxFake + 1, nullptr), rAtLeastPhi(maxFake + 1, nullptr);
  for (int k = 1; k <= maxFake; ++k) {
    rExactPt[k] = makeRatio(hExactPt[k], hAllPt, Form("fakeRate_exact%d_vs_pT", k),
                            Form("exactly %d fake hits;p_{T} [GeV];fraction of tracks", k));
    rExactEta[k] = makeRatio(hExactEta[k], hAllEta, Form("fakeRate_exact%d_vs_eta", k),
                             Form("exactly %d fake hits;#eta;fraction of tracks", k));
    rExactPhi[k] = makeRatio(hExactPhi[k], hAllPhi, Form("fakeRate_exact%d_vs_phi", k),
                             Form("exactly %d fake hits;#phi;fraction of tracks", k));
    rAtLeastPt[k] = makeRatio(hAtLeastPt[k], hAllPt, Form("fakeRate_atLeast%d_vs_pT", k),
                              Form("#geq %d fake hits;p_{T} [GeV];fraction of tracks", k));
    rAtLeastEta[k] = makeRatio(hAtLeastEta[k], hAllEta, Form("fakeRate_atLeast%d_vs_eta", k),
                               Form("#geq %d fake hits;#eta;fraction of tracks", k));
    rAtLeastPhi[k] = makeRatio(hAtLeastPhi[k], hAllPhi, Form("fakeRate_atLeast%d_vs_phi", k),
                               Form("#geq %d fake hits;#phi;fraction of tracks", k));
  }

  // Fake-origin rates: fraction of selected tracks with >=1 fake hit of origin o.
  std::vector<TH1D*> rOrigPt(4, nullptr), rOrigEta(4, nullptr), rOrigPhi(4, nullptr);
  if (drawFakeOrigins) {
    for (int o = 0; o < 4; ++o) {
      rOrigPt[o] = makeRatio(hOrigPt[o], hAllPt, Form("fakeRate_origin_%s_vs_pT", kOriginName[o]),
                             Form("#geq 1 fake from %s;p_{T} [GeV];fraction of tracks", kOriginName[o]));
      rOrigEta[o] = makeRatio(hOrigEta[o], hAllEta, Form("fakeRate_origin_%s_vs_eta", kOriginName[o]),
                              Form("#geq 1 fake from %s;#eta;fraction of tracks", kOriginName[o]));
      rOrigPhi[o] = makeRatio(hOrigPhi[o], hAllPhi, Form("fakeRate_origin_%s_vs_phi", kOriginName[o]),
                              Form("#geq 1 fake from %s;#phi;fraction of tracks", kOriginName[o]));
    }
  }

  // Pretty-label lookup + x-axis display order, both from kPrettyLayerOrder
  // (single source of truth): prettyMap = (name,rawId)->label, prettyRank =
  // (name,rawId)->position (used to sort the per-layer bins).
  std::map<std::pair<std::string, int>, std::string> prettyMap;
  std::map<std::pair<std::string, int>, int> prettyRank;
  for (std::size_t i = 0; i < kPrettyLayerOrder.size(); ++i) {
    const auto& e = kPrettyLayerOrder[i];
    prettyMap[{e.vol, e.layer}] = e.label;
    prettyRank[{e.vol, e.layer}] = (int)i;
  }

  // Build a per-layer rate histogram (numerator/denominator per merged layer).
  // The "layer" here is AFTER merging (mergeTriplets / mergeTripletsByName). A
  // per-volume layer counter (0,1,2,...) is kept for the fallback label (raw
  // ACTS layer ids are even-only). The bins are ordered to follow
  // kPrettyLayerOrder; bins not listed there are appended in natural order.
  // Layers with denominator < minTracksPerLayer are dropped (dropLowStatLayers)
  // or set to 0. stats: (vol,layer) -> (num, den).
  auto buildPerLayerHist =
      [&](const std::map<std::pair<int, int>, std::pair<long, long>>& stats,
          const char* hname, const char* htitle, const char* headerName,
          const char* denWord, const char* numWord, bool printRows = true) -> TH1D* {
    struct Row { int vol, layerId, counter; long num, den; bool enough; };
    std::vector<Row> rows;
    int curVol = 0, layerCounter = 0;
    bool firstVol = true;
    for (const auto& ls : stats) {
      const int vol = ls.first.first;
      if (firstVol || vol != curVol) {
        curVol = vol;
        layerCounter = 0;
        firstVol = false;
      } else {
        ++layerCounter;
      }
      const long num = ls.second.first;
      const long den = ls.second.second;
      rows.push_back({vol, ls.first.second, layerCounter, num, den,
                      (den >= minTracksPerLayer && den > 0)});
    }

    // Order the bins to follow kPrettyLayerOrder (by volume name + raw layer id);
    // bins not listed there keep their natural (volume,layer) order at the end.
    auto rankOf = [&](const Row& r) -> int {
      auto vit = volNames.find(r.vol);
      if (vit == volNames.end()) return std::numeric_limits<int>::max();
      auto rit = prettyRank.find({vit->second, r.layerId});
      return rit != prettyRank.end() ? rit->second : std::numeric_limits<int>::max();
    };
    std::stable_sort(rows.begin(), rows.end(),
                     [&](const Row& a, const Row& b) { return rankOf(a) < rankOf(b); });

    int nShown = 0;
    for (const auto& r : rows) {
      if (!(dropLowStatLayers && !r.enough)) ++nShown;
    }
    const int nBins = (nShown > 0) ? nShown : 1;

    TH1D* h = new TH1D(hname, htitle, nBins, 0.5, nBins + 0.5);
    h->SetDirectory(nullptr);
    h->GetXaxis()->SetLabelSize(0.05);  // bigger per-layer (bin) labels
    h->GetYaxis()->SetLabelSize(0.04);

    if (printPerLayer && printRows) {
      std::cout << "--- " << headerName << " (minTracksPerLayer = "
                << minTracksPerLayer << ", dropLowStatLayers = "
                << (dropLowStatLayers ? "true" : "false") << ") ---" << std::endl;
    }
    int bin = 1;
    for (const auto& r : rows) {
      // x-axis bin label for this (volume, raw layer_id, per-volume counter):
      //  - if the volume name is known (digiConfig loaded) and prettyLayerNames is
      //    on AND (name, RAW layer_id) is in kPrettyLayer  -> friendly name
      //    (VD0, ML0, BD3, FD0, ...);
      //  - else if the name is known                       -> "<name> L<counter>";
      //  - else (no name)                                  -> "v<volId> L<counter>".
      // Pretty names are keyed on the raw layer id, so a missing/unmatched id
      // just falls back to the generic label (never a wrong pretty name).
      auto nit = volNames.find(r.vol);
      TString label;
      if (nit != volNames.end()) {
        TString pretty;
        if (prettyLayerNames) {
          auto pit = prettyMap.find({nit->second, r.layerId});
          if (pit != prettyMap.end()) pretty = pit->second.c_str();
        }
        label = pretty.Length() ? pretty
                                : TString(Form("%s L%d", nit->second.c_str(), r.counter));
      } else {
        label = Form("v%d L%d", r.vol, r.counter);
      }
      const bool drop = dropLowStatLayers && !r.enough;

      if (printPerLayer && printRows) {
        const double rate = (r.den > 0) ? (double)r.num / (double)r.den : 0.0;
        std::cout << "  " << label << "  (vol " << r.vol << ", layer " << r.layerId
                  << "): " << denWord << "=" << r.den << " " << numWord << "=" << r.num
                  << " rate=" << (r.enough ? rate : 0.0)
                  << (r.enough ? "" : (drop ? "  [removed: too few tracks]"
                                            : "  [zeroed: too few tracks]"))
                  << std::endl;
      }

      if (drop) continue;  // low-stat layer removed from the histogram

      double rate = 0.0, err = 0.0;
      if (r.enough) {
        rate = (double)r.num / (double)r.den;
        err = std::sqrt(rate * (1.0 - rate) / r.den);
      }
      h->SetBinContent(bin, rate);
      h->SetBinError(bin, err);
      h->GetXaxis()->SetBinLabel(bin, label);
      ++bin;
    }
    return h;
  };

  // Variant of buildPerLayerHist for a per-layer MEAN (not a ratio): stats holds
  // {n, sum, sumSq} per merged layer; the bin value is the mean = sum/n and the
  // error is the error on the mean = RMS/sqrt(n). Same bin ordering, pretty
  // labels and low-stat dropping (n < minTracksPerLayer) as buildPerLayerHist.
  auto buildPerLayerMean =
      [&](const std::map<std::pair<int, int>, std::array<double, 3>>& stats,
          const char* hname, const char* htitle, const char* headerName) -> TH1D* {
    struct Row { int vol, layerId, counter; double n, sum, sumSq; bool enough; };
    std::vector<Row> rows;
    int curVol = 0, layerCounter = 0;
    bool firstVol = true;
    for (const auto& ls : stats) {
      const int vol = ls.first.first;
      if (firstVol || vol != curVol) { curVol = vol; layerCounter = 0; firstVol = false; }
      else ++layerCounter;
      const double n = ls.second[0];
      rows.push_back({vol, ls.first.second, layerCounter, n, ls.second[1], ls.second[2],
                      (n >= minTracksPerLayer && n > 0)});
    }
    auto rankOf = [&](const Row& r) -> int {
      auto vit = volNames.find(r.vol);
      if (vit == volNames.end()) return std::numeric_limits<int>::max();
      auto rit = prettyRank.find({vit->second, r.layerId});
      return rit != prettyRank.end() ? rit->second : std::numeric_limits<int>::max();
    };
    std::stable_sort(rows.begin(), rows.end(),
                     [&](const Row& a, const Row& b) { return rankOf(a) < rankOf(b); });

    int nShown = 0;
    for (const auto& r : rows) if (!(dropLowStatLayers && !r.enough)) ++nShown;
    const int nBins = (nShown > 0) ? nShown : 1;
    TH1D* h = new TH1D(hname, htitle, nBins, 0.5, nBins + 0.5);
    h->SetDirectory(nullptr);
    h->GetXaxis()->SetLabelSize(0.05);  // bigger per-layer (bin) labels
    h->GetYaxis()->SetLabelSize(0.04);

    if (printPerLayer)
      std::cout << "--- " << headerName << " (minTracksPerLayer = " << minTracksPerLayer
                << " hits) ---" << std::endl;
    int bin = 1;
    for (const auto& r : rows) {
      auto nit = volNames.find(r.vol);
      TString label;
      if (nit != volNames.end()) {
        TString pretty;
        if (prettyLayerNames) {
          auto pit = prettyMap.find({nit->second, r.layerId});
          if (pit != prettyMap.end()) pretty = pit->second.c_str();
        }
        label = pretty.Length() ? pretty
                                : TString(Form("%s L%d", nit->second.c_str(), r.counter));
      } else {
        label = Form("v%d L%d", r.vol, r.counter);
      }
      const bool drop = dropLowStatLayers && !r.enough;
      const double mean = (r.n > 0) ? r.sum / r.n : 0.0;
      double err = 0.0;
      if (r.n > 1) {
        double var = (r.sumSq - r.n * mean * mean) / (r.n - 1.0);
        if (var < 0) var = 0;
        err = std::sqrt(var / r.n);  // error on the mean
      }
      if (printPerLayer)
        std::cout << "  " << label << "  (vol " << r.vol << ", layer " << r.layerId
                  << "): hits=" << (long)r.n << " <contrib>=" << (r.enough ? mean : 0.0)
                  << (r.enough ? "" : (drop ? "  [removed: too few hits]"
                                            : "  [zeroed: too few hits]"))
                  << std::endl;
      if (drop) continue;
      h->SetBinContent(bin, r.enough ? mean : 0.0);
      h->SetBinError(bin, r.enough ? err : 0.0);
      h->GetXaxis()->SetBinLabel(bin, label);
      ++bin;
    }
    return h;
  };

  // Per-layer wrong-hit rate = wrong / crossing tracks (see definition above).
  TH1D* hLayerRate = buildPerLayerHist(
      layerStats, "fakeRate_perLayer", "per-layer wrong-hit rate;;wrong / crossing tracks",
      "per-layer wrong-hit rate", "cross", "wrong");

  // Per-layer wrong-hit rate split by origin (primary / material / decay). Same
  // crossing denominator as hLayerRate, so the bins/labels line up one-to-one
  // (printRows=false: don't repeat the per-layer table for each origin).
  std::vector<TH1D*> hLayerRateOrig(4, nullptr);
  // Material sub-process per-layer histograms: hLayerMatSub = wrong-hit rate of
  // each sub-process vs crossing tracks (so it overlays the material total);
  // hLayerMatFrac = fraction of the layer's material fakes from that sub-process
  // (numerator = sub count, denominator = total material count on the layer).
  std::vector<TH1D*> hLayerMatSub(kNMatSub, nullptr), hLayerMatFrac(kNMatSub, nullptr);
  // Fraction of each layer's wrong hits coming from each origin (primary /
  // material / decay): numerator = origin count, denominator = total wrong
  // count on the layer. Origins can overlap on a (merged) layer and "unknown"
  // is not shown, so these need not sum to 1.
  std::vector<TH1D*> hLayerOrigFrac(4, nullptr);
  if (drawFakeOrigins) {
    for (int o = 0; o < 4; ++o) {
      hLayerRateOrig[o] = buildPerLayerHist(
          layerStatsOrig[o], Form("fakeRate_perLayer_%s", kOriginName[o]),
          ";;wrong / crossing tracks", "per-layer wrong-hit rate", "cross", "wrong", false);
      // Fraction of the layer's total wrong hits from this origin.
      std::map<std::pair<int, int>, std::pair<long, long>> fracStats;
      for (const auto& kv : layerStatsOrig[o]) {
        const auto& key = kv.first;
        const long num = kv.second.first;  // tracks with >=1 wrong hit of origin o
        const long tot = layerStats.count(key) ? layerStats.at(key).first : 0;  // all wrong
        fracStats[key] = {num, tot};
      }
      hLayerOrigFrac[o] = buildPerLayerHist(
          fracStats, Form("fakeFrac_perLayer_%s", kOriginName[o]),
          ";;fraction of wrong hits", "per-layer origin fraction", "wrong", "origin", false);
    }
    for (int s = 0; s < kNMatSub; ++s) {
      hLayerMatSub[s] = buildPerLayerHist(
          layerStatsMatSub[s], Form("fakeRate_perLayer_mat_%s", kMatSubTag[s]),
          ";;wrong / crossing tracks", "per-layer material wrong-hit rate", "cross", "wrong",
          false);
      // Build (sub count, total material count) per layer for the fraction plot.
      std::map<std::pair<int, int>, std::pair<long, long>> fracStats;
      for (const auto& kv : layerStatsMatSub[s]) {
        const auto& key = kv.first;
        const long sub = kv.second.first;
        const long matTot =
            layerStatsOrig[1].count(key) ? layerStatsOrig[1].at(key).first : 0;
        fracStats[key] = {sub, matTot};
      }
      hLayerMatFrac[s] = buildPerLayerHist(
          fracStats, Form("fakeFrac_perLayer_mat_%s", kMatSubTag[s]),
          ";;fraction of material fakes", "per-layer material fraction", "material", "sub",
          false);
    }
  }

  // Per-layer hole rate = holes / tracks expected to cross the layer. The "layer"
  // is after merging (mergeTriplets / mergeTripletsByName), just like above.
  TH1D* hLayerHoleRate = buildPerLayerHist(
      holeStats, "holeRate_perLayer", "per-layer hole rate;;holes / tracks expected to cross",
      "per-layer hole rate", "expected", "holes");

  // Per-layer mean number of sim particles contributing to each on-track hit.
  TH1D* hLayerContrib = nullptr;
  if (drawContribMult)
    hLayerContrib = buildPerLayerMean(
        contribStats, "contribMult_perLayer",
        "per-layer mean contributors per hit;;#LT sim particles / hit #GT",
        "per-layer mean contributors per hit");

  // ------------------------------------------------------------------
  // 4b) Tracking efficiency vs pT / eta, split by the fake-hit count of the
  //     matching track. Denominator = selectable truth particles; numerator =
  //     those matched to a selected track (optionally with k / >=k fake hits).
  // ------------------------------------------------------------------
  const int effHitsMin = (nHitsMin < 0) ? (int)nMeasMin : nHitsMin;
  bool haveEff = false;
  TH1D* hDenPt = nullptr;
  TH1D* hDenEta = nullptr;
  TH1D* hDenPhi = nullptr;
  TH1D* hNumAllPt = nullptr;
  TH1D* hNumAllEta = nullptr;
  TH1D* hNumAllPhi = nullptr;
  std::vector<TH1D*> hNumExactPt(maxFake + 1, nullptr), hNumExactEta(maxFake + 1, nullptr);
  std::vector<TH1D*> hNumAtMostPt(maxFake + 1, nullptr), hNumAtMostEta(maxFake + 1, nullptr);
  std::vector<TH1D*> hNumExactPhi(maxFake + 1, nullptr), hNumAtMostPhi(maxFake + 1, nullptr);
  TH1D* effAllPt = nullptr;
  TH1D* effAllEta = nullptr;
  TH1D* effAllPhi = nullptr;
  std::vector<TH1D*> effExactPt(maxFake + 1, nullptr), effExactEta(maxFake + 1, nullptr);
  std::vector<TH1D*> effAtMostPt(maxFake + 1, nullptr), effAtMostEta(maxFake + 1, nullptr);
  std::vector<TH1D*> effExactPhi(maxFake + 1, nullptr), effAtMostPhi(maxFake + 1, nullptr);
  {
    std::string effPath = std::string(dirRec) + "/" + particlesEffFile;
    TFile* fEff = TFile::Open(effPath.c_str(), "READ");
    if (fEff == nullptr || fEff->IsZombie()) {
      std::cout << "WARNING: cannot open " << effPath
                << " -> efficiency plots skipped." << std::endl;
    } else {
      auto makeEtaHist = [&](const char* name) {
        TH1D* h = new TH1D(name, ";#eta;particles", etaNbins, etaMin, etaMax);
        h->SetDirectory(nullptr);
        return h;
      };
      hDenPt = makePtHist("hDenPt");
      hDenEta = makeEtaHist("hDenEta");
      hDenPhi = makePhiHist("hDenPhi");
      hNumAllPt = makePtHist("hNumAllPt");
      hNumAllEta = makeEtaHist("hNumAllEta");
      hNumAllPhi = makePhiHist("hNumAllPhi");
      hNumExactPt[0] = makePtHist("hNumExactPt_0");    // clean tracks (0 fake hits)
      hNumExactEta[0] = makeEtaHist("hNumExactEta_0");
      hNumExactPhi[0] = makePhiHist("hNumExactPhi_0");
      hNumAtMostPt[0] = makePtHist("hNumAtMostPt_0");   // "<=0" = clean (same set)
      hNumAtMostEta[0] = makeEtaHist("hNumAtMostEta_0");
      hNumAtMostPhi[0] = makePhiHist("hNumAtMostPhi_0");
      for (int k = 1; k <= maxFake; ++k) {
        hNumExactPt[k] = makePtHist(Form("hNumExactPt_%d", k));
        hNumExactEta[k] = makeEtaHist(Form("hNumExactEta_%d", k));
        hNumExactPhi[k] = makePhiHist(Form("hNumExactPhi_%d", k));
        hNumAtMostPt[k] = makePtHist(Form("hNumAtMostPt_%d", k));
        hNumAtMostEta[k] = makeEtaHist(Form("hNumAtMostEta_%d", k));
        hNumAtMostPhi[k] = makePhiHist(Form("hNumAtMostPhi_%d", k));
      }

      TTreeReader rEff(particlesEffTree, fEff);
      TTreeReaderValue<std::uint32_t> pEv(rEff, "event_id");
      TTreeReaderValue<std::vector<std::int32_t>> pType(rEff, "particle_type");
      TTreeReaderValue<std::vector<float>> pPt(rEff, "pt");
      TTreeReaderValue<std::vector<float>> pEta(rEff, "eta");
      TTreeReaderValue<std::vector<float>> pPhi(rEff, "phi");
      TTreeReaderValue<std::vector<float>> pVz(rEff, "vz");  // production-vertex z [mm]
      TTreeReaderValue<std::vector<std::int32_t>> pNhits(rEff, "number_of_hits");
      TTreeReaderValue<std::vector<std::uint32_t>> pVP(rEff, "vertex_primary");
      TTreeReaderValue<std::vector<std::uint32_t>> pVS(rEff, "vertex_secondary");
      TTreeReaderValue<std::vector<std::uint32_t>> pPa(rEff, "particle");
      TTreeReaderValue<std::vector<std::uint32_t>> pGe(rEff, "generation");
      TTreeReaderValue<std::vector<std::uint32_t>> pSu(rEff, "sub_particle");

      // (hitCounts / hitLayers were filled up-front from hits.root, if needed.)

      long nDen = 0, nNum = 0;
      long nEffEvents = 0, nEffEvSkipped = 0;
      while (rEff.Next()) {
        ++nEffEvents;
        // Only count truth particles from events that were actually reconstructed;
        // otherwise the denominator is inflated when reco processed fewer events.
        if (!recoEvents.count(*pEv)) { ++nEffEvSkipped; continue; }
        for (std::size_t i = 0; i < pType->size(); ++i) {
          if (primaryOnly && ((*pGe)[i] != 0 || (*pVS)[i] != 0))
            continue;  // primaries only (generation==0 && vertex_secondary==0)
          if (!passesPid((*pType)[i])) continue;
          if ((*pVz)[i] < zVtxMinMm || (*pVz)[i] > zVtxMaxMm) continue;  // collision-vertex z
          const float dpt = (*pPt)[i];
          const float deta = (*pEta)[i];
          const float dphi = wrapPhi((*pPhi)[i]);
          if (dpt < ptMin || dpt > ptMax) continue;
          if (deta < etaMin || deta > etaMax) continue;
          if (dphi < phiMin || dphi > phiMax) continue;

          PKey key{*pEv, (*pVP)[i], (*pVS)[i], (*pPa)[i], (*pGe)[i], (*pSu)[i]};

          // Required-layers cut (same as for tracks): must have sim hits on all.
          if (!passesRequired(key)) continue;

          // Cut quantity for the denominator. With per-layer counting, use the
          // number of distinct (merged, non-masked) layers the particle has sim
          // hits on; otherwise the (masked-excluded) sim-hit count. Masked layers
          // are always excluded; merging is applied to the layer count.
          int nh;
          if (usePerLayerCountsForMeasAndHits) {
            std::set<std::pair<int, int>> lays;
            auto hl = hitLayers.find(key);
            if (hl != hitLayers.end()) {
              for (const auto& raw : hl->second) {
                if (!masked.count(raw)) lays.insert(resolveLayer(raw));
              }
            }
            nh = (int)lays.size();
          } else if (!masked.empty()) {
            auto hc = hitCounts.find(key);
            const long all = (hc != hitCounts.end()) ? hc->second.first : 0;
            const long used = (hc != hitCounts.end()) ? hc->second.second : 0;
            if (all != (long)(*pNhits)[i]) {
              std::cout << (strictCountCheck ? "ERROR" : "Warning")
                        << ": event " << *pEv << " number_of_hits mismatch: hits.root="
                        << all << " particles=" << (*pNhits)[i] << std::endl;
              if (strictCountCheck) {
                std::cout << "Aborting (strictCountCheck=true)." << std::endl;
                return;
              }
            }
            nh = (int)used;
          } else {
            nh = (*pNhits)[i];  // no masking: summary hit count
          }
          if (nh < effHitsMin) continue;

          ++nDen;
          hDenPt->Fill(dpt);
          hDenEta->Fill(deta);
          hDenPhi->Fill(dphi);

          auto it = recoMap.find(key);
          if (it == recoMap.end()) continue;  // not reconstructed
          ++nNum;
          hNumAllPt->Fill(dpt);
          hNumAllEta->Fill(deta);
          hNumAllPhi->Fill(dphi);
          const int k = it->second.nFake;
          if (k >= 0 && k <= maxFake) {  // k==0 = clean (0 fake hits)
            hNumExactPt[k]->Fill(dpt);
            hNumExactEta[k]->Fill(deta);
            hNumExactPhi[k]->Fill(dphi);
          }
          // "at most j" fake hits: a track with k fakes counts for all j >= k
          // (j == 0 therefore collects only clean tracks).
          for (int j = k; j <= maxFake; ++j) {
            hNumAtMostPt[j]->Fill(dpt);
            hNumAtMostEta[j]->Fill(deta);
            hNumAtMostPhi[j]->Fill(dphi);
          }
        }
      }

      std::cout << "efficiency denominator: particle-file events=" << nEffEvents
                << ", reconstructed events=" << recoEvents.size()
                << " -> denominator restricted to reco events (skipped "
                << nEffEvSkipped << " non-reconstructed)." << std::endl;

      effAllPt = makeRatio(hNumAllPt, hDenPt, "eff_all_vs_pT",
                           "all reconstructed;p_{T} [GeV];efficiency");
      effAllEta = makeRatio(hNumAllEta, hDenEta, "eff_all_vs_eta",
                            "all reconstructed;#eta;efficiency");
      effAllPhi = makeRatio(hNumAllPhi, hDenPhi, "eff_all_vs_phi",
                            "all reconstructed;#phi;efficiency");
      effExactPt[0] = makeRatio(hNumExactPt[0], hDenPt, "eff_exact0_vs_pT",
                                "0 fake (clean);p_{T} [GeV];efficiency");
      effExactEta[0] = makeRatio(hNumExactEta[0], hDenEta, "eff_exact0_vs_eta",
                                 "0 fake (clean);#eta;efficiency");
      effExactPhi[0] = makeRatio(hNumExactPhi[0], hDenPhi, "eff_exact0_vs_phi",
                                 "0 fake (clean);#phi;efficiency");
      effAtMostPt[0] = makeRatio(hNumAtMostPt[0], hDenPt, "eff_atMost0_vs_pT",
                                 "#leq 0 fake (clean);p_{T} [GeV];efficiency");
      effAtMostEta[0] = makeRatio(hNumAtMostEta[0], hDenEta, "eff_atMost0_vs_eta",
                                  "#leq 0 fake (clean);#eta;efficiency");
      effAtMostPhi[0] = makeRatio(hNumAtMostPhi[0], hDenPhi, "eff_atMost0_vs_phi",
                                  "#leq 0 fake (clean);#phi;efficiency");
      for (int k = 1; k <= maxFake; ++k) {
        effExactPt[k] = makeRatio(hNumExactPt[k], hDenPt, Form("eff_exact%d_vs_pT", k),
                                  Form("exactly %d fake;p_{T} [GeV];efficiency", k));
        effExactEta[k] = makeRatio(hNumExactEta[k], hDenEta, Form("eff_exact%d_vs_eta", k),
                                   Form("exactly %d fake;#eta;efficiency", k));
        effExactPhi[k] = makeRatio(hNumExactPhi[k], hDenPhi, Form("eff_exact%d_vs_phi", k),
                                   Form("exactly %d fake;#phi;efficiency", k));
        effAtMostPt[k] = makeRatio(hNumAtMostPt[k], hDenPt, Form("eff_atMost%d_vs_pT", k),
                                    Form("#leq %d fake;p_{T} [GeV];efficiency", k));
        effAtMostEta[k] = makeRatio(hNumAtMostEta[k], hDenEta, Form("eff_atMost%d_vs_eta", k),
                                     Form("#leq %d fake;#eta;efficiency", k));
        effAtMostPhi[k] = makeRatio(hNumAtMostPhi[k], hDenPhi, Form("eff_atMost%d_vs_phi", k),
                                     Form("#leq %d fake;#phi;efficiency", k));
      }
      std::cout << "Efficiency: denominator particles = " << nDen
                << ", reconstructed = " << nNum;
      if (nDen > 0) std::cout << "  (" << 100.0 * nNum / nDen << " %)";
      std::cout << "  [denominator: " << particlesEffFile
                << ", number_of_hits >= " << effHitsMin << "]" << std::endl;
      haveEff = true;
    }
  }

  // ------------------------------------------------------------------
  // 4c) "Physics" efficiency: denominator = all CHARGED primaries from the
  //     generator file (particlesGenFile). Only the pid + pT/eta acceptance
  //     cuts are applied -- NO nHits / required-layer / primaryOnly cuts, since
  //     the generator file already defines the primaries (PYTHIA: particles with
  //     ctau0 > 1 cm are stable). Same numerator as the algorithmic efficiency
  //     (matched to a selected track), split by the matching track's fake count.
  //     Drawn as OPEN markers on the efficiency canvases.
  // ------------------------------------------------------------------
  bool haveGen = false;
  TH1D* effPhysAllPt = nullptr;
  TH1D* effPhysAllEta = nullptr;
  TH1D* effPhysAllPhi = nullptr;
  std::vector<TH1D*> effPhysExactPt(maxFake + 1, nullptr), effPhysExactEta(maxFake + 1, nullptr);
  std::vector<TH1D*> effPhysAtMostPt(maxFake + 1, nullptr), effPhysAtMostEta(maxFake + 1, nullptr);
  std::vector<TH1D*> effPhysExactPhi(maxFake + 1, nullptr), effPhysAtMostPhi(maxFake + 1, nullptr);
  {
    std::string genPath = std::string(dirSim) + "/" + particlesGenFile;
    TFile* fGen = TFile::Open(genPath.c_str(), "READ");
    if (fGen == nullptr || fGen->IsZombie()) {
      std::cout << "WARNING: cannot open " << genPath
                << " -> physics-efficiency curves skipped." << std::endl;
    } else {
      auto makeEtaHist = [&](const char* name) {
        TH1D* h = new TH1D(name, ";#eta;particles", etaNbins, etaMin, etaMax);
        h->SetDirectory(nullptr);
        return h;
      };
      TH1D* hGenDenPt = makePtHist("hGenDenPt");
      TH1D* hGenDenEta = makeEtaHist("hGenDenEta");
      TH1D* hGenDenPhi = makePhiHist("hGenDenPhi");
      TH1D* hGenNumAllPt = makePtHist("hGenNumAllPt");
      TH1D* hGenNumAllEta = makeEtaHist("hGenNumAllEta");
      TH1D* hGenNumAllPhi = makePhiHist("hGenNumAllPhi");
      std::vector<TH1D*> hGenExactPt(maxFake + 1, nullptr), hGenExactEta(maxFake + 1, nullptr);
      std::vector<TH1D*> hGenAtMostPt(maxFake + 1, nullptr), hGenAtMostEta(maxFake + 1, nullptr);
      std::vector<TH1D*> hGenExactPhi(maxFake + 1, nullptr), hGenAtMostPhi(maxFake + 1, nullptr);
      for (int k = 0; k <= maxFake; ++k) {
        hGenExactPt[k] = makePtHist(Form("hGenExactPt_%d", k));
        hGenExactEta[k] = makeEtaHist(Form("hGenExactEta_%d", k));
        hGenAtMostPt[k] = makePtHist(Form("hGenAtMostPt_%d", k));
        hGenAtMostEta[k] = makeEtaHist(Form("hGenAtMostEta_%d", k));
        hGenExactPhi[k] = makePhiHist(Form("hGenExactPhi_%d", k));
        hGenAtMostPhi[k] = makePhiHist(Form("hGenAtMostPhi_%d", k));
      }

      TTreeReader rGen(particlesGenTree, fGen);
      TTreeReaderValue<std::uint32_t> gEv(rGen, "event_id");
      TTreeReaderValue<std::vector<std::int32_t>> gType(rGen, "particle_type");
      TTreeReaderValue<std::vector<float>> gPt(rGen, "pt");
      TTreeReaderValue<std::vector<float>> gEta(rGen, "eta");
      TTreeReaderValue<std::vector<float>> gPhi(rGen, "phi");
      TTreeReaderValue<std::vector<float>> gVz(rGen, "vz");  // production-vertex z [mm]
      TTreeReaderValue<std::vector<float>> gQ(rGen, "q");
      TTreeReaderValue<std::vector<std::uint32_t>> gVP(rGen, "vertex_primary");
      TTreeReaderValue<std::vector<std::uint32_t>> gVS(rGen, "vertex_secondary");
      TTreeReaderValue<std::vector<std::uint32_t>> gPa(rGen, "particle");
      TTreeReaderValue<std::vector<std::uint32_t>> gGe(rGen, "generation");
      TTreeReaderValue<std::vector<std::uint32_t>> gSu(rGen, "sub_particle");

      long nGenDen = 0, nGenNum = 0;
      long nGenEvents = 0, nGenEvSkipped = 0;
      while (rGen.Next()) {
        ++nGenEvents;
        // Restrict to reconstructed events (see the algorithmic-efficiency loop).
        if (!recoEvents.count(*gEv)) { ++nGenEvSkipped; continue; }
        for (std::size_t i = 0; i < gType->size(); ++i) {
          if (std::abs((*gQ)[i]) < 0.5f) continue;  // charged only
          if (!passesPid((*gType)[i])) continue;
          if ((*gVz)[i] < zVtxMinMm || (*gVz)[i] > zVtxMaxMm) continue;  // collision-vertex z
          const float dpt = (*gPt)[i];
          const float deta = (*gEta)[i];
          const float dphi = wrapPhi((*gPhi)[i]);
          if (dpt < ptMin || dpt > ptMax) continue;
          if (deta < etaMin || deta > etaMax) continue;
          if (dphi < phiMin || dphi > phiMax) continue;

          ++nGenDen;
          hGenDenPt->Fill(dpt);
          hGenDenEta->Fill(deta);
          hGenDenPhi->Fill(dphi);

          PKey key{*gEv, (*gVP)[i], (*gVS)[i], (*gPa)[i], (*gGe)[i], (*gSu)[i]};
          auto it = recoMap.find(key);
          if (it == recoMap.end()) continue;
          ++nGenNum;
          hGenNumAllPt->Fill(dpt);
          hGenNumAllEta->Fill(deta);
          hGenNumAllPhi->Fill(dphi);
          const int k = it->second.nFake;
          if (k >= 0 && k <= maxFake) {
            hGenExactPt[k]->Fill(dpt);
            hGenExactEta[k]->Fill(deta);
            hGenExactPhi[k]->Fill(dphi);
          }
          for (int j = k; j <= maxFake; ++j) {
            hGenAtMostPt[j]->Fill(dpt);
            hGenAtMostEta[j]->Fill(deta);
            hGenAtMostPhi[j]->Fill(dphi);
          }
        }
      }

      std::cout << "physics-efficiency denominator: generator-file events=" << nGenEvents
                << ", reconstructed events=" << recoEvents.size()
                << " -> denominator restricted to reco events (skipped "
                << nGenEvSkipped << " non-reconstructed)." << std::endl;

      effPhysAllPt = makeRatio(hGenNumAllPt, hGenDenPt, "physeff_all_vs_pT",
                               "all reconstructed (phys);p_{T} [GeV];efficiency");
      effPhysAllEta = makeRatio(hGenNumAllEta, hGenDenEta, "physeff_all_vs_eta",
                                "all reconstructed (phys);#eta;efficiency");
      effPhysAllPhi = makeRatio(hGenNumAllPhi, hGenDenPhi, "physeff_all_vs_phi",
                                "all reconstructed (phys);#phi;efficiency");
      for (int k = 0; k <= maxFake; ++k) {
        effPhysExactPt[k] = makeRatio(hGenExactPt[k], hGenDenPt, Form("physeff_exact%d_vs_pT", k),
                                      Form("%d fake (phys);p_{T} [GeV];efficiency", k));
        effPhysExactEta[k] = makeRatio(hGenExactEta[k], hGenDenEta, Form("physeff_exact%d_vs_eta", k),
                                       Form("%d fake (phys);#eta;efficiency", k));
        effPhysExactPhi[k] = makeRatio(hGenExactPhi[k], hGenDenPhi, Form("physeff_exact%d_vs_phi", k),
                                       Form("%d fake (phys);#phi;efficiency", k));
        effPhysAtMostPt[k] = makeRatio(hGenAtMostPt[k], hGenDenPt, Form("physeff_atMost%d_vs_pT", k),
                                       Form("#leq %d fake (phys);p_{T} [GeV];efficiency", k));
        effPhysAtMostEta[k] = makeRatio(hGenAtMostEta[k], hGenDenEta, Form("physeff_atMost%d_vs_eta", k),
                                        Form("#leq %d fake (phys);#eta;efficiency", k));
        effPhysAtMostPhi[k] = makeRatio(hGenAtMostPhi[k], hGenDenPhi, Form("physeff_atMost%d_vs_phi", k),
                                        Form("#leq %d fake (phys);#phi;efficiency", k));
      }
      std::cout << "Physics efficiency: charged generator primaries = " << nGenDen
                << ", reconstructed = " << nGenNum;
      if (nGenDen > 0) std::cout << "  (" << 100.0 * nGenNum / nGenDen << " %)";
      std::cout << "  [denominator: " << particlesGenFile << "]" << std::endl;
      haveGen = true;
    }
  }

  // ------------------------------------------------------------------
  // 4d) Resolution & bias: per-x-bin Gaussian fit of each residual TH2 within
  //     +/-3 RMS -> sigma(X) (resolution) and mean(X) (bias). resWarn[p] is set
  //     if the configured residual range is narrower than 4*RMS.
  // ------------------------------------------------------------------
  auto fitSlices = [&](TH2D* h2, const char* base, TH1D*& hSig, TH1D*& hMean) {
    hSig = (TH1D*)h2->ProjectionX(Form("%s_sigma", base));
    hMean = (TH1D*)h2->ProjectionX(Form("%s_mean", base));
    hSig->Reset("ICES");
    hMean->Reset("ICES");
    hSig->SetDirectory(nullptr);
    hMean->SetDirectory(nullptr);
    for (int bx = 1; bx <= h2->GetNbinsX(); ++bx) {
      TH1D* py = h2->ProjectionY(Form("%s_py%d", base, bx), bx, bx);
      if (py->GetEntries() >= 20) {
        const double mu0 = py->GetMean();
        const double rms = py->GetRMS();
        if (rms > 0) {
          TF1 g("g", "gaus", mu0 - 3 * rms, mu0 + 3 * rms);
          py->Fit(&g, "QNR");
          hSig->SetBinContent(bx, g.GetParameter(2));
          hSig->SetBinError(bx, g.GetParError(2));
          hMean->SetBinContent(bx, g.GetParameter(1));
          hMean->SetBinError(bx, g.GetParError(1));
        }
      }
      delete py;
    }
  };

  std::vector<std::vector<std::vector<TH1D*>>> resSig(
      3, std::vector<std::vector<TH1D*>>(3, std::vector<TH1D*>(nResCat, nullptr)));
  std::vector<std::vector<std::vector<TH1D*>>> resMean(
      3, std::vector<std::vector<TH1D*>>(3, std::vector<TH1D*>(nResCat, nullptr)));
  bool resWarn[3] = {false, false, false};
  for (int p = 0; p < 3; ++p) {
    TH1D* py = res2[p][0][0]->ProjectionY(Form("warnpy_%d", p));
    if (py->GetEntries() > 0 && 2.0 * resParMax[p] < 4.0 * py->GetRMS()) resWarn[p] = true;
    delete py;
    for (int a = 0; a < 3; ++a) {
      for (int c = 0; c < nResCat; ++c) {
        fitSlices(res2[p][a][c], Form("res_%s_%s_c%d", resParName[p], resAxName[a], c),
                  resSig[p][a][c], resMean[p][a][c]);
        resSig[p][a][c]->GetXaxis()->SetTitle(resAxTitle[a]);
        resMean[p][a][c]->GetXaxis()->SetTitle(resAxTitle[a]);
        resSig[p][a][c]->SetTitle("");
        resMean[p][a][c]->SetTitle("");
      }
    }
  }

  // ------------------------------------------------------------------
  // 5) Draw and save.
  // ------------------------------------------------------------------
  const int colors[] = {kBlack, kRed + 1, kBlue + 1, kGreen + 2, kMagenta + 1, kOrange + 7};
  auto drawSet = [&](std::vector<TH1D*>& hs, const char* cname, const char* ctitle,
                     bool logx) {
    TCanvas* c = new TCanvas(cname, ctitle, 900, 650);
    c->SetGrid();
    if (logx) c->SetLogx();
    TLegend* leg = new TLegend(0.65, 0.70, 0.88, 0.88);
    bool first = true;
    for (int k = 1; k <= maxFake; ++k) {
      TH1D* h = hs[k];
      h->SetLineColor(colors[k % 6]);
      h->SetMarkerColor(colors[k % 6]);
      h->SetMarkerStyle(20);
      h->SetMinimum(0.0);
      if (logx) {
        h->GetXaxis()->SetMoreLogLabels();
        h->GetXaxis()->SetNoExponent();
      }
      h->Draw(first ? "E1" : "E1 SAME");
      leg->AddEntry(h, h->GetTitle(), "lp");
      first = false;
    }
    leg->Draw();
    return c;
  };

  // pT canvases honor the logPt option; eta canvases are always linear.
  TCanvas* cExPt = drawSet(rExactPt, "cExactPt", "exactly-k fake hits vs pT", logPt);
  TCanvas* cExEt = drawSet(rExactEta, "cExactEta", "exactly-k fake hits vs eta", false);
  TCanvas* cExPh = drawSet(rExactPhi, "cExactPhi", "exactly-k fake hits vs phi", false);
  TCanvas* cAlPt = drawSet(rAtLeastPt, "cAtLeastPt", "at-least-k fake hits vs pT", logPt);
  TCanvas* cAlEt = drawSet(rAtLeastEta, "cAtLeastEta", "at-least-k fake hits vs eta", false);
  TCanvas* cAlPh = drawSet(rAtLeastPhi, "cAtLeastPhi", "at-least-k fake hits vs phi", false);

  // Fake-origin canvases: overlay the primary / material / decay curves, with the
  // "#geq 1 fake (any)" curve as a black reference. ref may be nullptr.
  const int kOriginColor[4] = {kRed + 1, kBlue + 1, kGreen + 2, kOrange + 7};  // primary, material, decay, noise
  auto drawOrigin = [&](std::vector<TH1D*>& hs, TH1D* ref, const char* cname,
                        const char* ctitle, bool logx) -> TCanvas* {
    TCanvas* c = new TCanvas(cname, ctitle, 900, 650);
    c->SetGrid();
    if (logx) c->SetLogx();
    TLegend* leg = new TLegend(0.62, 0.70, 0.88, 0.88);
    bool first = true;
    if (ref != nullptr) {
      // Clone the reference so restyling it here does NOT mutate the shared
      // at-least-k histogram (which is drawn with its own colors elsewhere).
      TH1D* refC = (TH1D*)ref->Clone(Form("%s_ref", cname));
      refC->SetDirectory(nullptr);
      refC->SetLineColor(kBlack);
      refC->SetMarkerColor(kBlack);
      refC->SetMarkerStyle(24);  // open circle for the "any" reference
      refC->SetMinimum(0.0);
      if (logx) { refC->GetXaxis()->SetMoreLogLabels(); refC->GetXaxis()->SetNoExponent(); }
      refC->Draw("E1");
      leg->AddEntry(refC, "#geq 1 fake (any)", "lp");
      first = false;
    }
    for (int o = 0; o < 4; ++o) {
      if (o == 3 && !drawNoise) continue;  // 3 = injected noise
      TH1D* h = hs[o];
      h->SetLineColor(kOriginColor[o]);
      h->SetMarkerColor(kOriginColor[o]);
      h->SetMarkerStyle(20);
      h->SetMinimum(0.0);
      if (logx) { h->GetXaxis()->SetMoreLogLabels(); h->GetXaxis()->SetNoExponent(); }
      h->Draw(first ? "E1" : "E1 SAME");
      leg->AddEntry(h, Form("from %s", kOriginName[o]), "lp");
      first = false;
    }
    leg->Draw();
    return c;
  };

  TCanvas *cOrigPt = nullptr, *cOrigEt = nullptr, *cOrigPh = nullptr;
  if (drawFakeOrigins) {
    TH1D* refPt = (maxFake >= 1) ? rAtLeastPt[1] : nullptr;
    TH1D* refEta = (maxFake >= 1) ? rAtLeastEta[1] : nullptr;
    TH1D* refPhi = (maxFake >= 1) ? rAtLeastPhi[1] : nullptr;
    cOrigPt = drawOrigin(rOrigPt, refPt, "cOriginPt", "fake-origin rate vs pT", logPt);
    cOrigEt = drawOrigin(rOrigEta, refEta, "cOriginEta", "fake-origin rate vs eta", false);
    cOrigPh = drawOrigin(rOrigPhi, refPhi, "cOriginPhi", "fake-origin rate vs phi", false);
  }

  // Efficiency canvases: overlay 'all reconstructed' with the per-fake-count
  // curves (filled markers = reconstructable/algorithmic efficiency). If the
  // generator file is available, the physics-efficiency twins are overlaid with
  // the SAME colors but OPEN markers.
  auto effDraw = [&](TH1D* hAll, std::vector<TH1D*>& hs,
                     TH1D* hAllPhys, std::vector<TH1D*>& hsPhys, bool drawPhys,
                     const char* cname, const char* ctitle, bool logx, int kStart) {
    TCanvas* c = new TCanvas(cname, ctitle, 900, 650);
    c->SetGrid();
    if (logx) c->SetLogx();
    TLegend* leg = new TLegend(0.65, 0.18, 0.88, 0.42);
    hAll->SetLineColor(colors[0]);
    hAll->SetMarkerColor(colors[0]);
    hAll->SetMarkerStyle(20);
    hAll->SetMinimum(0.0);
    hAll->SetMaximum(1.05);
    if (logx) {
      hAll->GetXaxis()->SetMoreLogLabels();
      hAll->GetXaxis()->SetNoExponent();
    }
    hAll->Draw("E1");
    leg->AddEntry(hAll, "all reconstructed", "lp");
    for (int k = kStart; k <= maxFake; ++k) {
      TH1D* h = hs[k];
      const int col = colors[(k + 1) % 6];  // +1 so k==0 is not black like 'all'
      h->SetLineColor(col);
      h->SetMarkerColor(col);
      h->SetMarkerStyle(20);
      h->Draw("E1 SAME");
      leg->AddEntry(h, h->GetTitle(), "lp");
    }
    if (drawPhys && hAllPhys != nullptr) {
      hAllPhys->SetLineColor(colors[0]);
      hAllPhys->SetMarkerColor(colors[0]);
      hAllPhys->SetMarkerStyle(24);  // open circle
      hAllPhys->Draw("E1 SAME");
      for (int k = kStart; k <= maxFake; ++k) {
        TH1D* h = hsPhys[k];
        const int col = colors[(k + 1) % 6];
        h->SetLineColor(col);
        h->SetMarkerColor(col);
        h->SetMarkerStyle(24);
        h->Draw("E1 SAME");
      }
      TLatex tl;
      tl.SetNDC();
      tl.SetTextSize(0.03);
      tl.DrawLatex(0.13, 0.91,
                   "filled: reconstructable eff    open: physics eff");
    }
    leg->Draw();
    return c;
  };

  TCanvas *cEffExPt = nullptr, *cEffExEt = nullptr, *cEffExPh = nullptr;
  TCanvas *cEffAmPt = nullptr, *cEffAmEt = nullptr, *cEffAmPh = nullptr;
  TCanvas* cEffBreakPt = nullptr;  // example-style fake-track breakdown vs pT
  if (haveEff) {
    // exactly-k / at-most canvases include k==0 (clean tracks)
    cEffExPt = effDraw(effAllPt, effExactPt, effPhysAllPt, effPhysExactPt, haveGen,
                       "cEffExactPt", "efficiency (exactly-k) vs pT", logPt, 0);
    cEffExEt = effDraw(effAllEta, effExactEta, effPhysAllEta, effPhysExactEta, haveGen,
                       "cEffExactEta", "efficiency (exactly-k) vs eta", false, 0);
    cEffExPh = effDraw(effAllPhi, effExactPhi, effPhysAllPhi, effPhysExactPhi, haveGen,
                       "cEffExactPhi", "efficiency (exactly-k) vs phi", false, 0);
    cEffAmPt = effDraw(effAllPt, effAtMostPt, effPhysAllPt, effPhysAtMostPt, haveGen,
                       "cEffAtMostPt", "efficiency (<=k) vs pT", logPt, 0);
    cEffAmEt = effDraw(effAllEta, effAtMostEta, effPhysAllEta, effPhysAtMostEta, haveGen,
                       "cEffAtMostEta", "efficiency (<=k) vs eta", false, 0);
    cEffAmPh = effDraw(effAllPhi, effAtMostPhi, effPhysAllPhi, effPhysAtMostPhi, haveGen,
                       "cEffAtMostPhi", "efficiency (<=k) vs phi", false, 0);

    // ------------------------------------------------------------------
    // Example-style "fake-track breakdown" vs pT:
    //   Total       = all reconstructed (effAll, black line)
    //   Good tracks = 0 fake hits       (effExact[0], blue)
    //   1 fake      = exactly 1 fake    (effExact[1], red)
    //   2 fake      = exactly 2 fakes   (effExact[2], gray hatched fill)
    //   (higher k, if maxFake>2, added with extra colors)
    // ------------------------------------------------------------------
    cEffBreakPt = new TCanvas("cEffFakeBreakdownPt", "fake-track breakdown vs pT", 900, 650);
    cEffBreakPt->SetGrid();
    if (logPt) cEffBreakPt->SetLogx();
    TLegend* legB = new TLegend(0.30, 0.78, 0.88, 0.90);
    legB->SetNColumns(2);
    legB->SetBorderSize(0);
    legB->SetFillStyle(0);

    // Total (black) sets the frame.
    effAllPt->SetLineColor(kBlack);
    effAllPt->SetLineWidth(2);
    effAllPt->SetMarkerStyle(1);
    effAllPt->SetMinimum(0.0);
    effAllPt->SetMaximum(1.15);
    effAllPt->SetTitle(";p_{T} [GeV/#it{c}];Efficiency (fake-track rate)");
    if (logPt) { effAllPt->GetXaxis()->SetMoreLogLabels(); effAllPt->GetXaxis()->SetNoExponent(); }
    effAllPt->Draw("HIST");
    legB->AddEntry(effAllPt, "Total", "l");

    // 1-fake and 2-fake as hatched fills (drawn before the Good/Total lines so
    // the lines stay on top); 1-fake is larger -> drawn first, 2-fake on top.
    if (maxFake >= 1) {
      effExactPt[1]->SetLineColor(kRed + 1);
      effExactPt[1]->SetFillColor(kRed + 1);
      effExactPt[1]->SetFillStyle(3354);
      effExactPt[1]->SetMarkerStyle(1);
      effExactPt[1]->Draw("HIST SAME");
    }
    if (maxFake >= 2) {
      effExactPt[2]->SetLineColor(kGray + 2);
      effExactPt[2]->SetFillColor(kGray + 1);
      effExactPt[2]->SetFillStyle(3345);
      effExactPt[2]->SetMarkerStyle(1);
      effExactPt[2]->Draw("HIST SAME");
    }
    // Good tracks (0 fakes).
    effExactPt[0]->SetLineColor(kBlue + 1);
    effExactPt[0]->SetLineWidth(2);
    effExactPt[0]->SetMarkerStyle(1);
    effExactPt[0]->Draw("HIST E SAME");
    legB->AddEntry(effExactPt[0], "Good tracks (0 fake)", "l");
    if (maxFake >= 1) legB->AddEntry(effExactPt[1], "Tracks with 1 fake cluster", "f");
    if (maxFake >= 2) legB->AddEntry(effExactPt[2], "Tracks with 2 fake clusters", "f");
    // higher fake counts (if any).
    const int kBreakCol[6] = {kBlue + 1, kRed + 1, kGray + 2, kMagenta + 1, kOrange + 7, kAzure + 2};
    for (int k = 3; k <= maxFake; ++k) {
      effExactPt[k]->SetLineColor(kBreakCol[k % 6]);
      effExactPt[k]->SetLineWidth(2);
      effExactPt[k]->SetMarkerStyle(1);
      effExactPt[k]->Draw("HIST E SAME");
      legB->AddEntry(effExactPt[k], Form("Tracks with %d fake clusters", k), "l");
    }
    effAllPt->Draw("HIST SAME");  // keep the Total line on top
    legB->Draw();
  }

  TCanvas* cLayer = new TCanvas("cLayer", "per-layer wrong-hit rate", 1100, 650);
  cLayer->SetBottomMargin(0.4);  // room for the (vertical) per-layer labels
  cLayer->SetGridy();            // horizontal grid lines
  hLayerRate->SetMinimum(0.0);
  hLayerRate->SetFillColor(kAzure - 9);
  hLayerRate->LabelsOption("v");  // vertical bin labels
  hLayerRate->Draw("HIST E");
  // Overlay the per-layer origin breakdown (primary / material / decay) as
  // markers on top of the total wrong-hit bars. Bins line up because all four
  // histograms share the crossing denominator and the same drop/order logic.
  if (drawFakeOrigins) {
    TLegend* legL = new TLegend(0.70, 0.72, 0.88, 0.88);
    legL->AddEntry(hLayerRate, "all (any origin)", "f");
    for (int o = 0; o < 4; ++o) {
      if (o == 3 && !drawNoise) continue;  // 3 = injected noise
      TH1D* h = hLayerRateOrig[o];
      h->SetLineColor(kOriginColor[o]);
      h->SetMarkerColor(kOriginColor[o]);
      h->SetMarkerStyle(20 + o);  // 20=circle,21=square,22=triangle,23=down-triangle(noise)
      h->SetMarkerSize(0.9);
      h->Draw("E1 X0 SAME");
      legL->AddEntry(h, kOriginName[o], "lp");
    }
    legL->Draw();
  }

  // Per-layer origin FRACTIONS: share of each layer's wrong hits from primary /
  // material / decay (numerator = origin count, denominator = total wrong on the
  // layer). Origins can overlap and "unknown" is omitted, so they need not sum to 1.
  TCanvas* cLayerOrigFrac = nullptr;
  if (drawFakeOrigins) {
    cLayerOrigFrac = new TCanvas("cLayerOrigFrac", "per-layer origin fractions", 1100, 650);
    cLayerOrigFrac->SetBottomMargin(0.4);
    cLayerOrigFrac->SetGridy();
    TLegend* legOF = new TLegend(0.70, 0.72, 0.88, 0.88);
    bool firstOF = true;
    for (int o = 0; o < 4; ++o) {
      if (o == 3 && !drawNoise) continue;  // 3 = injected noise
      TH1D* h = hLayerOrigFrac[o];
      h->SetTitle("per-layer origin fractions;;fraction of wrong hits");
      h->SetMinimum(0.0);
      h->SetMaximum(1.05);
      h->SetLineColor(kOriginColor[o]);
      h->SetMarkerColor(kOriginColor[o]);
      h->SetMarkerStyle(20 + o);
      h->SetMarkerSize(0.9);
      h->LabelsOption("v");
      h->Draw(firstOF ? "E1 X0" : "E1 X0 SAME");
      legOF->AddEntry(h, kOriginName[o], "lp");
      firstOF = false;
    }
    legOF->Draw();
  }

  TCanvas* cLayerHole = new TCanvas("cLayerHole", "per-layer hole rate", 1100, 650);
  cLayerHole->SetBottomMargin(0.4);
  cLayerHole->SetGridy();  // horizontal grid lines
  hLayerHoleRate->SetMinimum(0.0);
  hLayerHoleRate->SetFillColor(kOrange - 9);
  hLayerHoleRate->LabelsOption("v");
  hLayerHoleRate->Draw("HIST E");

  // Per-layer mean contributor multiplicity.
  TCanvas* cLayerContrib = nullptr;
  if (drawContribMult) {
    cLayerContrib =
        new TCanvas("cLayerContrib", "per-layer mean contributors per hit", 1100, 650);
    cLayerContrib->SetBottomMargin(0.4);
    cLayerContrib->SetGridy();
    hLayerContrib->SetMinimum(0.99);
    hLayerContrib->SetMaximum(1.03);
    hLayerContrib->SetMarkerStyle(20);
    hLayerContrib->SetMarkerSize(0.9);
    hLayerContrib->SetMarkerColor(kAzure + 2);
    hLayerContrib->SetLineColor(kAzure + 2);
    hLayerContrib->LabelsOption("v");
    hLayerContrib->Draw("E1 X0");
  }

  // Per-layer MATERIAL breakdown. Two canvases:
  //  cLayerMat     : material total (bars) + the three sub-processes (markers),
  //                  all as wrong/crossing rates so they overlay directly;
  //  cLayerMatFrac : fraction of each layer's material fakes from each
  //                  sub-process (sub / total material on the layer).
  // conv, brems, nuclear, delta-ray, other
  const int kMatSubColor[5] = {kViolet + 1, kOrange + 7, kTeal + 2, kRed + 1, kGray + 2};
  TCanvas *cLayerMat = nullptr, *cLayerMatFrac = nullptr;
  // combined (rate top panel + fraction bottom panel) versions:
  TCanvas *cLayerRateFrac = nullptr, *cLayerMatRateFrac = nullptr;
  if (drawFakeOrigins) {
    cLayerMat = new TCanvas("cLayerMat", "per-layer material wrong-hit rate", 1100, 650);
    cLayerMat->SetBottomMargin(0.4);
    cLayerMat->SetGridy();
    // Clone the material-total hist so styling it as bars here does not disturb
    // its marker styling on the main per-layer canvas.
    TH1D* hMatTot = (TH1D*)hLayerRateOrig[1]->Clone("fakeRate_perLayer_materialTotal");
    hMatTot->SetDirectory(nullptr);
    hMatTot->SetTitle("per-layer material wrong-hit rate;;wrong / crossing tracks");
    hMatTot->SetMinimum(0.0);
    hMatTot->SetFillColor(kAzure - 9);
    hMatTot->SetLineColor(kAzure - 5);
    hMatTot->LabelsOption("v");
    hMatTot->Draw("HIST E");
    TLegend* legM = new TLegend(0.66, 0.70, 0.88, 0.88);
    legM->AddEntry(hMatTot, "material (total)", "f");
    for (int s = 0; s < kNMatSub; ++s) {
      if (s == 1) continue;  // brems folded into "other" - not drawn
      TH1D* h = hLayerMatSub[s];
      h->SetLineColor(kMatSubColor[s]);
      h->SetMarkerColor(kMatSubColor[s]);
      h->SetMarkerStyle(20 + s);
      h->SetMarkerSize(0.9);
      h->Draw("E1 X0 SAME");
      legM->AddEntry(h, kMatSubName[s], "lp");
    }
    legM->Draw();

    cLayerMatFrac = new TCanvas("cLayerMatFrac", "per-layer material fractions", 1100, 650);
    cLayerMatFrac->SetBottomMargin(0.4);
    cLayerMatFrac->SetGridy();
    TLegend* legMF = new TLegend(0.66, 0.70, 0.88, 0.88);
    bool firstMF = true;
    for (int s = 0; s < kNMatSub; ++s) {
      if (s == 1) continue;  // brems folded into "other" - not drawn
      TH1D* h = hLayerMatFrac[s];
      h->SetTitle("per-layer material fractions;;fraction of material fakes");
      h->SetMinimum(0.0);
      h->SetMaximum(1.05);
      h->SetLineColor(kMatSubColor[s]);
      h->SetMarkerColor(kMatSubColor[s]);
      h->SetMarkerStyle(20 + s);
      h->SetMarkerSize(0.9);
      h->LabelsOption("v");
      h->Draw(firstMF ? "E1 X0" : "E1 X0 SAME");
      legMF->AddEntry(h, kMatSubName[s], "lp");
      firstMF = false;
    }
    legMF->Draw();

    // ---- combined canvases: per-layer RATE (top panel) + FRACTIONS (bottom
    //      panel), i.e. the same content as the two separate canvases above but
    //      stacked as ratio-style panels sharing the x axis. ----
    auto makeRateFracCanvas =
        [&](const char* cname, const char* ctitle, TH1D* barsSrc,
            const char* barsLabel, std::vector<TH1D*>& overlays,
            std::vector<TH1D*>& fracs, const int* cols, const char** names,
            const char* topY, const char* botY) -> TCanvas* {
      TCanvas* c = new TCanvas(cname, ctitle, 1100, 850);
      TPad* pTop = new TPad(Form("%s_top", cname), "", 0, 0.42, 1, 1);
      TPad* pBot = new TPad(Form("%s_bot", cname), "", 0, 0.0, 1, 0.42);
      pTop->SetBottomMargin(0.02);
      pTop->SetGridy();
      pBot->SetTopMargin(0.02);
      pBot->SetBottomMargin(0.45);  // room for the vertical per-layer labels
      pBot->SetGridy();
      pTop->Draw();
      pBot->Draw();

      // top: total wrong-hit rate (bars) + per-origin markers, labels hidden.
      pTop->cd();
      TH1D* bars = (TH1D*)barsSrc->Clone(Form("%s_bars", cname));
      bars->SetDirectory(nullptr);
      bars->SetTitle(Form("%s;;%s", ctitle, topY));
      bars->SetMinimum(0.0);
      bars->SetFillColor(kAzure - 9);
      bars->SetLineColor(kAzure - 5);
      bars->GetXaxis()->SetLabelSize(0.0);
      bars->Draw("HIST E");
      TLegend* legT = new TLegend(0.70, 0.58, 0.88, 0.88);
      legT->AddEntry(bars, barsLabel, "f");
      for (std::size_t o = 0; o < overlays.size(); ++o) {
        TH1D* h = (TH1D*)overlays[o]->Clone(Form("%s_ov%zu", cname, o));
        h->SetDirectory(nullptr);
        h->SetLineColor(cols[o]);
        h->SetMarkerColor(cols[o]);
        h->SetMarkerStyle(20 + (int)o);
        h->SetMarkerSize(0.9);
        h->GetXaxis()->SetLabelSize(0.0);
        h->Draw("E1 X0 SAME");
        legT->AddEntry(h, names[o], "lp");
      }
      legT->Draw();

      // bottom: the corresponding fractions, with the (vertical) layer labels.
      pBot->cd();
      bool first = true;
      for (std::size_t o = 0; o < fracs.size(); ++o) {
        TH1D* h = (TH1D*)fracs[o]->Clone(Form("%s_fr%zu", cname, o));
        h->SetDirectory(nullptr);
        h->SetTitle(Form(";;%s", botY));
        h->SetMinimum(0.0);
        h->SetMaximum(1.05);
        h->SetLineColor(cols[o]);
        h->SetMarkerColor(cols[o]);
        h->SetMarkerStyle(20 + (int)o);
        h->SetMarkerSize(0.9);
        h->LabelsOption("v");
        // ratio-panel axis cosmetics
        h->GetXaxis()->SetLabelSize(0.10);
        h->GetYaxis()->SetLabelSize(0.06);
        h->GetYaxis()->SetTitleSize(0.06);
        h->GetYaxis()->SetTitleOffset(0.8);
        h->Draw(first ? "E1 X0" : "E1 X0 SAME");
        first = false;
      }
      return c;
    };

    // Reduced series for the combined canvases: drop the injected-noise origin
    // (unless drawNoise) and always drop bremsstrahlung (folded into "other").
    std::vector<TH1D*> ovOrigR, ovOrigF;
    std::vector<int> colOrig; std::vector<const char*> namOrig;
    for (int o = 0; o < 4; ++o) {
      if (o == 3 && !drawNoise) continue;  // 3 = noise
      ovOrigR.push_back(hLayerRateOrig[o]); ovOrigF.push_back(hLayerOrigFrac[o]);
      colOrig.push_back(kOriginColor[o]);   namOrig.push_back(kOriginName[o]);
    }
    std::vector<TH1D*> ovMatR, ovMatF;
    std::vector<int> colMat; std::vector<const char*> namMat;
    for (int s = 0; s < kNMatSub; ++s) {
      if (s == 1) continue;  // brems folded into "other"
      ovMatR.push_back(hLayerMatSub[s]); ovMatF.push_back(hLayerMatFrac[s]);
      colMat.push_back(kMatSubColor[s]);  namMat.push_back(kMatSubName[s]);
    }

    cLayerRateFrac = makeRateFracCanvas(
        "cLayerRateFrac", "per-layer wrong-hit rate + origin fractions",
        hLayerRate, "all (any origin)", ovOrigR, ovOrigF,
        colOrig.data(), namOrig.data(), "wrong / crossing tracks",
        "fraction of wrong hits");

    cLayerMatRateFrac = makeRateFracCanvas(
        "cLayerMatRateFrac", "per-layer material rate + sub-process fractions",
        hLayerRateOrig[1], "material (total)", ovMatR, ovMatF,
        colMat.data(), namMat.data(), "wrong / crossing tracks",
        "fraction of material fakes");
  }

  // Resolution / bias canvases: overlay 'all matched' + '<=k' curves.
  auto resDraw = [&](std::vector<TH1D*>& hs, const char* cname, const char* ctitle,
                     const char* ytitle, bool logx, bool warn, bool symRange) {
    TCanvas* c = new TCanvas(cname, ctitle, 900, 650);
    c->SetGrid();
    if (logx) c->SetLogx();
    TLegend* leg = new TLegend(0.65, 0.68, 0.88, 0.88);
    double ymax = 0.0, ymin = 0.0;
    for (int cc = 0; cc < nResCat; ++cc) {
      ymax = std::max(ymax, hs[cc]->GetMaximum());
      ymin = std::min(ymin, hs[cc]->GetMinimum());
    }
    bool first = true;
    for (int cc = 0; cc < nResCat; ++cc) {
      TH1D* h = hs[cc];
      const int col = (cc == 0) ? colors[0] : colors[cc % 6];
      h->SetLineColor(col);
      h->SetMarkerColor(col);
      h->SetMarkerStyle(20);
      h->GetYaxis()->SetTitle(ytitle);
      if (symRange) {
        double m = 1.2 * std::max(std::fabs(ymin), std::fabs(ymax));
        if (m <= 0) m = 1.0;
        h->GetYaxis()->SetRangeUser(-m, m);
      } else {
        h->SetMinimum(0.0);
        h->SetMaximum(ymax > 0 ? 1.3 * ymax : 1.0);
      }
      if (logx) {
        h->GetXaxis()->SetMoreLogLabels();
        h->GetXaxis()->SetNoExponent();
      }
      h->Draw(first ? "E1" : "E1 SAME");
      leg->AddEntry(h, (cc == 0) ? "all matched" : Form("#leq %d fake", cc - 1), "lp");
      first = false;
    }
    leg->Draw();
    if (warn) {
      TLatex t;
      t.SetNDC();
      t.SetTextColor(kRed);
      t.SetTextSize(0.035);
      t.DrawLatex(0.13, 0.91, "WARNING: residual range may be < 4#times RMS");
    }
    // On bias plots, show the average Y (over fitted bins) of the 'all matched'
    // curve to judge the overall deviation from 0.
    if (symRange) {
      double sum = 0.0;
      int n = 0;
      for (int b = 1; b <= hs[0]->GetNbinsX(); ++b)
        if (hs[0]->GetBinError(b) > 0) { sum += hs[0]->GetBinContent(b); ++n; }
      const double avg = (n > 0) ? sum / n : 0.0;
      TLatex t;
      t.SetNDC();
      t.SetTextColor(colors[0]);
      t.SetTextSize(0.035);
      t.DrawLatex(0.13, 0.85, Form("#LTall matched#GT = %.3g", avg));
    }
    return c;
  };

  std::vector<TCanvas*> resSigCanv, resMeanCanv;
  for (int p = 0; p < 3; ++p) {
    for (int a = 0; a < 3; ++a) {
      const bool lx = (a == 0 && logPt);
      resSigCanv.push_back(resDraw(resSig[p][a],
          Form("cResSig_%s_%s", resParName[p], resAxName[a]),
          Form("%s resolution vs %s", resParName[p], resAxName[a]),
          resSigTitle[p], lx, resWarn[p], false));
      if (drawBias)
        resMeanCanv.push_back(resDraw(resMean[p][a],
            Form("cResBias_%s_%s", resParName[p], resAxName[a]),
            Form("%s bias vs %s", resParName[p], resAxName[a]),
            resBiasTitle[p], lx, resWarn[p], true));
    }
  }

  // QA residual-distribution canvases (per parameter, per eta bin; pads = pT bins).
  std::vector<TCanvas*> qaCanv;
  if (drawResiduals) {
    const int ncol = 3;
    const int nrow = (nQaPt + ncol - 1) / ncol;
    for (int p = 0; p < 3; ++p) {
      for (int ie = 0; ie < nQaEta; ++ie) {
        TCanvas* c = new TCanvas(
            Form("cQA_%s_e%d", resParName[p], ie),
            Form("QA %s residual, eta [%.2f,%.2f]", resParName[p], qaEta[ie].first,
                 qaEta[ie].second),
            1200, 800);
        c->Divide(ncol, nrow > 0 ? nrow : 1);
        for (int ip = 0; ip < nQaPt; ++ip) {
          c->cd(ip + 1);
          TH1D* h = qaRes[p][ie][ip];
          h->GetXaxis()->SetTitle(resResTitle[p]);
          h->SetTitle(Form("p_{T} [%.2f,%.2f]", qaPt[ip].first, qaPt[ip].second));
          h->Draw("E");
          if (h->GetEntries() >= 20) {
            const double mu0 = h->GetMean(), rms = h->GetRMS();
            if (rms > 0) {
              TF1* g = new TF1(Form("gqa_%s_%d_%d", resParName[p], ie, ip), "gaus",
                               mu0 - 3 * rms, mu0 + 3 * rms);
              h->Fit(g, "QNR");
              g->SetLineColor(kRed);
              g->Draw("SAME");
              TLatex t;
              t.SetNDC();
              t.SetTextSize(0.06);
              t.DrawLatex(0.15, 0.80, Form("#sigma = %.3g", g->GetParameter(2)));
            }
          }
        }
        qaCanv.push_back(c);
      }
    }
  }

  // Tag encoding the pT / eta acceptance and species, appended to the output
  // file names (e.g. "_pt_0.2_1_eta_-2.5_2.5_pid_211"; pid 999 = all species).
  // PID tag: codes joined by '-', e.g. "211-321-2212-11-13"; "999"/empty => "999".
  std::string pidTag;
  if (pidCodes.empty()) {
    pidTag = "999";
  } else {
    for (std::size_t i = 0; i < pidCodes.size(); ++i) {
      if (i) pidTag += "-";
      pidTag += std::to_string(pidCodes[i]);
    }
  }
  std::string rangeTag = Form("_pt_%g_%g_eta_%g_%g_pid_%s", ptMin, ptMax, etaMin, etaMax,
                              pidTag.c_str());

  // QA canvas (3x2): measurement/hit counts, layer-with-measurement/hit counts,
  // and the per-track (measurements - layers) overlap/merge multiplicity.
  TCanvas* cQAcounts =
      new TCanvas("cQAcounts", "QA: measurement & layer counts", 1500, 800);
  cQAcounts->Divide(3, 2);
  {
    TH1D* qh[5] = {hQaNMeas, hQaNLayMeas, hQaMeasMinusLay, hQaNHits, hQaNLayHits};
    const int qc[5] = {kBlue + 1, kAzure + 1, kGreen + 2, kRed + 1, kOrange + 7};
    for (int i = 0; i < 5; ++i) {
      cQAcounts->cd(i + 1);
      gPad->SetGrid();
      qh[i]->SetLineColor(qc[i]);
      qh[i]->SetLineWidth(2);
      qh[i]->SetFillColorAlpha(qc[i], 0.25);
      qh[i]->Draw("hist");
      TLatex tl;
      tl.SetNDC();
      tl.SetTextSize(0.045);
      tl.DrawLatex(0.50, 0.84, Form("mean %.2f", qh[i]->GetMean()));
      tl.DrawLatex(0.50, 0.78, Form("entries %.0f", qh[i]->GetEntries()));
    }
  }

  // Save canvases into one PDF (fake-rate, then efficiency, then per-layer).
  std::string pdf = std::string(dirRec) + "/EffResFakes_plots" + rangeTag + ".pdf";
  std::vector<TCanvas*> pages = {cExPt, cExEt, cExPh, cAlPt, cAlEt, cAlPh};
  if (drawFakeOrigins) {
    pages.push_back(cOrigPt);
    pages.push_back(cOrigEt);
    pages.push_back(cOrigPh);
  }
  if (haveEff) {
    pages.push_back(cEffBreakPt);
    pages.push_back(cEffExPt);
    pages.push_back(cEffExEt);
    pages.push_back(cEffExPh);
    pages.push_back(cEffAmPt);
    pages.push_back(cEffAmEt);
    pages.push_back(cEffAmPh);
  }
  pages.push_back(cLayer);
  if (drawFakeOrigins) {
    pages.push_back(cLayerOrigFrac);
    pages.push_back(cLayerMat);
    pages.push_back(cLayerMatFrac);
    pages.push_back(cLayerRateFrac);
    pages.push_back(cLayerMatRateFrac);
  }
  pages.push_back(cLayerHole);
  if (drawContribMult) pages.push_back(cLayerContrib);
  for (auto* c : resSigCanv) pages.push_back(c);
  for (auto* c : resMeanCanv) pages.push_back(c);
  for (auto* c : qaCanv) pages.push_back(c);
  pages.push_back(cQAcounts);
  for (std::size_t i = 0; i < pages.size(); ++i) {
    std::string suffix = (i == 0) ? "(" : (i + 1 == pages.size() ? ")" : "");
    pages[i]->Print((pdf + suffix).c_str());
  }
  std::cout << "Wrote plots to " << pdf << std::endl;

  // Save all histograms into a ROOT file for later inspection. Insert the
  // pT/eta range tag before the ".root" extension of outFile.
  std::string outName = outFile;
  const std::string rootExt = ".root";
  if (outName.size() >= rootExt.size() &&
      outName.compare(outName.size() - rootExt.size(), rootExt.size(), rootExt) == 0)
    outName = outName.substr(0, outName.size() - rootExt.size());
  std::string outPath = std::string(dirRec) + "/" + outName + rangeTag + rootExt;
  TFile* fOut = TFile::Open(outPath.c_str(), "RECREATE");
  hAllPt->Write();
  hAllEta->Write();
  hAllPhi->Write();
  hQaNMeas->Write();
  hQaNLayMeas->Write();
  hQaMeasMinusLay->Write();
  hQaNHits->Write();
  hQaNLayHits->Write();
  cQAcounts->Write();
  for (int k = 1; k <= maxFake; ++k) {
    rExactPt[k]->Write();
    rExactEta[k]->Write();
    rExactPhi[k]->Write();
    rAtLeastPt[k]->Write();
    rAtLeastEta[k]->Write();
    rAtLeastPhi[k]->Write();
  }
  hLayerRate->Write();
  hLayerHoleRate->Write();
  if (drawContribMult) {
    hLayerContrib->Write();
    cLayerContrib->Write();
  }
  if (drawFakeOrigins) {
    for (int o = 0; o < 4; ++o) {
      rOrigPt[o]->Write();
      rOrigEta[o]->Write();
      rOrigPhi[o]->Write();
      hLayerRateOrig[o]->Write();
      hLayerOrigFrac[o]->Write();
    }
    for (int s = 0; s < kNMatSub; ++s) {
      hLayerMatSub[s]->Write();
      hLayerMatFrac[s]->Write();
    }
    cOrigPt->Write();
    cOrigEt->Write();
    cOrigPh->Write();
    cLayerOrigFrac->Write();
    cLayerMat->Write();
    cLayerMatFrac->Write();
    cLayerRateFrac->Write();
    cLayerMatRateFrac->Write();
  }
  cLayerHole->Write();
  cExPt->Write();
  cExEt->Write();
  cExPh->Write();
  cAlPt->Write();
  cAlEt->Write();
  cAlPh->Write();
  cLayer->Write();
  if (haveEff) {
    hDenPt->Write();
    hDenEta->Write();
    hDenPhi->Write();
    effAllPt->Write();
    effAllEta->Write();
    effAllPhi->Write();
    effExactPt[0]->Write();
    effExactEta[0]->Write();
    effExactPhi[0]->Write();
    effAtMostPt[0]->Write();
    effAtMostEta[0]->Write();
    effAtMostPhi[0]->Write();
    for (int k = 1; k <= maxFake; ++k) {
      effExactPt[k]->Write();
      effExactEta[k]->Write();
      effExactPhi[k]->Write();
      effAtMostPt[k]->Write();
      effAtMostEta[k]->Write();
      effAtMostPhi[k]->Write();
    }
    cEffBreakPt->Write();
    cEffExPt->Write();
    cEffExEt->Write();
    cEffExPh->Write();
    cEffAmPt->Write();
    cEffAmEt->Write();
    cEffAmPh->Write();
  }
  if (haveGen) {
    effPhysAllPt->Write();
    effPhysAllEta->Write();
    effPhysAllPhi->Write();
    for (int k = 0; k <= maxFake; ++k) {
      effPhysExactPt[k]->Write();
      effPhysExactEta[k]->Write();
      effPhysExactPhi[k]->Write();
      effPhysAtMostPt[k]->Write();
      effPhysAtMostEta[k]->Write();
      effPhysAtMostPhi[k]->Write();
    }
  }
  // resolution / bias / QA outputs
  for (int p = 0; p < 3; ++p)
    for (int a = 0; a < 3; ++a)
      for (int c = 0; c < nResCat; ++c) {
        resSig[p][a][c]->Write();
        resMean[p][a][c]->Write();
      }
  for (auto* c : resSigCanv) c->Write();
  for (auto* c : resMeanCanv) c->Write();
  for (auto* c : qaCanv) c->Write();
  fOut->Close();
  std::cout << "Wrote histograms to " << outPath << std::endl;
}
