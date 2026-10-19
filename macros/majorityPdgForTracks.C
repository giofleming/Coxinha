// majorityPdgForTracks.C
//
// PDG code of the majority (truth) particle for EVERY reconstructed track,
// looping over all events and all tracks.
//
// A track stores its truth particle only as an ACTS barcode, split across five
// component branches. particles.root stores the same five components plus the
// PDG (particle_type). We match on all five, within the same event.
//
// Branch schema (from the ACTS Root writers):
//   tracksummary : event_nr (scalar); per-track vectors
//                  majorityParticleId_{vertex_primary,vertex_secondary,
//                                      particle,generation,sub_particle}
//   particles    : event_id (scalar); per-particle vectors
//                  {vertex_primary,vertex_secondary,particle,generation,
//                   sub_particle}, particle_type (= PDG)
//
// The two trees may be in different event orders (files are written in
// completion order), so the particles tree is indexed by event_id and looked
// up per event with GetEntryWithIndex - no assumption that entry i is event i.
//
// Usage:
//   root -l -b -q 'majorityPdgForTracks.C("tracksummary_ambi.root","particles.root")'

#include <TFile.h>
#include <TTree.h>
#include <vector>
#include <iostream>

int majorityPdgForTracks(const char* recFile = "tracksummary_ambi.root",
                const char* simFile = "particles.root") {

  // ---- tracksummary: majority barcode components per track (bind once) ---
  TFile* fr = TFile::Open(recFile);
  TTree* tr = (TTree*)fr->Get("tracksummary");   // adjust name if different

  unsigned int recEvent = 0;
  std::vector<unsigned int> *mVP = 0, *mVS = 0, *mP = 0, *mG = 0, *mS = 0;
  tr->SetBranchAddress("event_nr", &recEvent);
  tr->SetBranchAddress("majorityParticleId_vertex_primary",   &mVP);
  tr->SetBranchAddress("majorityParticleId_vertex_secondary", &mVS);
  tr->SetBranchAddress("majorityParticleId_particle",         &mP);
  tr->SetBranchAddress("majorityParticleId_generation",       &mG);
  tr->SetBranchAddress("majorityParticleId_sub_particle",     &mS);

  // ---- particles: PDG + barcode components per particle (bind once) ------
  TFile* fs = TFile::Open(simFile);
  TTree* ts = (TTree*)fs->Get("particles");

  unsigned int simEvent = 0;
  std::vector<int>* pdg = 0;                                 // particle_type
  std::vector<unsigned int> *pVP = 0, *pVS = 0, *pP = 0, *pG = 0, *pS = 0;
  ts->SetBranchAddress("event_id",         &simEvent);
  ts->SetBranchAddress("particle_type",    &pdg);
  ts->SetBranchAddress("vertex_primary",   &pVP);
  ts->SetBranchAddress("vertex_secondary", &pVS);
  ts->SetBranchAddress("particle",         &pP);
  ts->SetBranchAddress("generation",       &pG);
  ts->SetBranchAddress("sub_particle",     &pS);
  ts->BuildIndex("event_id");   // fast per-event lookup below

  long nTracks = 0, nMatched = 0;

  // ---- loop over all events (rec tree entries) --------------------------
  for (Long64_t ie = 0; ie < tr->GetEntries(); ++ie) {
    tr->GetEntry(ie);

    // load this event's particles once (not per track)
    if (ts->GetEntryWithIndex(recEvent) <= 0) {
      std::cout << "event " << recEvent << ": not found in " << simFile << "\n";
      continue;
    }

    // ---- loop over all tracks in the event ------------------------------
    for (size_t t = 0; t < mVP->size(); ++t) {
      ++nTracks;
      unsigned int vp = (*mVP)[t], vs = (*mVS)[t], pa = (*mP)[t],
                   ge = (*mG)[t], su = (*mS)[t];

      int foundPdg = 0;
      bool ok = false;
      for (size_t j = 0; j < pVP->size(); ++j) {
        if ((*pVP)[j] == vp && (*pVS)[j] == vs && (*pP)[j] == pa &&
            (*pG)[j] == ge && (*pS)[j] == su) {
          foundPdg = (*pdg)[j];
          ok = true;
          break;
        }
      }

      if (ok) {
        ++nMatched;
        std::cout << "ev " << recEvent << " trk " << t
                  << " : PDG = " << foundPdg << "\n";
      } else {
        std::cout << "ev " << recEvent << " trk " << t
                  << " : no truth match (fake / majority not written)\n";
      }
    }
  }

  std::cout << "\n" << nMatched << " / " << nTracks
            << " tracks matched to a truth particle\n";
  return 0;
}
