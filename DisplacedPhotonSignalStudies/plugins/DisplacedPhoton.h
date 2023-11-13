// -*- C++ -*-
//
// Package:    DisplacedPhoton/DisplacedPhoton
// Class:      DisplacedPhoton
//
/**\class DisplacedPhoton DisplacedPhoton.cc DisplacedPhoton/DisplacedPhoton/plugins/DisplacedPhoton.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Sunil Manohar Dogra
//         Created:  Fri, 09 Jun 2023 04:44:47 GMT
//
//

#ifndef _DisplacedPhoton_h
#define _DisplacedPhoton_h

#include <memory>
#include <string>
#include <vector>
#include <tuple>
#include <fstream>
#include <TRandom3.h>
//ROOT includes
#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"
// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

//#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
//#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"

#define OBJECTARRAYSIZE 5000
#define CSCRECHITARRAYSIZE 1000000
#define RECHITARRAYSIZE 20000
#define HORECHITARRAYSIZE 2000
#define GENPARTICLEARRAYSIZE 2000
#define MAX_NPV 1000
#define MAX_NPFCAND 5000
#define MAX_NPU 1000
#define MAX_NBX 1000
#define LLP_ARRAY_SIZE 2
#define LLP_DAUGHTER_ARRAY_SIZE 4
#define LLP_GRAND_DAUGHTER_ARRAY_SIZE 4
#define STRIP_DIGI_THRESHOLD 13.3

int nGenParticle;
int gParticleMotherId[GENPARTICLEARRAYSIZE];
int gParticleMotherIndex[GENPARTICLEARRAYSIZE];
int gParticleId[GENPARTICLEARRAYSIZE];
int gParticleStatus[GENPARTICLEARRAYSIZE];
float gParticleE[GENPARTICLEARRAYSIZE];
float gParticlePt[GENPARTICLEARRAYSIZE];
float gParticlePx[GENPARTICLEARRAYSIZE];
float gParticlePy[GENPARTICLEARRAYSIZE];
float gParticlePz[GENPARTICLEARRAYSIZE];
float gParticleEta[GENPARTICLEARRAYSIZE];
float gParticlePhi[GENPARTICLEARRAYSIZE];

float gParticleProdVertexX[GENPARTICLEARRAYSIZE];
float gParticleProdVertexY[GENPARTICLEARRAYSIZE];
float gParticleProdVertexZ[GENPARTICLEARRAYSIZE];

float gParticleDecayVertexX[GENPARTICLEARRAYSIZE];
float gParticleDecayVertexY[GENPARTICLEARRAYSIZE];
float gParticleDecayVertexZ[GENPARTICLEARRAYSIZE];
float gLLP_prod_vertex_x[LLP_ARRAY_SIZE];
float gLLP_prod_vertex_y[LLP_ARRAY_SIZE];
float gLLP_prod_vertex_z[LLP_ARRAY_SIZE];
float gLLP_decay_vertex_x[LLP_ARRAY_SIZE];
float gLLP_decay_vertex_y[LLP_ARRAY_SIZE];
float gLLP_decay_vertex_z[LLP_ARRAY_SIZE];
float gLLP_beta[LLP_ARRAY_SIZE];
float gLLP_travel_time[LLP_ARRAY_SIZE];
float gLLP_pt[LLP_ARRAY_SIZE];
float gLLP_e[LLP_ARRAY_SIZE];
float gLLP_eta[LLP_ARRAY_SIZE];
float gLLP_phi[LLP_ARRAY_SIZE];
bool gLLP_csc[LLP_ARRAY_SIZE];
bool gLLP_dt[LLP_ARRAY_SIZE];
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.

using reco::TrackCollection;

class DisplacedPhoton : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit DisplacedPhoton(const edm::ParameterSet&);
  ~DisplacedPhoton() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  const reco::Candidate* findFirstMotherWithDifferentID(const reco::Candidate *particle);
  const reco::Candidate* findOriginalMotherWithSameID(const reco::Candidate *particle);
private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;
  
  edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken_;
  //  edm::EDGetTokenT<pat::PackedGenParticleCollection> packedGenToken_;
  // ----------member data ---------------------------
  const double pt_cut = 0.0;
  bool    isData_;
  bool    useGen_;
  bool    isRECO_;
  bool    isRAW_;
  bool    isFastsim_;

  bool    isBParkAOD_;
  bool enableTriggerInfo_;
  bool enableGenLLPInfo_;
  bool enableEcalRechits_;
  bool readGenVertexTime_;

  bool enableAK8Jets_;
  bool readMuonDigis_;
  int  llpId_;
  float gen_time[LLP_DAUGHTER_ARRAY_SIZE];
  float gen_time_pv[LLP_DAUGHTER_ARRAY_SIZE];
  float gLLP_daughter_travel_time[LLP_DAUGHTER_ARRAY_SIZE];
  int   gLLP_daughter_id[LLP_DAUGHTER_ARRAY_SIZE];
  float gLLP_daughter_pt[LLP_DAUGHTER_ARRAY_SIZE];
  float gLLP_daughter_eta[LLP_DAUGHTER_ARRAY_SIZE];
  float gLLP_daughter_phi[LLP_DAUGHTER_ARRAY_SIZE];
  float gLLP_daughter_eta_ecalcorr[LLP_DAUGHTER_ARRAY_SIZE];
  float gLLP_daughter_phi_ecalcorr[LLP_DAUGHTER_ARRAY_SIZE];
  float gLLP_daughter_e[LLP_DAUGHTER_ARRAY_SIZE];
  float gLLP_daughter_mass[LLP_DAUGHTER_ARRAY_SIZE];

  //grandaughters
  float gen_time_dau[LLP_GRAND_DAUGHTER_ARRAY_SIZE];
  float gen_time_dau_pv[LLP_GRAND_DAUGHTER_ARRAY_SIZE];
  float photon_travel_time_dau[LLP_DAUGHTER_ARRAY_SIZE];
  float photon_travel_time_dau_pv[LLP_DAUGHTER_ARRAY_SIZE];
  float gLLP_grandaughter_travel_time[LLP_DAUGHTER_ARRAY_SIZE];
  float photon_travel_time[LLP_DAUGHTER_ARRAY_SIZE];
  float photon_travel_time_pv[LLP_DAUGHTER_ARRAY_SIZE];
  int   gLLP_grandaughter_id[LLP_GRAND_DAUGHTER_ARRAY_SIZE];
  float gLLP_grandaughter_pt[LLP_GRAND_DAUGHTER_ARRAY_SIZE];
  float gLLP_grandaughter_eta[LLP_GRAND_DAUGHTER_ARRAY_SIZE];
  float gLLP_grandaughter_phi[LLP_GRAND_DAUGHTER_ARRAY_SIZE];
  float gLLP_grandaughter_eta_ecalcorr[LLP_GRAND_DAUGHTER_ARRAY_SIZE];
  float gLLP_grandaughter_phi_ecalcorr[LLP_GRAND_DAUGHTER_ARRAY_SIZE];
  float gLLP_grandaughter_e[LLP_GRAND_DAUGHTER_ARRAY_SIZE];
  float gLLP_grandaughter_mass[LLP_GRAND_DAUGHTER_ARRAY_SIZE];

  int tot_events=0;
  int csc1_events=0;
  int dt1_events=0;

  int csc2_events=0;
  int dt2_events=0;
};
#endif
