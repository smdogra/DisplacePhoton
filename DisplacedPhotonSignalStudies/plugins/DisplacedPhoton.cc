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

// system include files
#include <memory>
#include "DisplacedPhoton.h"
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
//
// class declaration
//
// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.
//
// constants, enums and typedefs
//
//
// static data member definitions
//
//
// constructors and destructor
//
DisplacedPhoton::DisplacedPhoton(const edm::ParameterSet& iConfig)
  :
  //  enableGenLLPInfo_(iConfig.getParameter<bool> ("enableGenLLPInfo")),
  genParticlesToken_(consumes<reco::GenParticleCollection>(iConfig.getParameter < edm::InputTag > ("genParticles")))
  //packedGenToken_(consumes<pat::PackedGenParticleCollection>(iConfig.getParameter <edm::InputTag> ("packedGenParticles")))
{  
  //now do what ever initialization is needed
}

    
DisplacedPhoton::~DisplacedPhoton() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  //
  // please remove this method altogether if it would be left empty
}

//
// member functions
//

// ------------ method called for each event  ------------
void DisplacedPhoton::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using std::vector;
  using namespace edm;
  using namespace reco;
  using namespace std;
  std::vector<const reco::Candidate*> prunedV;//Allows easier comparison for mother finding
  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByToken(genParticlesToken_,genParticles);

  tot_events++;
  for ( int i = 0; i < LLP_ARRAY_SIZE; i++ )    {
    gLLP_prod_vertex_x[i] = -666.;
    gLLP_prod_vertex_y[i] = -666.;
    gLLP_prod_vertex_z[i] = -666.;
    gLLP_decay_vertex_x[i] = -666.;
    gLLP_decay_vertex_y[i] = -666.;
    gLLP_decay_vertex_z[i] = -666.;
    gLLP_beta[i] = -666.;
    gLLP_pt[i] = -666.;
    gLLP_eta[i] = -666.;
    gLLP_e[i] = -666.;
    gLLP_phi[i] = -666.;
    gLLP_csc[i] = false;
    gLLP_dt[i] = false;
    gLLP_travel_time[i] = -666.;    

    
    }
  vector<int> llpIDs;
  /*  llpIDs.push_back(9000006);
  llpIDs.push_back(9000007);
  llpIDs.push_back(1023);
  llpIDs.push_back(1000023);
  llpIDs.push_back(1000025);
  llpIDs.push_back(6000113);
  llpIDs.push_back(9900012);
  llpIDs.push_back(9900014);
  llpIDs.push_back(9900016);
  */
  llpIDs.push_back(1000022);
  llpIDs.push_back(1000023);

  for(size_t i=0; i<genParticles->size(); i++){
    bool alreadySaved = false;
    if(
       (abs((*genParticles)[i].pdgId()) >= 1 && abs((*genParticles)[i].pdgId()) <= 6 && ( (*genParticles)[i].status() < 30 ))
       || (abs((*genParticles)[i].pdgId()) >= 11 && abs((*genParticles)[i].pdgId()) <= 16)
       || (abs((*genParticles)[i].pdgId()) == 21 && (*genParticles)[i].status() < 30)
       || (abs((*genParticles)[i].pdgId()) == 22 && (*genParticles)[i].pt() > 10.0 )
       || (abs((*genParticles)[i].pdgId()) >= 23 && abs((*genParticles)[i].pdgId()) <= 25)
       || (abs((*genParticles)[i].pdgId()) >= 32 && abs((*genParticles)[i].pdgId()) <= 42)
       //|| (abs((*genParticles)[i].pdgId()) >= 100 && abs((*genParticles)[i].pdgId()) <= 350)                                                                                                                          
       || (abs((*genParticles)[i].pdgId()) == 1023)
       || (abs((*genParticles)[i].pdgId()) >= 1000001 && abs((*genParticles)[i].pdgId()) <= 1000039)
       || (abs((*genParticles)[i].pdgId()) == 9000006 || abs((*genParticles)[i].pdgId()) == 9000007)
       || (abs((*genParticles)[i].pdgId()) == 6000113 )
       || (abs((*genParticles)[i].pdgId()) == 9900012 || abs((*genParticles)[i].pdgId()) == 9900014 || abs((*genParticles)[i].pdgId()) == 9900016)
       )
      {
        if ((*genParticles)[i].pt()>pt_cut){
	  prunedV.push_back(&(*genParticles)[i]);
	  alreadySaved = true;
	}
      }
    
    //if particle is a daughter of a tau, then save it 
    if (!alreadySaved) {
      if((*genParticles)[i].numberOfMothers() > 0) {
	const reco::Candidate* firstMotherWithDifferentID = findFirstMotherWithDifferentID(&(*genParticles)[i]);
        if (firstMotherWithDifferentID && abs(firstMotherWithDifferentID->pdgId()) == 15 ) {
	  prunedV.push_back(&(*genParticles)[i]);
          alreadySaved = true;
        }
      }
    }
  } //loop over all gen particles
  //Total number of gen particles                                                                                                                   
  nGenParticle = prunedV.size();
  if (nGenParticle > GENPARTICLEARRAYSIZE) {
    cout << "ERROR: nGenParticle exceeded maximum array size: " << GENPARTICLEARRAYSIZE << "\n";
    assert(false);
  }
  //  cout<< "genParticles size " <<  genParticles->size() <<"  " << prunedV.size() << endl;
  bool _found_first_llp = false;
  //Look for mother particle and Fill gen variables                                                                                                 
  for(unsigned int i = 0; i < prunedV.size(); i++) {
    gParticleId[i] = prunedV[i]->pdgId();
    gParticleStatus[i] = prunedV[i]->status();
    gParticleE[i] = prunedV[i]->energy();
    gParticlePt[i] = prunedV[i]->pt();
    gParticlePx[i] = prunedV[i]->px();
    gParticlePy[i] = prunedV[i]->py();
    gParticlePz[i] = prunedV[i]->pz();
    gParticleEta[i] = prunedV[i]->eta();
    gParticlePhi[i] = prunedV[i]->phi();
    gParticleProdVertexX[i] = prunedV[i]->vx();
    gParticleProdVertexY[i] = prunedV[i]->vy();
    gParticleProdVertexZ[i] = prunedV[i]->vz();
    gParticleMotherId[i] = 0;
    gParticleMotherIndex[i] = -1;

    if(prunedV[i]->numberOfMothers() > 0) {
      //find the ID of the first mother that has a different ID than the particle itself                                                      
      const reco::Candidate* firstMotherWithDifferentID = findFirstMotherWithDifferentID(prunedV[i]);
      if (firstMotherWithDifferentID){
	gParticleMotherId[i] = firstMotherWithDifferentID->pdgId();
	//gParticleDecayVertexX[i] = firstMotherWithDifferentID->vx();                                                                        
	//gParticleDecayVertexY[i] = firstMotherWithDifferentID->vy();                                                                        
	//	gParticleDecayVertexZ[i] = firstMotherWithDifferentID->vz();                                                                        
      }
      
      //find the mother and keep going up the mother chain if the ID's are the same                                                                 
      const reco::Candidate* originalMotherWithSameID = findOriginalMotherWithSameID(prunedV[i]);
      for(unsigned int j = 0; j < prunedV.size(); j++) {
        if(prunedV[j] == originalMotherWithSameID) {
          gParticleMotherIndex[i] = j;
          break;
        }
      }
    } else {
      gParticleMotherIndex[i] = -1;
    }
    //---------------------------------------                                                                                                       
    //Find LLPs production and decay vertices                                                                                                       
    //---------------------------------------                                                                                                       
    //    cout<<  " " <<  enableGenLLPInfo_ << endl;
    //    if (enableGenLLPInfo_) {
    {  
    //match with one of the entries in the llpIDs List                                                                                            
      bool matchedLLP = false;
      int matchedLLPID = 0;
      if (gParticleStatus[i] == 22) {
        for (uint d=0 ; d < llpIDs.size() ; d++) {
          if ( abs(gParticleId[i]) == llpIDs[d] ) {
            matchedLLPID = gParticleId[i];
            matchedLLP = true;
          }
        }
      }
      
      if ( matchedLLP ) {
        if (!_found_first_llp) {
          gLLP_prod_vertex_x[0] = prunedV[i]->vx();
          gLLP_prod_vertex_y[0] = prunedV[i]->vy();
          gLLP_prod_vertex_z[0] = prunedV[i]->vz();
        }
        else {
          gLLP_prod_vertex_x[1] = prunedV[i]->vx();
          gLLP_prod_vertex_y[1] = prunedV[i]->vy();
          gLLP_prod_vertex_z[1] = prunedV[i]->vz();
        }
	
        const reco::Candidate *dau = 0;
        bool foundDaughter = false;
        bool noDaughter = false;
        const reco::Candidate *tmpParticle = prunedV[i];
	
        while (!foundDaughter && !noDaughter) {
          if (tmpParticle->numberOfDaughters() > 0) {
            dau = tmpParticle->daughter(0);
            if (dau && (dau->pdgId() != matchedLLPID)) {
              foundDaughter = true;
            } else {
              tmpParticle = dau;
            }
          } else {
            noDaughter = true;
          }
        }
	if (foundDaughter) {
          gParticleDecayVertexX[i] = dau->vx();
          gParticleDecayVertexY[i] = dau->vy();
          gParticleDecayVertexZ[i] = dau->vz();

          if (!_found_first_llp) {
            _found_first_llp = true;
            gLLP_decay_vertex_x[0] = dau->vx();
            gLLP_decay_vertex_y[0] = dau->vy();
            gLLP_decay_vertex_z[0] = dau->vz();
	    //cout<<"Decay vertex " <<  dau->vx() << " " <<  dau->vy() << " " << dau->vz() << endl;;
	    //	    cout<< gLLP_decay_vertex_x[0] <<" " <<  gLLP_decay_vertex_y[0] <<" " <<  gLLP_decay_vertex_z[0]<< endl;
            gLLP_pt[0] = sqrt(gParticlePx[i]*gParticlePx[i]+gParticlePy[i]*gParticlePy[i]);
            gLLP_e[0] = gParticleE[i];
            gLLP_eta[0] = gParticleEta[i];
            gLLP_phi[0] = gParticlePhi[i];
            gLLP_beta[0] = sqrt(gParticlePx[i]*gParticlePx[i]+gParticlePy[i]*gParticlePy[i]+gParticlePz[i]*gParticlePz[i])/gParticleE[i];
            gLLP_travel_time[0] = sqrt(pow(gLLP_decay_vertex_x[0]-gLLP_prod_vertex_x[0],2)
				       +pow(gLLP_decay_vertex_y[0]-gLLP_prod_vertex_y[0],2)
				       +pow(gLLP_decay_vertex_z[0]-gLLP_prod_vertex_z[0],2))/(30. * gLLP_beta[0]);//1/30 is to convert cm to ns                                                                                                           
	    
            double radius = sqrt( pow(gLLP_decay_vertex_x[0],2) + pow(gLLP_decay_vertex_y[0],2) );
            if (abs(gLLP_eta[0]) < 2.4 && abs(gLLP_eta[0]) > 0.9
		&& abs(gLLP_decay_vertex_z[0])<1100 && abs(gLLP_decay_vertex_z[0])>568
		&& radius < 695.5) gLLP_csc[0] = true;
	    if (radius < 740 && radius > 400 && abs(gLLP_decay_vertex_z[0])< 650 ) gLLP_dt[0] = true;
	    cout<<gLLP_pt[0] << " " << gLLP_eta[0] <<  "  radius  " << radius<<  " gLLP_decay_vertex_z    " <<  gLLP_decay_vertex_z[0] <<  endl; 
	  }//loop over all gLLP daughters 
          
	  } //end if first llp                                                        
	  else {
            gLLP_decay_vertex_x[1] = dau->vx();
            gLLP_decay_vertex_y[1] = dau->vy();
            gLLP_decay_vertex_z[1] = dau->vz();
            gLLP_pt[1] = sqrt(gParticlePx[i]*gParticlePx[i]+gParticlePy[i]*gParticlePy[i]);
            gLLP_e[1] = gParticleE[i];
            gLLP_eta[1] = gParticleEta[i];
            gLLP_phi[1] = gParticlePhi[i];
            gLLP_beta[1] = sqrt(gParticlePx[i]*gParticlePx[i]+gParticlePy[i]*gParticlePy[i]+gParticlePz[i]*gParticlePz[i])/gParticleE[i];
            gLLP_travel_time[1] = sqrt(pow(gLLP_decay_vertex_x[1]-gLLP_prod_vertex_x[1],2)
				       +pow(gLLP_decay_vertex_y[1]-gLLP_prod_vertex_y[1],2)
				       +pow(gLLP_decay_vertex_z[1]-gLLP_prod_vertex_z[1],2))/(30. * gLLP_beta[1]);//1/30 is to convert cm to ns  
            double radius = sqrt( pow(gLLP_decay_vertex_x[1],2) + pow(gLLP_decay_vertex_y[1],2) );
            if (abs(gLLP_eta[1]) < 2.4 && abs(gLLP_eta[1]) > 0.9
                && abs(gLLP_decay_vertex_z[1])<1100 && abs(gLLP_decay_vertex_z[1])>568
                && radius < 695.5) gLLP_csc[1] = true;
            if (radius < 740 && radius > 400 && abs(gLLP_decay_vertex_z[1])< 650 ) gLLP_dt[1] = true;
	  } // end of else
       	cout<<"CSC info for First and 2nd Photon "<<  gLLP_csc[0]  << " " << gLLP_csc[1] <<" DT info for Ist and 2nd Photon " << gLLP_dt[0] <<  " " << gLLP_dt[1] << endl;
	if (gLLP_csc[0] == 1 || gLLP_csc[1] == 1)csc1_eventsx++;
	if (gLLP_csc[1] == 0 && gLLP_csc[1] == 0)csc2_events++;	
	if (gLLP_dt[0] == 1)dt1_events++;
	if (gLLP_dt[1] == 1)dt2_events++;
	
      } // end of matchP
    } //end if enableGenLLPInfo

    
 }// for loop of genParticles         
  cout<<"Event number = "<<  tot_events<<"  CSC=  " <<  csc1_events << "  CSC2 =  " << csc2_events <<"  DT1 =  "<< dt1_events <<"  DT2=  "<< dt2_events<< endl;
}
  
// ------------ method called once each job just before starting event loop  ------------
void DisplacedPhoton::beginJob() {
  // please remove this method if not needed
}

// ------------ method called once each job just after ending the event loop  ------------
void DisplacedPhoton::endJob() {
  // please remove this method if not needed
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void DisplacedPhoton::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

const reco::Candidate* DisplacedPhoton::findFirstMotherWithDifferentID(const reco::Candidate *particle){
  if( particle == 0 ){
    printf("ERROR! null candidate pointer, this should never happen\n");
    return 0;
  }
  
  // Is this the first parent with a different ID? If yes, return, otherwise                                                                                                                                            
  // go deeper into recursion                                                                                                                                                                                           
  if (particle->numberOfMothers() > 0 && particle->pdgId() != 0) {
    if (particle->pdgId() == particle->mother(0)->pdgId()
        && particle->mother(0)->status() != 11  // prevent infinite loop for sherpa documentation gluons                                                                                                                
        ) {
      return findFirstMotherWithDifferentID(particle->mother(0));
    } else {
      return particle->mother(0);
    }
  }

  return 0;
};

const reco::Candidate* DisplacedPhoton::findOriginalMotherWithSameID(const reco::Candidate *particle){

  if( particle == 0 ){
    printf("ERROR! null candidate pointer, this should never happen\n");
    return 0;
  }

  // Is there another parent with the same ID? If yes, go deeper into recursion
  if (particle->numberOfMothers() > 0 && particle->pdgId() != 0) {
    if (particle->mother(0)->numberOfMothers() == 0 ||
	particle->mother(0)->status() == 11 ||  // prevent infinite loop for sherpa documentation gluons
	(particle->mother(0)->numberOfMothers() > 0 && particle->mother(0)->mother(0)->pdgId() != particle->mother(0)->pdgId())
	) {
      return particle->mother(0);
    } else {
      return findOriginalMotherWithSameID(particle->mother(0));
    }
  }

  return 0;
}

//define this as a plug-in
DEFINE_FWK_MODULE(DisplacedPhoton);
