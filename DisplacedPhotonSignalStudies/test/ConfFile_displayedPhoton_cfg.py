import FWCore.ParameterSet.Config as cms

process = cms.Process("DisplacedPhoton")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = 10000
process.source = cms.Source("PoolSource",
                                # replace 'myfile.root' with the source file you want to use
                                fileNames = cms.untracked.vstring(
                                    #'file:/afs/cern.ch/user/s/sdogra/work/DarkMatter/DelayedPhoton/TriggerStudies4Run3/re-emulate/CMSSW_13_0_3/src/step_MiniAOD_CSCVLoosePhoton_30GeV.root',
                                    #'file:/afs/cern.ch/user/s/sdogra/work/DarkMatter/DelayedPhoton/TriggerStudies4Run3/re-emulate/CMSSW_13_0_3/src/step1_CSCVLoosePhoton_30GeV.root'
                                    #'file:/eos/cms/store/user/sdogra/DelayedPhoton/TriggerStudies/Signal/output.root'
                                    #'file:/afs/cern.ch/user/s/sdogra/work/DarkMatter/TriggerStudies4Run3/re-emulate/CMSSW_13_2_4/src/TChiGG_mass400_pl1000_1.root'
                                    'file:/eos/user/s/sdogra/TChiGG_mass200_pl5000_1.root'

                )
                            )

#TFileService for output
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('displacedPhoton_TChiGG_mass400_pl1000_ntuple.root'),
    closeFileFast = cms.untracked.bool(True)
)

process.DisplacedPhoton = cms.EDAnalyzer('DisplacedPhoton',
                              genParticles = cms.InputTag("genParticles"),
                              tracks    = cms.untracked.InputTag('generalTracks'),
                              packedGenParticles = cms.InputTag("packedGenParticles"),
                              )

process.p = cms.Path(process.DisplacedPhoton)
