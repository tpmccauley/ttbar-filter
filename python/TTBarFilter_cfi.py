import FWCore.ParameterSet.Config as cms

TTBarFilter = cms.EDFilter('TTBarFilter',
                           muonInputTag = cms.InputTag("muons"),
                           pfJetInputTag = cms.InputTag("ak5PFJets"),
                           csvFileName = cms.string("ttbar.csv"),
                           minMuonPt = cms.double(20.0),
                           maxMuonEta = cms.double(2.1)
                           )
