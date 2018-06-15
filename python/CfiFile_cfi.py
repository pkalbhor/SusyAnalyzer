import FWCore.ParameterSet.Config as cms

demo = cms.EDAnalyzer('SusyAnalyzer'
     ,tracks = cms.untracked.InputTag('ctfWithMaterialTracks')
)
