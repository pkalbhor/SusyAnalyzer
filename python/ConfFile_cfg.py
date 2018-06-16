import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'root://eoscms.cern.ch//eos/cms/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/80000/F283191C-11C4-E611-973D-00215E2EB74E.root'
    )
)

process.demo = cms.EDAnalyzer('SusyAnalyzer', 
       jettag = cms.untracked.InputTag("slimmedJets"),
       prunedGenParticles = cms.untracked.InputTag("prunedGenParticles"),
       slimmedElectrons = cms.untracked.InputTag("slimmedElectrons"),
       slimmedMuons = cms.untracked.InputTag("slimmedMuons"),
       slimmedPhotons = cms.untracked.InputTag("slimmedPhotons")
)


# Define output file name
import os
process.TFileService = cms.Service("TFileService", fileName = cms.string("Susy_Tree.root"))#(os.getenv('CMSSW_BASE') + '/src/Demo/DemoAnalyzer/test/Pt_ditribution.root'))

process.p = cms.Path(process.demo)
