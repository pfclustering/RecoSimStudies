import FWCore.ParameterSet.Config as cms

import Geometry.CaloEventSetup.caloTowerConstituents_cfi 
CaloTowerConstituentsMapBuilder = cms.ESProducer("CaloTowerConstituentsMapBuilder",
   MapFile = cms.untracked.string('Geometry/CaloTopology/data/CaloTowerEEGeometric.map.gz')
   )

recosimdumper = cms.EDAnalyzer("RecoSimDumperPF",

    rhoCollection                   = cms.InputTag("fixedGridRhoAll"),
    vertexCollection                = cms.InputTag("offlinePrimaryVertices"),
    genParticleCollection             = cms.InputTag("genParticles",""),
    caloParticleCollection            = cms.InputTag("mix","MergedCaloTruth"),
    ebRechitCollection                = cms.InputTag("ecalRecHit","EcalRecHitsEB","RECO"),
    eeRechitCollection                = cms.InputTag("ecalRecHit","EcalRecHitsEE","RECO"),
    pfRechitCollection                = cms.InputTag("particleFlowRecHitECAL","","RECO"),
    pfClusterCollection               = cms.InputTag("particleFlowClusterECAL","","RECO"),
    
    doCompression                     = cms.bool(True),  #do the compression of floats
    nBits                             = cms.int32(23),   #nbits for float compression (<=23)

    saveGenParticles                  = cms.bool(True),  #save genParticles information   
    saveCaloParticles                 = cms.bool(True),  #save caloParticles information
    saveSimhits                       = cms.bool(True),  #save simHits information
    savePFRechits                     = cms.bool(True),  #save pfRecHits information
    savePFCluster                     = cms.bool(True),  #save pfClusters information
    saveShowerShapes                  = cms.bool(False),  #save showerShapes information
    saveScores                        = cms.bool(False),  #save scores information
    scoreType                         = cms.string("sim_fraction"),  #score to be used for caloParticle matching
    genID                             = cms.vint32(22,11, -11,), #save only caloParticles with this pdgId 
)
