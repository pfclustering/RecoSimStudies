import FWCore.ParameterSet.Config as cms

deepclusteringdumper = cms.EDAnalyzer("DeepClusteringDumper",

    #genParticleCollection             = cms.InputTag("genParticles",""),
    genParticleCollection             = cms.InputTag("genParticles","", "HLT"),
    caloParticleCollection            = cms.InputTag("mix","MergedCaloTruth"),
    ebRechitCollection                = cms.InputTag("ecalRecHit","EcalRecHitsEB","RECO"),
    eeRechitCollection                = cms.InputTag("ecalRecHit","EcalRecHitsEE","RECO"),
    pfRechitCollection                = cms.InputTag("particleFlowRecHitECAL","","RECO"),
    pfClusterCollection               = cms.InputTag("particleFlowClusterECAL","","RECO"),
    ebSuperClusterCollection          = cms.InputTag("particleFlowSuperClusterECAL","particleFlowSuperClusterECALBarrel","RECO"), 
    eeSuperClusterCollection          = cms.InputTag("particleFlowSuperClusterECAL","particleFlowSuperClusterECALEndcapWithPreshower","RECO"), 
    
    doCompression                     = cms.bool(True),  #do the compression of floats
    nBits                             = cms.int32(23),   #nbits for float compression (<=23)

    saveGenParticles                  = cms.bool(True),  #save genParticles information   
    saveCaloParticles                 = cms.bool(True),  #save caloParticles information
    saveSimhits                       = cms.bool(True),  #save simHits information
    saveRechits                       = cms.bool(True),  #save recHits information
    savePFRechits                     = cms.bool(True),  #save pfRecHits information
    savePFCluster                     = cms.bool(True),  #save pfClusters information
    savePFClusterhits                 = cms.bool(True),  #save pfClustershits information
    saveSuperCluster                  = cms.bool(False),  #save superClusters information
    saveShowerShapes                  = cms.bool(False),  #save saveShowerShapes information
    saveScores                        = cms.bool(True),  #save saveScores information
    genID                             = cms.vint32(22,11), #save only caloParticles with this pdgId 
    #genID                            = cms.vdouble(0),  #save only caloParticles with this pdgId 

)
