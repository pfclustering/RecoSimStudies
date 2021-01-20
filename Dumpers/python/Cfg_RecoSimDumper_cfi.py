import FWCore.ParameterSet.Config as cms

import Geometry.CaloEventSetup.caloTowerConstituents_cfi 
CaloTowerConstituentsMapBuilder = cms.ESProducer("CaloTowerConstituentsMapBuilder",
   MapFile = cms.untracked.string('Geometry/CaloTopology/data/CaloTowerEEGeometric.map.gz')
   )

recosimdumper = cms.EDAnalyzer("RecoSimDumper",

    rhoCollection                   = cms.InputTag("fixedGridRhoAll"),
    vertexCollection                = cms.InputTag("offlinePrimaryVertices"),
    genParticleCollection             = cms.InputTag("genParticles",""),
    caloParticleCollection            = cms.InputTag("mix","MergedCaloTruth"),
    ebRechitCollection                = cms.InputTag("ecalRecHit","EcalRecHitsEB","RECO"),
    eeRechitCollection                = cms.InputTag("ecalRecHit","EcalRecHitsEE","RECO"),
    pfRechitCollection                = cms.InputTag("particleFlowRecHitECAL","","RECO"),
    pfClusterCollection               = cms.InputTag("particleFlowClusterECAL","","RECO"),
    ebSuperClusterCollection          = cms.InputTag("particleFlowSuperClusterECAL","particleFlowSuperClusterECALBarrel","RECO"), 
    eeSuperClusterCollection          = cms.InputTag("particleFlowSuperClusterECAL","particleFlowSuperClusterECALEndcapWithPreshower","RECO"), 
    useHcalTowers                   = cms.bool(False),  #compute HoE
    hcalTowersCollection            = cms.InputTag("towerMaker"), #compute HoE
    useRetunedSC                    = cms.bool(False),  #run on new RetunedSCs
    useDeepSC                       = cms.bool(False),  #run on new DeepSCs
    #ebRetunedSuperClusterCollection = cms.InputTag("particleFlowSuperClusterECALNewParams","particleFlowSuperClusterECALBarrelMustacheNewParams","RECO"), 
    #eeRetunedSuperClusterCollection = cms.InputTag("particleFlowSuperClusterECALNewParams","particleFlowSuperClusterECALEndcapWithPreshowerMustacheNewParams","RECO"),
    #ebDeepSuperClusterCollection    = cms.InputTag("particleFlowDeepSuperClusterECAL","particleFlowDeepSuperClusterECALBarrel","RECO"), 
    #eeDeepSuperClusterCollection    = cms.InputTag("particleFlowDeepSuperClusterECAL","particleFlowDeepSuperClusterECALEndcapWithPreshower","RECO"),
    #ebDeepSuperClusterLWPCollection = cms.InputTag("particleFlowDeepSuperClusterLWPECAL","particleFlowDeepSuperClusterLWPECALBarrel","RECO"), 
    #eeDeepSuperClusterLWPCollection = cms.InputTag("particleFlowDeepSuperClusterLWPECAL","particleFlowDeepSuperClusterLWPECALEndcapWithPreshower","RECO"),
    #ebDeepSuperClusterTWPCollection = cms.InputTag("particleFlowDeepSuperClusterTWPECAL","particleFlowDeepSuperClusterTWPECALBarrel","RECO"), 
    #eeDeepSuperClusterTWPCollection = cms.InputTag("particleFlowDeepSuperClusterTWPECAL","particleFlowDeepSuperClusterTWPECALEndcapWithPreshower","RECO"),    
    
    doCompression                     = cms.bool(True),  #do the compression of floats
    nBits                             = cms.int32(23),   #nbits for float compression (<=23)

    saveGenParticles                  = cms.bool(True),  #save genParticles information   
    saveCaloParticles                 = cms.bool(True),  #save caloParticles information
    saveSimhits                       = cms.bool(True),  #save simHits information
    saveRechits                       = cms.bool(True),  #save recHits information
    savePFRechits                     = cms.bool(True),  #save pfRecHits information
    saveEBPFRechits                   = cms.bool(False),  
    savePFCluster                     = cms.bool(True),  #save pfClusters information
    savePFClusterhits                 = cms.bool(True),  #save pfClustershits information
    saveSuperCluster                  = cms.bool(False),  #save superClusters information
    saveShowerShapes                  = cms.bool(False),  #save showerShapes information
    saveScores                        = cms.bool(True),  #save scores information
    scoreType                         = cms.string("sim_fraction"),  #score to be used for caloParticle matching
    genID                             = cms.vint32(22,11, -11,), #save only caloParticles with this pdgId 
    #genID                            = cms.vdouble(0),  #save only caloParticles with this pdgId 

)
