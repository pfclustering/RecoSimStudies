#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Utils/interface/TFileDirectory.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDConsumerBase.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/LooperFactory.h"
#include "FWCore/Framework/interface/ESProducerLooper.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/ESProducts.h"
#include "FWCore/Framework/interface/Event.h"

#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/FWLite/interface/Handle.h"

#include "FWCore/FWLite/interface/FWLiteEnabler.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
//#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"

#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
#include "Geometry/CaloTopology/interface/CaloSubdetectorTopology.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "Geometry/EcalAlgo/interface/EcalBarrelGeometry.h"
#include "Geometry/EcalAlgo/interface/EcalEndcapGeometry.h"

#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "SimDataFormats/CaloHit/interface/PCaloHitContainer.h"
#include "SimDataFormats/CaloAnalysis/interface/SimCluster.h"
#include "SimDataFormats/CaloAnalysis/interface/CaloParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHit.h"
#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/Math/interface/libminifloat.h"

#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/EgammaCandidates/interface/PhotonCore.h"
#include "DataFormats/EgammaReco/interface/ElectronSeed.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"
#include "FWCore/Utilities/interface/isFinite.h"
#include "CondFormats/EcalObjects/interface/EcalIntercalibConstants.h"
#include "CondFormats/DataRecord/interface/EcalIntercalibConstantsRcd.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "RecoSimStudies/Dumpers/plugins/RecoSimDumperPF.h"
#include "RecoEgamma/EgammaIsolationAlgos/interface/EgammaHcalIsolation.h"
#include "RecoEgamma/EgammaIsolationAlgos/interface/EgammaTowerIsolation.h"
#include "RecoEgamma/EgammaIsolationAlgos/interface/EgammaHadTower.h"
#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"

#include "TSystem.h"
#include "TFile.h"
#include "TProfile.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TProfile2D.h"
#include "TChain.h"
#include "TLorentzVector.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TTree.h"

#include <memory>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <iostream>
#include <map>
#include <vector>
#include <algorithm>
#include <functional>
#include <set>
#include <assert.h>
#include <time.h>

#include <TMath.h>
#include <Math/VectorUtil.h>
//#include <boost/tokenizer.hpp>

using namespace cms;
using namespace edm;
using namespace std;
using namespace reco;

//
// constructors and destructor
//
RecoSimDumperPF::RecoSimDumperPF(const edm::ParameterSet& iConfig)
{

   vtxToken_                      = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertexCollection"));
   rhoToken_                      = consumes<double>(iConfig.getParameter<edm::InputTag>("rhoCollection"));
   genToken_                      = consumes<std::vector<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("genParticleCollection"));
   caloPartToken_                 = consumes<std::vector<CaloParticle> >(iConfig.getParameter<edm::InputTag>("caloParticleCollection"));
   ebRechitToken_                 = consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("ebRechitCollection"));
   eeRechitToken_                 = consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("eeRechitCollection"));
   pfRecHitToken_                 = consumes<std::vector<reco::PFRecHit> >(iConfig.getParameter<edm::InputTag>("pfRechitCollection")); 
   pfClusterToken_                = consumes<std::vector<reco::PFCluster> >(iConfig.getParameter<edm::InputTag>("pfClusterCollection")); 
   ebSuperClusterToken_           = consumes<std::vector<reco::SuperCluster> >(iConfig.getParameter<edm::InputTag>("ebSuperClusterCollection"));
   eeSuperClusterToken_           = consumes<std::vector<reco::SuperCluster> >(iConfig.getParameter<edm::InputTag>("eeSuperClusterCollection"));
   useHcalTowers_                 = iConfig.getParameter<bool>("useHcalTowers");  
   hcalTowersToken_               = consumes<CaloTowerCollection>(iConfig.getParameter<edm::InputTag>("hcalTowersCollection"));
   doCompression_                 = iConfig.getParameter<bool>("doCompression");
   nBits_                         = iConfig.getParameter<int>("nBits");
   saveGenParticles_              = iConfig.getParameter<bool>("saveGenParticles");
   saveCaloParticles_             = iConfig.getParameter<bool>("saveCaloParticles");
   saveSimhits_             	  = iConfig.getParameter<bool>("saveSimhits");
   saveEBPFRechits_               = iConfig.getParameter<bool>("saveEBPFRechits"); 
   savePFCluster_                 = iConfig.getParameter<bool>("savePFCluster");
   savePFClusterhits_             = iConfig.getParameter<bool>("savePFClusterhits");
   saveSuperCluster_              = iConfig.getParameter<bool>("saveSuperCluster");
   saveShowerShapes_              = iConfig.getParameter<bool>("saveShowerShapes");
   scoreType_                     = iConfig.getParameter<std::string>("scoreType");
   saveScores_                    = iConfig.getParameter<bool>("saveScores");
   genID_                         = iConfig.getParameter<std::vector<int>>("genID");
   
   if(nBits_>23 && doCompression_){
      cout << "WARNING: float compression bits > 23 ---> Using 23 (i.e. no compression) instead!" << endl;
      nBits_=23;
   }

   //output file, historgrams and trees
   tree = iFile->make<TTree>("caloTree","caloTree"); 
   tree->Branch("eventId", &eventId, "eventId/L");
   tree->Branch("lumiId", &lumiId, "lumiId/I");
   tree->Branch("runId", &runId, "runId/I");
   tree->Branch("nVtx", &nVtx, "nVtx/I");
   tree->Branch("rho", &rho, "rho/F"); 
   if(saveGenParticles_){
      tree->Branch("genParticle_id","std::vector<int>",&genParticle_id);
      tree->Branch("genParticle_isGammaFromMeson","std::vector<bool>",&genParticle_isGammaFromMeson);
      tree->Branch("genParticle_energy","std::vector<float>",&genParticle_energy);
      tree->Branch("genParticle_pt","std::vector<float>",&genParticle_pt);
      tree->Branch("genParticle_eta","std::vector<float>",&genParticle_eta);
      tree->Branch("genParticle_phi","std::vector<float>",&genParticle_phi);
      if(savePFCluster_) tree->Branch("genParticle_pfCluster_dR_genScore_MatchedIndex","std::vector<std::vector<int> >",&genParticle_pfCluster_dR_genScore_MatchedIndex);
   }
   if(saveCaloParticles_){
      tree->Branch("caloParticle_id","std::vector<int>",&caloParticle_id); 
      tree->Branch("caloParticle_genEnergy","std::vector<float>",&caloParticle_genEnergy);
      tree->Branch("caloParticle_simEnergy","std::vector<float>",&caloParticle_simEnergy); 
      tree->Branch("caloParticle_genPt","std::vector<float>",&caloParticle_genPt);
      tree->Branch("caloParticle_simPt","std::vector<float>",&caloParticle_simPt);
      tree->Branch("caloParticle_genEta","std::vector<float>",&caloParticle_genEta);
      tree->Branch("caloParticle_simEta","std::vector<float>",&caloParticle_simEta);
      tree->Branch("caloParticle_genPhi","std::vector<float>",&caloParticle_genPhi);
      tree->Branch("caloParticle_simPhi","std::vector<float>",&caloParticle_simPhi);
      tree->Branch("caloParticle_simIeta","std::vector<int>",&caloParticle_simIeta);
      tree->Branch("caloParticle_simIphi","std::vector<int>",&caloParticle_simIphi);
      tree->Branch("caloParticle_simIz","std::vector<int>",&caloParticle_simIz);
      if(savePFCluster_){   
         tree->Branch("caloParticle_pfCluster_dR_simScore_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_pfCluster_dR_simScore_MatchedIndex);
         tree->Branch("caloParticle_pfCluster_sim_fraction_old_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_pfCluster_sim_fraction_old_MatchedIndex);
         if(!saveScores_) tree->Branch("caloParticle_pfCluster_simScore_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_pfCluster_simScore_MatchedIndex);
         if(saveScores_){
            tree->Branch("caloParticle_pfCluster_n_shared_xtals_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_pfCluster_n_shared_xtals_MatchedIndex);
            tree->Branch("caloParticle_pfCluster_sim_fraction_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_pfCluster_sim_fraction_MatchedIndex);
            tree->Branch("caloParticle_pfCluster_sim_fraction_1MeVCut_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_pfCluster_sim_fraction_1MeVCut_MatchedIndex);
            tree->Branch("caloParticle_pfCluster_sim_fraction_5MeVCut_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_pfCluster_sim_fraction_5MeVCut_MatchedIndex);
            tree->Branch("caloParticle_pfCluster_sim_fraction_10MeVCut_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_pfCluster_sim_fraction_10MeVCut_MatchedIndex);
            tree->Branch("caloParticle_pfCluster_sim_fraction_50MeVCut_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_pfCluster_sim_fraction_50MeVCut_MatchedIndex);
            tree->Branch("caloParticle_pfCluster_sim_fraction_100MeVCut_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_pfCluster_sim_fraction_100MeVCut_MatchedIndex);
            tree->Branch("caloParticle_pfCluster_sim_fraction_500MeVCut_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_pfCluster_sim_fraction_500MeVCut_MatchedIndex);
            tree->Branch("caloParticle_pfCluster_sim_fraction_1GeVCut_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_pfCluster_sim_fraction_1GeVCut_MatchedIndex);
            tree->Branch("caloParticle_pfCluster_sim_rechit_diff_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_pfCluster_sim_rechit_diff_MatchedIndex);
            tree->Branch("caloParticle_pfCluster_sim_rechit_fraction_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_pfCluster_sim_rechit_fraction_MatchedIndex);   
            tree->Branch("caloParticle_pfCluster_global_sim_rechit_fraction_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_pfCluster_global_sim_rechit_fraction_MatchedIndex); 
            tree->Branch("caloParticle_pfCluster_hgcal_caloToCluster_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_pfCluster_hgcal_caloToCluster_MatchedIndex);  
            tree->Branch("caloParticle_pfCluster_hgcal_clusterToCalo_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_pfCluster_hgcal_clusterToCalo_MatchedIndex); 
            tree->Branch("caloParticle_pfCluster_sim_rechit_combined_fraction_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_pfCluster_sim_rechit_combined_fraction_MatchedIndex);  
            tree->Branch("caloParticle_pfCluster_rechit_sim_combined_fraction_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_pfCluster_rechit_sim_combined_fraction_MatchedIndex);        
         }   
      }
   }
   if(saveSimhits_){
      tree->Branch("simHit_energy","std::vector<float>",&simHit_energy);
      tree->Branch("simHit_eta","std::vector<float>",&simHit_eta);
      tree->Branch("simHit_phi","std::vector<float>",&simHit_phi);
      tree->Branch("simHit_ieta","std::vector<int>",&simHit_ieta);
      tree->Branch("simHit_iphi","std::vector<int>",&simHit_iphi);
      tree->Branch("simHit_iz","std::vector<int>",&simHit_iz);
      tree->Branch("simHit_icP","std::vector<int>",&simHit_icP);
   }
   if(saveEBPFRechits_){ 
      tree->Branch("pfRecHit_energy","std::vector<float>",&pfRecHit_energy);
      tree->Branch("pfRecHit_eta","std::vector<float>",&pfRecHit_eta); 
      tree->Branch("pfRecHit_phi","std::vector<float>",&pfRecHit_phi);
      tree->Branch("pfRecHit_ieta","std::vector<int>",&pfRecHit_ieta); 
      tree->Branch("pfRecHit_iphi","std::vector<int>",&pfRecHit_iphi);
      tree->Branch("pfRecHit_iz","std::vector<int>",&pfRecHit_iz);     
   }
   if(savePFCluster_){
      tree->Branch("pfCluster_energy","std::vector<float>",&pfCluster_energy);
      tree->Branch("pfCluster_eta","std::vector<float>",&pfCluster_eta);
      tree->Branch("pfCluster_phi","std::vector<float>",&pfCluster_phi);   
      tree->Branch("pfCluster_ieta","std::vector<int>",&pfCluster_ieta);
      tree->Branch("pfCluster_iphi","std::vector<int>",&pfCluster_iphi);   
      tree->Branch("pfCluster_iz","std::vector<int>",&pfCluster_iz);
      tree->Branch("pfCluster_nXtals","std::vector<int>",&pfCluster_nXtals);  
      if(saveSuperCluster_) tree->Branch("pfCluster_superClustersIndex","std::vector<std::vector<int> >",&pfCluster_superClustersIndex); 
      if(saveCaloParticles_){ 
         tree->Branch("pfCluster_dR_genScore_MatchedIndex","std::vector<int>",&pfCluster_dR_genScore_MatchedIndex);
         tree->Branch("pfCluster_dR_simScore_MatchedIndex","std::vector<int>",&pfCluster_dR_simScore_MatchedIndex);
         tree->Branch("pfCluster_sim_fraction_old_MatchedIndex","std::vector<int>",&pfCluster_sim_fraction_old_MatchedIndex);
         if(!saveScores_) tree->Branch("pfCluster_simScore_MatchedIndex","std::vector<int>",&pfCluster_simScore_MatchedIndex);
         if(saveScores_){
            tree->Branch("pfCluster_n_shared_xtals_MatchedIndex","std::vector<int>",&pfCluster_n_shared_xtals_MatchedIndex);
            tree->Branch("pfCluster_sim_fraction_MatchedIndex","std::vector<int>",&pfCluster_sim_fraction_MatchedIndex);
            tree->Branch("pfCluster_sim_fraction_1MeVCut_MatchedIndex","std::vector<int>",&pfCluster_sim_fraction_1MeVCut_MatchedIndex); 
            tree->Branch("pfCluster_sim_fraction_5MeVCut_MatchedIndex","std::vector<int>",&pfCluster_sim_fraction_5MeVCut_MatchedIndex);
            tree->Branch("pfCluster_sim_fraction_10MeVCut_MatchedIndex","std::vector<int>",&pfCluster_sim_fraction_10MeVCut_MatchedIndex);
            tree->Branch("pfCluster_sim_fraction_50MeVCut_MatchedIndex","std::vector<int>",&pfCluster_sim_fraction_50MeVCut_MatchedIndex);
            tree->Branch("pfCluster_sim_fraction_100MeVCut_MatchedIndex","std::vector<int>",&pfCluster_sim_fraction_100MeVCut_MatchedIndex);
            tree->Branch("pfCluster_sim_fraction_500MeVCut_MatchedIndex","std::vector<int>",&pfCluster_sim_fraction_500MeVCut_MatchedIndex);
            tree->Branch("pfCluster_sim_fraction_1GeVCut_MatchedIndex","std::vector<int>",&pfCluster_sim_fraction_1GeVCut_MatchedIndex);
            tree->Branch("pfCluster_sim_rechit_diff_MatchedIndex","std::vector<int>",&pfCluster_sim_rechit_diff_MatchedIndex);
            tree->Branch("pfCluster_sim_rechit_fraction_MatchedIndex","std::vector<int>",&pfCluster_sim_rechit_fraction_MatchedIndex);   
            tree->Branch("pfCluster_global_sim_rechit_fraction_MatchedIndex","std::vector<int>",&pfCluster_global_sim_rechit_fraction_MatchedIndex); 
            tree->Branch("pfCluster_hgcal_caloToCluster_MatchedIndex","std::vector<int>",&pfCluster_hgcal_caloToCluster_MatchedIndex); 
            tree->Branch("pfCluster_hgcal_clusterToCalo_MatchedIndex","std::vector<int>",&pfCluster_hgcal_clusterToCalo_MatchedIndex);   
            tree->Branch("pfCluster_sim_rechit_combined_fraction_MatchedIndex","std::vector<int>",&pfCluster_sim_rechit_combined_fraction_MatchedIndex); 
            tree->Branch("pfCluster_rechit_sim_combined_fraction_MatchedIndex","std::vector<int>",&pfCluster_rechit_sim_combined_fraction_MatchedIndex);     
         } 
      } 
      if(saveCaloParticles_){
         tree->Branch("pfCluster_dR_genScore","std::vector<std::vector<double> >",&pfCluster_dR_genScore);
         tree->Branch("pfCluster_dR_simScore","std::vector<std::vector<double> >",&pfCluster_dR_simScore);
         tree->Branch("pfCluster_sim_fraction_old","std::vector<std::vector<double> >",&pfCluster_sim_fraction_old);
         if(!saveScores_) tree->Branch("pfCluster_simScore","std::vector<std::vector<double> >",&pfCluster_simScore);
         if(saveScores_){
            tree->Branch("pfCluster_n_shared_xtals","std::vector<std::vector<double> >",&pfCluster_n_shared_xtals);
            tree->Branch("pfCluster_sim_fraction","std::vector<std::vector<double> >",&pfCluster_sim_fraction);
            tree->Branch("pfCluster_sim_fraction_1MeVCut","std::vector<std::vector<double> >",&pfCluster_sim_fraction_1MeVCut);
            tree->Branch("pfCluster_sim_fraction_5MeVCut","std::vector<std::vector<double> >",&pfCluster_sim_fraction_5MeVCut);
            tree->Branch("pfCluster_sim_fraction_10MeVCut","std::vector<std::vector<double> >",&pfCluster_sim_fraction_10MeVCut);
            tree->Branch("pfCluster_sim_fraction_50MeVCut","std::vector<std::vector<double> >",&pfCluster_sim_fraction_50MeVCut);
            tree->Branch("pfCluster_sim_fraction_100MeVCut","std::vector<std::vector<double> >",&pfCluster_sim_fraction_100MeVCut); 
            tree->Branch("pfCluster_sim_fraction_500MeVCut","std::vector<std::vector<double> >",&pfCluster_sim_fraction_500MeVCut); 
            tree->Branch("pfCluster_sim_fraction_1GeVCut","std::vector<std::vector<double> >",&pfCluster_sim_fraction_1GeVCut);      
            tree->Branch("pfCluster_sim_rechit_diff","std::vector<std::vector<double> >",&pfCluster_sim_rechit_diff);
            tree->Branch("pfCluster_sim_rechit_fraction","std::vector<std::vector<double> >",&pfCluster_sim_rechit_fraction);   
            tree->Branch("pfCluster_global_sim_rechit_fraction","std::vector<std::vector<double> >",&pfCluster_global_sim_rechit_fraction); 
            tree->Branch("pfCluster_hgcal_caloToCluster","std::vector<std::vector<double> >",&pfCluster_hgcal_caloToCluster); 
            tree->Branch("pfCluster_hgcal_clusterToCalo","std::vector<std::vector<double> >",&pfCluster_hgcal_clusterToCalo); 
            tree->Branch("pfCluster_sim_rechit_combined_fraction","std::vector<std::vector<double> >",&pfCluster_sim_rechit_combined_fraction); 
            tree->Branch("pfCluster_rechit_sim_combined_fraction","std::vector<std::vector<double> >",&pfCluster_rechit_sim_combined_fraction);   
         }
      }
      if(savePFClusterhits_){ 
         tree->Branch("pfClusterHit_energy","std::vector<std::vector<float> >",&pfClusterHit_energy);
         tree->Branch("pfClusterHit_rechitEnergy","std::vector<std::vector<float> >",&pfClusterHit_rechitEnergy);
         tree->Branch("pfClusterHit_eta","std::vector<std::vector<float> >",&pfClusterHit_eta);
         tree->Branch("pfClusterHit_phi","std::vector<std::vector<float> >",&pfClusterHit_phi);
         tree->Branch("pfClusterHit_ieta","std::vector<std::vector<int> >",&pfClusterHit_ieta);
         tree->Branch("pfClusterHit_iphi","std::vector<std::vector<int> >",&pfClusterHit_iphi);
         tree->Branch("pfClusterHit_iz","std::vector<std::vector<int> >",&pfClusterHit_iz);
      }   
   }
   if(saveSuperCluster_){
      tree->Branch("superCluster_energy","std::vector<float> ",&superCluster_energy);
      tree->Branch("superCluster_eta","std::vector<float>",&superCluster_eta);
      tree->Branch("superCluster_phi","std::vector<float>",&superCluster_phi);  
      tree->Branch("superCluster_etaWidth","std::vector<float>",&superCluster_etaWidth);
      tree->Branch("superCluster_phiWidth","std::vector<float>",&superCluster_phiWidth);  
      tree->Branch("superCluster_R","std::vector<float>",&superCluster_R);   
      tree->Branch("superCluster_nPFClusters","std::vector<int>",&superCluster_nPFClusters);   
      tree->Branch("superCluster_ieta","std::vector<int>",&superCluster_ieta);
      tree->Branch("superCluster_iphi","std::vector<int>",&superCluster_iphi);  
      tree->Branch("superCluster_iz","std::vector<int>",&superCluster_iz);    
      if(savePFCluster_) tree->Branch("superCluster_seedIndex","std::vector<int>",&superCluster_seedIndex);     
      if(savePFCluster_) tree->Branch("superCluster_pfClustersIndex","std::vector<std::vector<int> >",&superCluster_pfClustersIndex); 
      tree->Branch("superCluster_psCluster_energy", "std::vector<std::vector<float> >", &superCluster_psCluster_energy);
      tree->Branch("superCluster_psCluster_eta", "std::vector<std::vector<float> >", &superCluster_psCluster_eta);
      tree->Branch("superCluster_psCluster_phi", "std::vector<std::vector<float> >", &superCluster_psCluster_phi); 
   }
   if(savePFCluster_ && saveShowerShapes_){  
      tree->Branch("pfCluster_e5x5","std::vector<float>",&pfCluster_e5x5);
      tree->Branch("pfCluster_e2x2Ratio","std::vector<float>",&pfCluster_e2x2Ratio);
      tree->Branch("pfCluster_e3x3Ratio","std::vector<float>",&pfCluster_e3x3Ratio);
      tree->Branch("pfCluster_eMaxRatio","std::vector<float>",&pfCluster_eMaxRatio);
      tree->Branch("pfCluster_e2ndRatio","std::vector<float>",&pfCluster_e2ndRatio);
      tree->Branch("pfCluster_eTopRatio","std::vector<float>",&pfCluster_eTopRatio);
      tree->Branch("pfCluster_eRightRatio","std::vector<float>",&pfCluster_eRightRatio);
      tree->Branch("pfCluster_eBottomRatio","std::vector<float>",&pfCluster_eBottomRatio);
      tree->Branch("pfCluster_eLeftRatio","std::vector<float>",&pfCluster_eLeftRatio);
      tree->Branch("pfCluster_e2x5MaxRatio","std::vector<float>",&pfCluster_e2x5MaxRatio);
      tree->Branch("pfCluster_e2x5TopRatio","std::vector<float>",&pfCluster_e2x5TopRatio);
      tree->Branch("pfCluster_e2x5RightRatio","std::vector<float>",&pfCluster_e2x5RightRatio);
      tree->Branch("pfCluster_e2x5BottomRatio","std::vector<float>",&pfCluster_e2x5BottomRatio); 
      tree->Branch("pfCluster_e2x5LeftRatio","std::vector<float>",&pfCluster_e2x5LeftRatio); 
      tree->Branch("pfCluster_swissCross","std::vector<float>",&pfCluster_swissCross); 
      tree->Branch("pfCluster_r9","std::vector<float>",&pfCluster_r9);
      tree->Branch("pfCluster_sigmaIetaIeta","std::vector<float>",&pfCluster_sigmaIetaIeta);
      tree->Branch("pfCluster_sigmaIetaIphi","std::vector<float>",&pfCluster_sigmaIetaIphi);
      tree->Branch("pfCluster_sigmaIphiIphi","std::vector<float>",&pfCluster_sigmaIphiIphi);
      tree->Branch("pfCluster_full5x5_e5x5","std::vector<float>",&pfCluster_full5x5_e5x5);
      tree->Branch("pfCluster_full5x5_e2x2Ratio","std::vector<float>",&pfCluster_full5x5_e2x2Ratio);
      tree->Branch("pfCluster_full5x5_e3x3Ratio","std::vector<float>",&pfCluster_full5x5_e3x3Ratio);
      tree->Branch("pfCluster_full5x5_eMaxRatio","std::vector<float>",&pfCluster_full5x5_eMaxRatio);
      tree->Branch("pfCluster_full5x5_e2ndRatio","std::vector<float>",&pfCluster_full5x5_e2ndRatio);
      tree->Branch("pfCluster_full5x5_eTopRatio","std::vector<float>",&pfCluster_full5x5_eTopRatio);
      tree->Branch("pfCluster_full5x5_eRightRatio","std::vector<float>",&pfCluster_full5x5_eRightRatio);
      tree->Branch("pfCluster_full5x5_eBottomRatio","std::vector<float>",&pfCluster_full5x5_eBottomRatio);
      tree->Branch("pfCluster_full5x5_eLeftRatio","std::vector<float>",&pfCluster_full5x5_eLeftRatio);
      tree->Branch("pfCluster_full5x5_e2x5MaxRatio","std::vector<float>",&pfCluster_full5x5_e2x5MaxRatio);
      tree->Branch("pfCluster_full5x5_e2x5TopRatio","std::vector<float>",&pfCluster_full5x5_e2x5TopRatio);
      tree->Branch("pfCluster_full5x5_e2x5RightRatio","std::vector<float>",&pfCluster_full5x5_e2x5RightRatio);
      tree->Branch("pfCluster_full5x5_e2x5BottomRatio","std::vector<float>",&pfCluster_full5x5_e2x5BottomRatio); 
      tree->Branch("pfCluster_full5x5_e2x5LeftRatio","std::vector<float>",&pfCluster_full5x5_e2x5LeftRatio); 
      tree->Branch("pfCluster_full5x5_swissCross","std::vector<float>",&pfCluster_full5x5_swissCross); 
      tree->Branch("pfCluster_full5x5_r9","std::vector<float>",&pfCluster_full5x5_r9);
      tree->Branch("pfCluster_full5x5_sigmaIetaIeta","std::vector<float>",&pfCluster_full5x5_sigmaIetaIeta);
      tree->Branch("pfCluster_full5x5_sigmaIetaIphi","std::vector<float>",&pfCluster_full5x5_sigmaIetaIphi);
      tree->Branch("pfCluster_full5x5_sigmaIphiIphi","std::vector<float>",&pfCluster_full5x5_sigmaIphiIphi);
   }
   if(saveSuperCluster_ && saveShowerShapes_){  
      tree->Branch("superCluster_e5x5","std::vector<float>",&superCluster_e5x5);
      tree->Branch("superCluster_e2x2Ratio","std::vector<float>",&superCluster_e2x2Ratio);
      tree->Branch("superCluster_e3x3Ratio","std::vector<float>",&superCluster_e3x3Ratio);
      tree->Branch("superCluster_eMaxRatio","std::vector<float>",&superCluster_eMaxRatio);
      tree->Branch("superCluster_e2ndRatio","std::vector<float>",&superCluster_e2ndRatio);
      tree->Branch("superCluster_eTopRatio","std::vector<float>",&superCluster_eTopRatio);
      tree->Branch("superCluster_eRightRatio","std::vector<float>",&superCluster_eRightRatio);
      tree->Branch("superCluster_eBottomRatio","std::vector<float>",&superCluster_eBottomRatio);
      tree->Branch("superCluster_eLeftRatio","std::vector<float>",&superCluster_eLeftRatio);
      tree->Branch("superCluster_e2x5MaxRatio","std::vector<float>",&superCluster_e2x5MaxRatio);
      tree->Branch("superCluster_e2x5TopRatio","std::vector<float>",&superCluster_e2x5TopRatio);
      tree->Branch("superCluster_e2x5RightRatio","std::vector<float>",&superCluster_e2x5RightRatio);
      tree->Branch("superCluster_e2x5BottomRatio","std::vector<float>",&superCluster_e2x5BottomRatio); 
      tree->Branch("superCluster_e2x5LeftRatio","std::vector<float>",&superCluster_e2x5LeftRatio); 
      tree->Branch("superCluster_swissCross","std::vector<float>",&superCluster_swissCross); 
      tree->Branch("superCluster_r9","std::vector<float>",&superCluster_r9);
      tree->Branch("superCluster_sigmaIetaIeta","std::vector<float>",&superCluster_sigmaIetaIeta);
      tree->Branch("superCluster_sigmaIetaIphi","std::vector<float>",&superCluster_sigmaIetaIphi);
      tree->Branch("superCluster_sigmaIphiIphi","std::vector<float>",&superCluster_sigmaIphiIphi);
      tree->Branch("superCluster_full5x5_e5x5","std::vector<float>",&superCluster_full5x5_e5x5);
      tree->Branch("superCluster_full5x5_e2x2Ratio","std::vector<float>",&superCluster_full5x5_e2x2Ratio);
      tree->Branch("superCluster_full5x5_e3x3Ratio","std::vector<float>",&superCluster_full5x5_e3x3Ratio);
      tree->Branch("superCluster_full5x5_eMaxRatio","std::vector<float>",&superCluster_full5x5_eMaxRatio);
      tree->Branch("superCluster_full5x5_e2ndRatio","std::vector<float>",&superCluster_full5x5_e2ndRatio);
      tree->Branch("superCluster_full5x5_eTopRatio","std::vector<float>",&superCluster_full5x5_eTopRatio);
      tree->Branch("superCluster_full5x5_eRightRatio","std::vector<float>",&superCluster_full5x5_eRightRatio);
      tree->Branch("superCluster_full5x5_eBottomRatio","std::vector<float>",&superCluster_full5x5_eBottomRatio);
      tree->Branch("superCluster_full5x5_eLeftRatio","std::vector<float>",&superCluster_full5x5_eLeftRatio);
      tree->Branch("superCluster_full5x5_e2x5MaxRatio","std::vector<float>",&superCluster_full5x5_e2x5MaxRatio);
      tree->Branch("superCluster_full5x5_e2x5TopRatio","std::vector<float>",&superCluster_full5x5_e2x5TopRatio);
      tree->Branch("superCluster_full5x5_e2x5RightRatio","std::vector<float>",&superCluster_full5x5_e2x5RightRatio);
      tree->Branch("superCluster_full5x5_e2x5BottomRatio","std::vector<float>",&superCluster_full5x5_e2x5BottomRatio); 
      tree->Branch("superCluster_full5x5_e2x5LeftRatio","std::vector<float>",&superCluster_full5x5_e2x5LeftRatio); 
      tree->Branch("superCluster_full5x5_swissCross","std::vector<float>",&superCluster_full5x5_swissCross); 
      tree->Branch("superCluster_full5x5_r9","std::vector<float>",&superCluster_full5x5_r9);
      tree->Branch("superCluster_full5x5_sigmaIetaIeta","std::vector<float>",&superCluster_full5x5_sigmaIetaIeta);
      tree->Branch("superCluster_full5x5_sigmaIetaIphi","std::vector<float>",&superCluster_full5x5_sigmaIetaIphi);
      tree->Branch("superCluster_full5x5_sigmaIphiIphi","std::vector<float>",&superCluster_full5x5_sigmaIphiIphi);
      if(useHcalTowers_ ){
         tree->Branch("superCluster_HoEraw","std::vector<float>",&superCluster_HoEraw);
         tree->Branch("superCluster_HoErawBC","std::vector<float>",&superCluster_HoErawBC);
        }
   }
}

RecoSimDumperPF::~RecoSimDumperPF()
{
        // do anything here that needs to be done at desctruction time
        // (e.g. close files, deallocate resources etc.)
}


//
// member functions
//

// ------------ method called to for each event  ------------
void RecoSimDumperPF::analyze(const edm::Event& ev, const edm::EventSetup& iSetup)
{

   //calo geometry
   edm::ESHandle<CaloGeometry> caloGeometry;
   iSetup.get<CaloGeometryRecord>().get(caloGeometry);
   const CaloGeometry *geometry = caloGeometry.product();
   _ebGeom = caloGeometry->getSubdetectorGeometry(DetId::Ecal, EcalBarrel);
   _eeGeom = caloGeometry->getSubdetectorGeometry(DetId::Ecal, EcalEndcap);
   _esGeom = caloGeometry->getSubdetectorGeometry(DetId::Ecal, EcalPreshower);
   if (_esGeom) {
    for (uint32_t ic = 0; ic < _esGeom->getValidDetIds().size() && (!_esPlus || !_esMinus); ++ic) {
      const double z = _esGeom->getGeometry(_esGeom->getValidDetIds()[ic])->getPosition().z();
      _esPlus = _esPlus || (0 < z);
      _esMinus = _esMinus || (0 > z);
    }
   }

   edm::ESHandle<CaloTopology> caloTopology;
   iSetup.get<CaloTopologyRecord>().get(caloTopology);
   const CaloTopology* topology = caloTopology.product();
   
   edm::Handle<double> rhos;
   ev.getByToken(rhoToken_,rhos);
   if (!rhos.isValid()) {
       std::cerr << "Analyze --> rhos not found" << std::endl; 
       return;
   }

   edm::Handle<reco::VertexCollection> vertices;
   ev.getByToken(vtxToken_,vertices);
   if (!vertices.isValid()) {
       std::cerr << "Analyze --> vertices not found" << std::endl; 
       return;
   }

   edm::Handle<std::vector<reco::GenParticle> > genParticles;
   ev.getByToken(genToken_,genParticles);
   if (!genParticles.isValid()) {
       std::cerr << "Analyze --> genParticles not found" << std::endl; 
       return;
   }

   edm::Handle<std::vector<CaloParticle> > caloParticles;
   ev.getByToken(caloPartToken_,caloParticles);
   if (!caloParticles.isValid()) {
       std::cerr << "Analyze --> caloParticles not found" << std::endl; 
       return;
   }

   edm::Handle<EcalRecHitCollection> recHitsEB;
   ev.getByToken(ebRechitToken_, recHitsEB);
   if (!recHitsEB.isValid()) {
       std::cerr << "Analyze --> recHitsEB not found" << std::endl; 
       return;
   }


   edm::Handle<EcalRecHitCollection> recHitsEE;
   ev.getByToken(eeRechitToken_, recHitsEE);
   if (!recHitsEE.isValid()) {
        std::cerr << "Analyze --> recHitsEE not found" << std::endl; 
        return;
   }

   edm::Handle<std::vector<reco::PFRecHit> > pfRecHits;
   ev.getByToken(pfRecHitToken_, pfRecHits);
   if(saveEBPFRechits_) {
      if (!pfRecHits.isValid()) {
          std::cerr << "Analyze --> pfRecHits not found" << std::endl; 
          return;
      }
   } 

   edm::Handle<std::vector<reco::PFCluster> > pfClusters;
   ev.getByToken(pfClusterToken_, pfClusters);
   if(savePFCluster_) {
      if (!pfClusters.isValid()) {
          std::cerr << "Analyze --> pfClusters not found" << std::endl; 
          return;
      }
   } 

   edm::Handle<std::vector<reco::SuperCluster> > superClusterEB;
   ev.getByToken(ebSuperClusterToken_, superClusterEB);
   if(saveSuperCluster_) {
      if (!superClusterEB.isValid()) {
          std::cerr << "Analyze --> superClusterEB not found" << std::endl; 
          return;
      }
   } 

   edm::Handle<std::vector<reco::SuperCluster> > superClusterEE;
   ev.getByToken(eeSuperClusterToken_, superClusterEE);
   if(saveSuperCluster_) {
      if (!superClusterEE.isValid()) {
          std::cerr << "Analyze --> superClusterEE not found" << std::endl; 
          return;
      }
   }

   //compute EgammaTowers;
   Handle<CaloTowerCollection> hcalTowers;
   ev.getByToken(hcalTowersToken_, hcalTowers);
   if(useHcalTowers_){
      if (!hcalTowers.isValid()) {
          std::cerr << "Analyze --> hcalTowers not found" << std::endl; 
          return;
      } 
      //hcalTowersColl = hcalTowers.product();
      towerIso1_ = new EgammaTowerIsolation(0.15, 0., 0., 1, hcalTowers.product());
      towerIso2_ = new EgammaTowerIsolation(0.15, 0., 0., 2, hcalTowers.product());
      egammaHadTower_ = new EgammaHadTower(iSetup);
      egammaHadTower_->setTowerCollection(hcalTowers.product()); 
   } 

   runId = ev.id().run();
   lumiId = ev.luminosityBlock();
   eventId = ev.id().event();

   nVtx = vertices->size();
   rho = *(rhos.product());

   genParticle_id.clear();
   genParticle_isGammaFromMeson.clear();
   genParticle_energy.clear();
   genParticle_pt.clear();
   genParticle_eta.clear();
   genParticle_phi.clear();
   std::vector<GenParticle> genParts;
   int counter=0;
   for(const auto& iGen : *(genParticles.product()))
   {
       bool isGoodParticle = false; 
       for(unsigned int id=0; id<genID_.size(); id++)
           if((iGen.pdgId()==genID_.at(id) || genID_.at(id)==0) && iGen.status()==1) isGoodParticle=true;
      
       if(!isGoodParticle) continue; 
       genParticle_id.push_back(iGen.pdgId()); 
       genParticle_energy.push_back(iGen.energy()); 
       genParticle_pt.push_back(iGen.pt());
       genParticle_eta.push_back(iGen.eta());
       genParticle_phi.push_back(iGen.phi());
       
       genParts.push_back(iGen); 

       // try to determine if it comes from a meson => go two generations up and check pdg id...
       bool isGammaFromMeson=false;
       if(iGen.pdgId()!=22) continue;
       if(iGen.pt()<5) continue;
       counter++;
       //std::vector<const GenParticle*>::iterator itGen;
       GenParticle *iGenClone = iGen.clone(); 
       for(unsigned int iM=0; iM<iGenClone->numberOfMothers(); iM++){
         if(isGammaFromMeson==true) break;
         const Candidate *this_mother = iGenClone->mother(iM);
         if(this_mother->pdgId()==111 || this_mother->pdgId()==221) {
           isGammaFromMeson=true;
         }else{
           for(unsigned int iMU=0; iMU<this_mother->numberOfMothers(); iMU++){
             if(isGammaFromMeson==true) break;
             const Candidate *this_motherUp = this_mother->mother(iM);
             if(this_motherUp->pdgId()==111 || this_motherUp->pdgId()==221) {
               isGammaFromMeson=true;
             }
           }
         }
       } // end loop over mothers 
       genParticle_isGammaFromMeson.push_back(isGammaFromMeson);
   } // loop over gen particles 
   
   int nGenParticles = genParts.size(); 
   
   std::vector<CaloParticle> caloParts;
   for(const auto& iCalo : *(caloParticles.product()))
   {
       bool isGoodParticle = false; 
       for(unsigned int id=0; id<genID_.size(); id++) 
           if(iCalo.pdgId()==genID_.at(id) || genID_.at(id)==0) isGoodParticle=true;

       if(!isGoodParticle) continue;     

       caloParts.push_back(iCalo); 
   }

   int nCaloParticles = caloParts.size(); 
   //std::cout << "CaloParticles size  : " << nCaloParticles << std::endl;
  
   genParticle_pfCluster_dR_genScore_MatchedIndex.clear();
   genParticle_pfCluster_dR_genScore_MatchedIndex.resize(nGenParticles);

   caloParticle_id.clear();
   caloParticle_genEnergy.clear();
   caloParticle_simEnergy.clear();
   caloParticle_genPt.clear();
   caloParticle_simPt.clear();
   caloParticle_genEta.clear();
   caloParticle_simEta.clear();
   caloParticle_genPhi.clear();
   caloParticle_simPhi.clear();
   caloParticle_simIeta.clear();
   caloParticle_simIphi.clear();
   caloParticle_simIz.clear();
   caloParticle_pfCluster_dR_simScore_MatchedIndex.clear();
   caloParticle_pfCluster_sim_fraction_old_MatchedIndex.clear();
   caloParticle_pfCluster_simScore_MatchedIndex.clear();
   caloParticle_pfCluster_n_shared_xtals_MatchedIndex.clear();
   caloParticle_pfCluster_sim_fraction_MatchedIndex.clear();
   caloParticle_pfCluster_sim_fraction_1MeVCut_MatchedIndex.clear();
   caloParticle_pfCluster_sim_fraction_5MeVCut_MatchedIndex.clear();
   caloParticle_pfCluster_sim_fraction_10MeVCut_MatchedIndex.clear();
   caloParticle_pfCluster_sim_fraction_50MeVCut_MatchedIndex.clear();
   caloParticle_pfCluster_sim_fraction_100MeVCut_MatchedIndex.clear();  
   caloParticle_pfCluster_sim_fraction_500MeVCut_MatchedIndex.clear();  
   caloParticle_pfCluster_sim_fraction_1GeVCut_MatchedIndex.clear();  
   caloParticle_pfCluster_sim_rechit_diff_MatchedIndex.clear();
   caloParticle_pfCluster_sim_rechit_fraction_MatchedIndex.clear();
   caloParticle_pfCluster_global_sim_rechit_fraction_MatchedIndex.clear();
   caloParticle_pfCluster_hgcal_caloToCluster_MatchedIndex.clear();
   caloParticle_pfCluster_hgcal_clusterToCalo_MatchedIndex.clear();
   caloParticle_pfCluster_sim_rechit_combined_fraction_MatchedIndex.clear();
   caloParticle_pfCluster_rechit_sim_combined_fraction_MatchedIndex.clear();
   caloParticle_pfCluster_dR_simScore_MatchedIndex.resize(nCaloParticles);
   caloParticle_pfCluster_sim_fraction_old_MatchedIndex.resize(nCaloParticles);
   caloParticle_pfCluster_simScore_MatchedIndex.resize(nCaloParticles);
   caloParticle_pfCluster_n_shared_xtals_MatchedIndex.resize(nCaloParticles);
   caloParticle_pfCluster_sim_fraction_MatchedIndex.resize(nCaloParticles);
   caloParticle_pfCluster_sim_fraction_1MeVCut_MatchedIndex.resize(nCaloParticles);
   caloParticle_pfCluster_sim_fraction_5MeVCut_MatchedIndex.resize(nCaloParticles);
   caloParticle_pfCluster_sim_fraction_10MeVCut_MatchedIndex.resize(nCaloParticles);
   caloParticle_pfCluster_sim_fraction_50MeVCut_MatchedIndex.resize(nCaloParticles);
   caloParticle_pfCluster_sim_fraction_100MeVCut_MatchedIndex.resize(nCaloParticles);
   caloParticle_pfCluster_sim_fraction_500MeVCut_MatchedIndex.resize(nCaloParticles);
   caloParticle_pfCluster_sim_fraction_1GeVCut_MatchedIndex.resize(nCaloParticles);
   caloParticle_pfCluster_sim_rechit_diff_MatchedIndex.resize(nCaloParticles);
   caloParticle_pfCluster_sim_rechit_fraction_MatchedIndex.resize(nCaloParticles);
   caloParticle_pfCluster_global_sim_rechit_fraction_MatchedIndex.resize(nCaloParticles);
   caloParticle_pfCluster_hgcal_caloToCluster_MatchedIndex.resize(nCaloParticles);
   caloParticle_pfCluster_hgcal_clusterToCalo_MatchedIndex.resize(nCaloParticles);  
   caloParticle_pfCluster_sim_rechit_combined_fraction_MatchedIndex.resize(nCaloParticles);
   caloParticle_pfCluster_rechit_sim_combined_fraction_MatchedIndex.resize(nCaloParticles);  
   
   simHit_energy.clear();
   simHit_eta.clear();
   simHit_phi.clear();
   simHit_ieta.clear();
   simHit_iphi.clear();
   simHit_iz.clear();
   simHit_icP.clear();

   pfRecHit_energy.clear();
   pfRecHit_eta.clear();
   pfRecHit_phi.clear();
   pfRecHit_ieta.clear();
   pfRecHit_iphi.clear();
   pfRecHit_iz.clear();

   int nPFClusters = (pfClusters.product())->size();
   pfCluster_energy.clear();
   pfCluster_eta.clear();
   pfCluster_phi.clear();
   pfCluster_ieta.clear();
   pfCluster_iphi.clear();
   pfCluster_iz.clear();
   pfCluster_superClustersIndex.clear();
   pfCluster_e5x5.clear();
   pfCluster_e2x2Ratio.clear();
   pfCluster_e3x3Ratio.clear();
   pfCluster_eMaxRatio.clear();
   pfCluster_e2ndRatio.clear();
   pfCluster_eTopRatio.clear();
   pfCluster_eRightRatio.clear();
   pfCluster_eBottomRatio.clear();
   pfCluster_eLeftRatio.clear();
   pfCluster_e2x5MaxRatio.clear();
   pfCluster_e2x5TopRatio.clear();
   pfCluster_e2x5RightRatio.clear();
   pfCluster_e2x5BottomRatio.clear();
   pfCluster_e2x5LeftRatio.clear();
   pfCluster_swissCross.clear();
   pfCluster_r9.clear();
   pfCluster_sigmaIetaIeta.clear(); 
   pfCluster_sigmaIetaIphi.clear(); 
   pfCluster_sigmaIphiIphi.clear(); 
   pfCluster_full5x5_e5x5.clear();
   pfCluster_full5x5_e2x2Ratio.clear();
   pfCluster_full5x5_e3x3Ratio.clear();
   pfCluster_full5x5_eMaxRatio.clear();
   pfCluster_full5x5_e2ndRatio.clear();
   pfCluster_full5x5_eTopRatio.clear();
   pfCluster_full5x5_eRightRatio.clear();
   pfCluster_full5x5_eBottomRatio.clear();
   pfCluster_full5x5_eLeftRatio.clear();
   pfCluster_full5x5_e2x5MaxRatio.clear();
   pfCluster_full5x5_e2x5TopRatio.clear();
   pfCluster_full5x5_e2x5RightRatio.clear();
   pfCluster_full5x5_e2x5BottomRatio.clear();
   pfCluster_full5x5_e2x5LeftRatio.clear();
   pfCluster_full5x5_swissCross.clear();
   pfCluster_full5x5_r9.clear();
   pfCluster_full5x5_sigmaIetaIeta.clear(); 
   pfCluster_full5x5_sigmaIetaIphi.clear(); 
   pfCluster_full5x5_sigmaIphiIphi.clear();    
   pfCluster_dR_genScore_MatchedIndex.clear();
   pfCluster_dR_simScore_MatchedIndex.clear();
   pfCluster_sim_fraction_old_MatchedIndex.clear();
   pfCluster_simScore_MatchedIndex.clear();
   pfCluster_n_shared_xtals_MatchedIndex.clear();
   pfCluster_sim_fraction_MatchedIndex.clear();
   pfCluster_sim_fraction_1MeVCut_MatchedIndex.clear();
   pfCluster_sim_fraction_5MeVCut_MatchedIndex.clear();
   pfCluster_sim_fraction_10MeVCut_MatchedIndex.clear();
   pfCluster_sim_fraction_50MeVCut_MatchedIndex.clear();
   pfCluster_sim_fraction_100MeVCut_MatchedIndex.clear();
   pfCluster_sim_fraction_500MeVCut_MatchedIndex.clear();
   pfCluster_sim_fraction_1GeVCut_MatchedIndex.clear();
   pfCluster_sim_rechit_diff_MatchedIndex.clear();
   pfCluster_sim_rechit_fraction_MatchedIndex.clear();
   pfCluster_global_sim_rechit_fraction_MatchedIndex.clear();  
   pfCluster_hgcal_caloToCluster_MatchedIndex.clear();  
   pfCluster_hgcal_clusterToCalo_MatchedIndex.clear();   
   pfCluster_sim_rechit_combined_fraction_MatchedIndex.clear();  
   pfCluster_rechit_sim_combined_fraction_MatchedIndex.clear();   
   pfCluster_dR_genScore.clear();
   pfCluster_dR_simScore.clear();
   pfCluster_sim_fraction_old.clear();
   pfCluster_simScore.clear();
   pfCluster_n_shared_xtals.clear();
   pfCluster_sim_fraction.clear();
   pfCluster_sim_fraction_1MeVCut.clear();
   pfCluster_sim_fraction_5MeVCut.clear();
   pfCluster_sim_fraction_10MeVCut.clear();
   pfCluster_sim_fraction_50MeVCut.clear();
   pfCluster_sim_fraction_100MeVCut.clear();
   pfCluster_sim_fraction_500MeVCut.clear();
   pfCluster_sim_fraction_1GeVCut.clear();
   pfCluster_sim_rechit_diff.clear();
   pfCluster_sim_rechit_fraction.clear();
   pfCluster_global_sim_rechit_fraction.clear();  
   pfCluster_hgcal_caloToCluster.clear();  
   pfCluster_hgcal_clusterToCalo.clear();    
   pfCluster_sim_rechit_combined_fraction.clear();  
   pfCluster_rechit_sim_combined_fraction.clear();    
   pfCluster_dR_genScore.resize(nPFClusters);
   pfCluster_dR_simScore.resize(nPFClusters);
   pfCluster_sim_fraction_old.resize(nPFClusters);
   pfCluster_simScore.resize(nPFClusters);
   pfCluster_n_shared_xtals.resize(nPFClusters);
   pfCluster_sim_fraction.resize(nPFClusters);
   pfCluster_sim_fraction_1MeVCut.resize(nPFClusters);
   pfCluster_sim_fraction_5MeVCut.resize(nPFClusters);
   pfCluster_sim_fraction_10MeVCut.resize(nPFClusters);
   pfCluster_sim_fraction_50MeVCut.resize(nPFClusters);
   pfCluster_sim_fraction_100MeVCut.resize(nPFClusters);
   pfCluster_sim_fraction_500MeVCut.resize(nPFClusters);
   pfCluster_sim_fraction_1GeVCut.resize(nPFClusters);
   pfCluster_sim_rechit_diff.resize(nPFClusters);
   pfCluster_sim_rechit_fraction.resize(nPFClusters);
   pfCluster_global_sim_rechit_fraction.resize(nPFClusters); 
   pfCluster_hgcal_caloToCluster.resize(nPFClusters);   
   pfCluster_hgcal_clusterToCalo.resize(nPFClusters);  
   pfCluster_sim_rechit_combined_fraction.resize(nPFClusters);   
   pfCluster_rechit_sim_combined_fraction.resize(nPFClusters);  
   pfCluster_superClustersIndex.resize(nPFClusters); 

   pfClusterHit_energy.clear();
   pfClusterHit_rechitEnergy.clear();
   pfClusterHit_eta.clear();
   pfClusterHit_phi.clear();   
   pfClusterHit_ieta.clear();
   pfClusterHit_iphi.clear(); 
   pfClusterHit_iz.clear();           
   pfClusterHit_energy.resize(nPFClusters);  
   pfClusterHit_rechitEnergy.resize(nPFClusters);  
   pfClusterHit_eta.resize(nPFClusters);  
   pfClusterHit_phi.resize(nPFClusters);     
   pfClusterHit_ieta.resize(nPFClusters);  
   pfClusterHit_iphi.resize(nPFClusters);   
   pfClusterHit_iz.resize(nPFClusters);   
  
   superCluster_energy.clear(); 
   superCluster_eta.clear(); 
   superCluster_phi.clear();  
   superCluster_etaWidth.clear();     
   superCluster_phiWidth.clear(); 
   superCluster_R.clear(); 
   superCluster_nPFClusters.clear(); 
   superCluster_ieta.clear(); 
   superCluster_iphi.clear();
   superCluster_iz.clear(); 
   superCluster_seedIndex.clear(); 
   superCluster_pfClustersIndex.clear(); 
   superCluster_e5x5.clear();
   superCluster_e2x2Ratio.clear();
   superCluster_e3x3Ratio.clear();
   superCluster_eMaxRatio.clear();
   superCluster_e2ndRatio.clear();
   superCluster_eTopRatio.clear();
   superCluster_eRightRatio.clear();
   superCluster_eBottomRatio.clear();
   superCluster_eLeftRatio.clear();
   superCluster_e2x5MaxRatio.clear();
   superCluster_e2x5TopRatio.clear();
   superCluster_e2x5RightRatio.clear();
   superCluster_e2x5BottomRatio.clear();
   superCluster_e2x5LeftRatio.clear();
   superCluster_swissCross.clear();
   superCluster_r9.clear();
   superCluster_sigmaIetaIeta.clear(); 
   superCluster_sigmaIetaIphi.clear(); 
   superCluster_sigmaIphiIphi.clear(); 
   superCluster_full5x5_e5x5.clear();
   superCluster_full5x5_e2x2Ratio.clear();
   superCluster_full5x5_e3x3Ratio.clear();
   superCluster_full5x5_eMaxRatio.clear();
   superCluster_full5x5_e2ndRatio.clear();
   superCluster_full5x5_eTopRatio.clear();
   superCluster_full5x5_eRightRatio.clear();
   superCluster_full5x5_eBottomRatio.clear();
   superCluster_full5x5_eLeftRatio.clear();
   superCluster_full5x5_e2x5MaxRatio.clear();
   superCluster_full5x5_e2x5TopRatio.clear();
   superCluster_full5x5_e2x5RightRatio.clear();
   superCluster_full5x5_e2x5BottomRatio.clear();
   superCluster_full5x5_e2x5LeftRatio.clear();
   superCluster_full5x5_swissCross.clear();
   superCluster_full5x5_r9.clear();
   superCluster_full5x5_sigmaIetaIeta.clear(); 
   superCluster_full5x5_sigmaIetaIphi.clear(); 
   superCluster_full5x5_sigmaIphiIphi.clear();    
   int nSuperClusters = (superClusterEB.product())->size() + (superClusterEE.product())->size();
   superCluster_seedIndex.resize(nSuperClusters); 
   superCluster_pfClustersIndex.resize(nSuperClusters);
   hitsAndEnergies_CaloPart.clear();
   hitsAndEnergies_CaloPart_1MeVCut.clear();
   hitsAndEnergies_CaloPart_5MeVCut.clear();
   hitsAndEnergies_CaloPart_10MeVCut.clear();
   hitsAndEnergies_CaloPart_50MeVCut.clear();
   hitsAndEnergies_CaloPart_100MeVCut.clear();
   hitsAndEnergies_PFCluster.clear();

   GlobalPoint caloParticle_position;
   GlobalPoint cell;

   for(unsigned int iCalo=0; iCalo<caloParts.size(); iCalo++){
   
       const auto& genParticles_caloPart = caloParts.at(iCalo).genParticles();
       caloParticle_id.push_back(caloParts.at(iCalo).pdgId());
       if(genParticles_caloPart.empty()){
          cout << "WARNING: no associated genParticle found, making standard dR matching" << endl;
          float dR=999.;
          int igen_tmp=-1; 
          int igen=0; 
          for(const auto& iGen : *(genParticles.product()))
          {
              float dR_tmp = deltaR(caloParts.at(iCalo).eta(),caloParts.at(iCalo).phi(),iGen.eta(),iGen.phi());  
              if(dR_tmp<dR && iGen.status()==1){
                 dR=dR_tmp;
                 igen_tmp=igen;
              }  
              igen++;
          } 
          const auto& genParticles_tmp = *(genParticles.product());
          auto genParticle = genParticles_tmp[igen_tmp]; 
          caloParticle_genEnergy.push_back(reduceFloat(genParticle.energy(),nBits_));
          caloParticle_genPt.push_back(reduceFloat(genParticle.pt(),nBits_));
          caloParticle_genEta.push_back(reduceFloat(genParticle.eta(),nBits_));
          caloParticle_genPhi.push_back(reduceFloat(genParticle.phi(),nBits_));
       }else{
          caloParticle_genEnergy.push_back(reduceFloat((*genParticles_caloPart.begin())->energy(),nBits_));
          caloParticle_genPt.push_back(reduceFloat((*genParticles_caloPart.begin())->pt(),nBits_));
          caloParticle_genEta.push_back(reduceFloat((*genParticles_caloPart.begin())->eta(),nBits_));
          caloParticle_genPhi.push_back(reduceFloat((*genParticles_caloPart.begin())->phi(),nBits_));
       }
      
       hitsAndEnergies_CaloPart.push_back(*getHitsAndEnergiesCaloPart(&(caloParts.at(iCalo)),-1.));
       hitsAndEnergies_CaloPart_1MeVCut.push_back(*getHitsAndEnergiesCaloPart(&(caloParts.at(iCalo)),0.001)); 
       hitsAndEnergies_CaloPart_5MeVCut.push_back(*getHitsAndEnergiesCaloPart(&(caloParts.at(iCalo)),0.005));    
       hitsAndEnergies_CaloPart_10MeVCut.push_back(*getHitsAndEnergiesCaloPart(&(caloParts.at(iCalo)),0.01)); 
       hitsAndEnergies_CaloPart_50MeVCut.push_back(*getHitsAndEnergiesCaloPart(&(caloParts.at(iCalo)),0.05)); 
       hitsAndEnergies_CaloPart_100MeVCut.push_back(*getHitsAndEnergiesCaloPart(&(caloParts.at(iCalo)),0.1));     
       GlobalPoint caloParticle_position = calculateAndSetPositionActual(&hitsAndEnergies_CaloPart.at(iCalo), 7.4, 3.1, 1.2, 4.2, 0.89, 0.,false);
       if (caloParticle_position == GlobalPoint(-999999., -999999., -999999.)) {
          std::cout << "Invalid position for caloparticle, skipping event" << std::endl;
          return;
       }
       caloParticle_simEta.push_back(reduceFloat(caloParticle_position.eta(),nBits_));
       caloParticle_simPhi.push_back(reduceFloat(caloParticle_position.phi(),nBits_));
       caloParticle_simPt.push_back(reduceFloat(sqrt(caloParts.at(iCalo).px()*caloParts.at(iCalo).px() + caloParts.at(iCalo).py()*caloParts.at(iCalo).py() + caloParts.at(iCalo).pz()*caloParts.at(iCalo).pz())/TMath::CosH(caloParticle_position.eta()),nBits_));   
       if(std::abs(caloParticle_position.eta()) < 1.479){  
          EBDetId eb_id(_ebGeom->getClosestCell(caloParticle_position));  
          caloParticle_simIeta.push_back(eb_id.ieta());
          caloParticle_simIphi.push_back(eb_id.iphi());
          caloParticle_simIz.push_back(0); 
       }else{            
          int iz=-99;
          EEDetId ee_id(_eeGeom->getClosestCell(caloParticle_position));   
          caloParticle_simIeta.push_back(ee_id.ix());
          caloParticle_simIphi.push_back(ee_id.iy());
          if(ee_id.zside()<0) iz=-1;
          if(ee_id.zside()>0) iz=1;  
          caloParticle_simIz.push_back(iz); 
       }   
   }
   
   //save hitsAndEnergies for each CaloParticle, PFcluster and SuperCluster
   if(saveCaloParticles_){
      hits_CaloPart.resize(nCaloParticles);
      for(unsigned int iCalo = 0; iCalo < hitsAndEnergies_CaloPart.size(); iCalo++){
          for(unsigned int i = 0; i < hitsAndEnergies_CaloPart.at(iCalo).size(); i++)
              hits_CaloPart[iCalo].push_back(hitsAndEnergies_CaloPart.at(iCalo).at(i).first);
      }
   } 

   if(savePFCluster_){
      hits_PFCluster.resize(nPFClusters);
      for(const auto& iPFCluster : *(pfClusters.product())){  
          reco::CaloCluster caloBC(iPFCluster);
          hitsAndEnergies_PFCluster.push_back(*getHitsAndEnergiesBC(&caloBC,recHitsEB,recHitsEE));
      }
      for(unsigned int iPFCl = 0; iPFCl < hitsAndEnergies_PFCluster.size(); iPFCl++){
          for(unsigned int i = 0; i < hitsAndEnergies_PFCluster.at(iPFCl).size(); i++)
              hits_PFCluster[iPFCl].push_back(hitsAndEnergies_PFCluster.at(iPFCl).at(i).first);
      }
   }
        
   //save simhits information
   for(unsigned int iCaloCount=0; iCaloCount<hitsAndEnergies_CaloPart.size(); iCaloCount++) 
   {
       float calo_simEnergy=0.;

       for(auto const& hit: hitsAndEnergies_CaloPart[iCaloCount])
       {
           DetId id(hit.first);
           if(id.subdetId()!=EcalBarrel && id.subdetId()!=EcalEndcap) continue;
               
           calo_simEnergy += hit.second; 

           cell = geometry->getPosition(id);
           float eta = cell.eta();  
           float phi = cell.phi();  
           int ieta = -99; 
           int iphi = -99;
           int iz = -99;  
           if(id.subdetId()==EcalBarrel){
              EBDetId eb_id(id);
              ieta = eb_id.ieta(); 
              iphi = eb_id.iphi();  
              iz = 0;   
           }
           if(id.subdetId()==EcalEndcap){
              EEDetId ee_id(id);
              ieta = ee_id.ix(); 
              iphi = ee_id.iy();
              if(ee_id.zside()<0) iz=-1;
              if(ee_id.zside()>0) iz=1;   
           }

           if(saveSimhits_){
              simHit_energy.push_back(reduceFloat(hit.second,nBits_));
              simHit_eta.push_back(reduceFloat(eta,nBits_));
              simHit_phi.push_back(reduceFloat(phi,nBits_));
              simHit_ieta.push_back(ieta);
              simHit_iphi.push_back(iphi);
              simHit_iz.push_back(iz); 
              simHit_icP.push_back(iCaloCount); 
           }
       } 
       caloParticle_simEnergy.push_back(reduceFloat(calo_simEnergy,nBits_));
   }
  
   // Save all EB PFRechits
   if(saveEBPFRechits_){
     for(const auto& iPFRechit : *(pfRecHits.product())){

       DetId pf_id(iPFRechit.detId());
       //pfRechit.push_back(pf_id); 
       cell = geometry->getPosition(pf_id); 
       if(pf_id.subdetId()==EcalBarrel){ 
          pfRecHit_energy.push_back(reduceFloat(iPFRechit.energy(),nBits_));    
          pfRecHit_eta.push_back(reduceFloat(cell.eta(),nBits_));  
          pfRecHit_phi.push_back(reduceFloat(cell.phi(),nBits_)); 
          EBDetId eb_id(pf_id);  
          pfRecHit_ieta.push_back(eb_id.ieta());  
          pfRecHit_iphi.push_back(eb_id.iphi());  
          pfRecHit_iz.push_back(0);     
       }
       /*else if(pf_id.subdetId()==EcalEndcap){
          int iz=-99;
          EEDetId ee_id(pf_id);  
          if(ee_id.zside()<0) iz=-1;
          if(ee_id.zside()>0) iz=1; 
          pfRecHit_ieta.push_back(ee_id.ix());  
          pfRecHit_iphi.push_back(ee_id.iy());  
          pfRecHit_iz.push_back(iz);    
       }*/ 
     } 
   }

   //Save PFClusters 
   if(savePFCluster_){
      
      int iPFCl=0;
      //std::cout << "PFClusters size     : " << (pfClusters.product())->size() << std::endl;
      for(const auto& iPFCluster : *(pfClusters.product())){  

          dR_genScore.clear();
          dR_simScore.clear();
          sim_fraction_old.clear();
          simScore.clear();
          n_shared_xtals.clear();
          sim_fraction.clear();
          sim_fraction_1MeVCut.clear();
          sim_fraction_5MeVCut.clear();
          sim_fraction_10MeVCut.clear();
          sim_fraction_50MeVCut.clear();
          sim_fraction_100MeVCut.clear();  
          sim_fraction_500MeVCut.clear(); 
          sim_fraction_1GeVCut.clear();    
          sim_rechit_diff.clear();
          sim_rechit_fraction.clear();
          global_sim_rechit_fraction.clear();
          hgcal_caloToCluster.clear();
          hgcal_clusterToCalo.clear();    
          sim_rechit_combined_fraction.clear();
          rechit_sim_combined_fraction.clear();

          pfCluster_energy.push_back(reduceFloat(iPFCluster.energy(),nBits_));
          pfCluster_eta.push_back(reduceFloat(iPFCluster.eta(),nBits_));
          pfCluster_phi.push_back(reduceFloat(iPFCluster.phi(),nBits_));
          reco::CaloCluster caloBC(iPFCluster);

          math::XYZPoint caloPos = caloBC.position();
          if(iPFCluster.layer() == PFLayer::ECAL_BARREL){  
             EBDetId eb_id(_ebGeom->getClosestCell(GlobalPoint(caloPos.x(),caloPos.y(),caloPos.z())));  
             pfCluster_ieta.push_back(eb_id.ieta());
             pfCluster_iphi.push_back(eb_id.iphi());
             pfCluster_iz.push_back(0); 
          }else if(iPFCluster.layer() == PFLayer::ECAL_ENDCAP){ 
             int iz=-99;
             EEDetId ee_id(_eeGeom->getClosestCell(GlobalPoint(caloPos.x(),caloPos.y(),caloPos.z())));  
             if(ee_id.zside()<0) iz=-1;
             if(ee_id.zside()>0) iz=1;     
             pfCluster_ieta.push_back(ee_id.ix());
             pfCluster_iphi.push_back(ee_id.iy());
             pfCluster_iz.push_back(iz); 
          } 
           
          if(saveShowerShapes_ && iPFCluster.layer() == PFLayer::ECAL_BARREL){
             showerShapes_.clear();
             showerShapes_ = getShowerShapes(&caloBC, &(*(recHitsEB.product())), &(*topology));  
             pfCluster_e5x5.push_back(reduceFloat(showerShapes_[0],nBits_));
             pfCluster_e2x2Ratio.push_back(reduceFloat(showerShapes_[1],nBits_));
             pfCluster_e3x3Ratio.push_back(reduceFloat(showerShapes_[2],nBits_));
      	     pfCluster_eMaxRatio.push_back(reduceFloat(showerShapes_[3],nBits_));
             pfCluster_e2ndRatio.push_back(reduceFloat(showerShapes_[4],nBits_));
             pfCluster_eTopRatio.push_back(reduceFloat(showerShapes_[5],nBits_));
             pfCluster_eRightRatio.push_back(reduceFloat(showerShapes_[6],nBits_));
             pfCluster_eBottomRatio.push_back(reduceFloat(showerShapes_[7],nBits_));
             pfCluster_eLeftRatio.push_back(reduceFloat(showerShapes_[8],nBits_));
             pfCluster_e2x5MaxRatio.push_back(reduceFloat(showerShapes_[9],nBits_));
             pfCluster_e2x5TopRatio.push_back(reduceFloat(showerShapes_[10],nBits_));
             pfCluster_e2x5RightRatio.push_back(reduceFloat(showerShapes_[11],nBits_));
             pfCluster_e2x5BottomRatio.push_back(reduceFloat(showerShapes_[12],nBits_));
             pfCluster_e2x5LeftRatio.push_back(reduceFloat(showerShapes_[13],nBits_));
             pfCluster_swissCross.push_back(reduceFloat(showerShapes_[14],nBits_));
             pfCluster_r9.push_back(reduceFloat(showerShapes_[15],nBits_));
             pfCluster_sigmaIetaIeta.push_back(reduceFloat(showerShapes_[16],nBits_)); 
             pfCluster_sigmaIetaIphi.push_back(reduceFloat(showerShapes_[17],nBits_)); 
             pfCluster_sigmaIphiIphi.push_back(reduceFloat(showerShapes_[18],nBits_)); 
             pfCluster_full5x5_e5x5.push_back(reduceFloat(showerShapes_[19],nBits_));
             pfCluster_full5x5_e2x2Ratio.push_back(reduceFloat(showerShapes_[20],nBits_));
             pfCluster_full5x5_e3x3Ratio.push_back(reduceFloat(showerShapes_[21],nBits_));
             pfCluster_full5x5_eMaxRatio.push_back(reduceFloat(showerShapes_[22],nBits_));
             pfCluster_full5x5_e2ndRatio.push_back(reduceFloat(showerShapes_[23],nBits_));
             pfCluster_full5x5_eTopRatio.push_back(reduceFloat(showerShapes_[24],nBits_));
             pfCluster_full5x5_eRightRatio.push_back(reduceFloat(showerShapes_[25],nBits_));
             pfCluster_full5x5_eBottomRatio.push_back(reduceFloat(showerShapes_[26],nBits_));
             pfCluster_full5x5_eLeftRatio.push_back(reduceFloat(showerShapes_[27],nBits_));
             pfCluster_full5x5_e2x5MaxRatio.push_back(reduceFloat(showerShapes_[28],nBits_));
             pfCluster_full5x5_e2x5TopRatio.push_back(reduceFloat(showerShapes_[29],nBits_));
             pfCluster_full5x5_e2x5RightRatio.push_back(reduceFloat(showerShapes_[30],nBits_));
             pfCluster_full5x5_e2x5BottomRatio.push_back(reduceFloat(showerShapes_[31],nBits_));
             pfCluster_full5x5_e2x5LeftRatio.push_back(reduceFloat(showerShapes_[32],nBits_));
             pfCluster_full5x5_swissCross.push_back(reduceFloat(showerShapes_[33],nBits_));
             pfCluster_full5x5_r9.push_back(reduceFloat(showerShapes_[34],nBits_));
             pfCluster_full5x5_sigmaIetaIeta.push_back(reduceFloat(showerShapes_[35],nBits_)); 
             pfCluster_full5x5_sigmaIetaIphi.push_back(reduceFloat(showerShapes_[36],nBits_)); 
             pfCluster_full5x5_sigmaIphiIphi.push_back(reduceFloat(showerShapes_[37],nBits_)); 
          }else if(saveShowerShapes_ && iPFCluster.layer() == PFLayer::ECAL_ENDCAP){ 
             showerShapes_.clear();
             showerShapes_ = getShowerShapes(&caloBC, &(*(recHitsEE.product())), &(*topology));  
             pfCluster_e5x5.push_back(reduceFloat(showerShapes_[0],nBits_));
             pfCluster_e2x2Ratio.push_back(reduceFloat(showerShapes_[1],nBits_));
             pfCluster_e3x3Ratio.push_back(reduceFloat(showerShapes_[2],nBits_));
      	     pfCluster_eMaxRatio.push_back(reduceFloat(showerShapes_[3],nBits_));
             pfCluster_e2ndRatio.push_back(reduceFloat(showerShapes_[4],nBits_));
             pfCluster_eTopRatio.push_back(reduceFloat(showerShapes_[5],nBits_));
             pfCluster_eRightRatio.push_back(reduceFloat(showerShapes_[6],nBits_));
             pfCluster_eBottomRatio.push_back(reduceFloat(showerShapes_[7],nBits_));
             pfCluster_eLeftRatio.push_back(reduceFloat(showerShapes_[8],nBits_));
             pfCluster_e2x5MaxRatio.push_back(reduceFloat(showerShapes_[9],nBits_));
             pfCluster_e2x5TopRatio.push_back(reduceFloat(showerShapes_[10],nBits_));
             pfCluster_e2x5RightRatio.push_back(reduceFloat(showerShapes_[11],nBits_));
             pfCluster_e2x5BottomRatio.push_back(reduceFloat(showerShapes_[12],nBits_));
             pfCluster_e2x5LeftRatio.push_back(reduceFloat(showerShapes_[13],nBits_));
             pfCluster_swissCross.push_back(reduceFloat(showerShapes_[14],nBits_));
             pfCluster_r9.push_back(reduceFloat(showerShapes_[15],nBits_));
             pfCluster_sigmaIetaIeta.push_back(reduceFloat(showerShapes_[16],nBits_)); 
             pfCluster_sigmaIetaIphi.push_back(reduceFloat(showerShapes_[17],nBits_)); 
             pfCluster_sigmaIphiIphi.push_back(reduceFloat(showerShapes_[18],nBits_)); 
             pfCluster_full5x5_e5x5.push_back(reduceFloat(showerShapes_[19],nBits_));
             pfCluster_full5x5_e2x2Ratio.push_back(reduceFloat(showerShapes_[20],nBits_));
             pfCluster_full5x5_e3x3Ratio.push_back(reduceFloat(showerShapes_[21],nBits_));
             pfCluster_full5x5_eMaxRatio.push_back(reduceFloat(showerShapes_[22],nBits_));
             pfCluster_full5x5_e2ndRatio.push_back(reduceFloat(showerShapes_[23],nBits_));
             pfCluster_full5x5_eTopRatio.push_back(reduceFloat(showerShapes_[24],nBits_));
             pfCluster_full5x5_eRightRatio.push_back(reduceFloat(showerShapes_[25],nBits_));
             pfCluster_full5x5_eBottomRatio.push_back(reduceFloat(showerShapes_[26],nBits_));
             pfCluster_full5x5_eLeftRatio.push_back(reduceFloat(showerShapes_[27],nBits_));
             pfCluster_full5x5_e2x5MaxRatio.push_back(reduceFloat(showerShapes_[28],nBits_));
             pfCluster_full5x5_e2x5TopRatio.push_back(reduceFloat(showerShapes_[29],nBits_));
             pfCluster_full5x5_e2x5RightRatio.push_back(reduceFloat(showerShapes_[30],nBits_));
             pfCluster_full5x5_e2x5BottomRatio.push_back(reduceFloat(showerShapes_[31],nBits_));
             pfCluster_full5x5_e2x5LeftRatio.push_back(reduceFloat(showerShapes_[32],nBits_));
             pfCluster_full5x5_swissCross.push_back(reduceFloat(showerShapes_[33],nBits_));
             pfCluster_full5x5_r9.push_back(reduceFloat(showerShapes_[34],nBits_));
             pfCluster_full5x5_sigmaIetaIeta.push_back(reduceFloat(showerShapes_[35],nBits_)); 
             pfCluster_full5x5_sigmaIetaIphi.push_back(reduceFloat(showerShapes_[36],nBits_)); 
             pfCluster_full5x5_sigmaIphiIphi.push_back(reduceFloat(showerShapes_[37],nBits_)); 
          }  
          
          if(savePFClusterhits_){
             //for save PFClusterHit      
             for(unsigned int i = 0; i < hitsAndEnergies_PFCluster.at(iPFCl).size(); i++){      
                 cell = geometry->getPosition(hitsAndEnergies_PFCluster.at(iPFCl).at(i).first);
                 pfClusterHit_energy[iPFCl].push_back(reduceFloat(hitsAndEnergies_PFCluster.at(iPFCl).at(i).second,nBits_));
                 pfClusterHit_eta[iPFCl].push_back(reduceFloat(cell.eta(),nBits_));
                 pfClusterHit_phi[iPFCl].push_back(reduceFloat(cell.phi(),nBits_));
                 if(hitsAndEnergies_PFCluster.at(iPFCl).at(i).first.subdetId()==EcalBarrel){ 
                    EBDetId eb_id(hitsAndEnergies_PFCluster.at(iPFCl).at(i).first); 
                    pfClusterHit_rechitEnergy[iPFCl].push_back(reduceFloat((*(recHitsEB.product())->find(hitsAndEnergies_PFCluster.at(iPFCl).at(i).first)).energy(),nBits_)); 
                    pfClusterHit_ieta[iPFCl].push_back(eb_id.ieta());
                    pfClusterHit_iphi[iPFCl].push_back(eb_id.iphi());
                    pfClusterHit_iz[iPFCl].push_back(0); 
                 }else if(hitsAndEnergies_PFCluster.at(iPFCl).at(i).first.subdetId()==EcalEndcap){  
                    int iz=-99;
                    EEDetId ee_id(hitsAndEnergies_PFCluster.at(iPFCl).at(i).first);  
                    pfClusterHit_rechitEnergy[iPFCl].push_back(reduceFloat((*(recHitsEE.product())->find(hitsAndEnergies_PFCluster.at(iPFCl).at(i).first)).energy(),nBits_)); 
                    pfClusterHit_ieta[iPFCl].push_back(ee_id.ix());
                    pfClusterHit_iphi[iPFCl].push_back(ee_id.iy());
                    if(ee_id.zside()<0) iz=-1;
                    if(ee_id.zside()>0) iz=1;   
                    pfClusterHit_iz[iPFCl].push_back(iz); 
                 } 
             }
          }

          //compute scores     
          if(saveGenParticles_){
             for(unsigned int iGen=0; iGen<genParts.size(); iGen++){
                 if(deltaR(genParts.at(iGen).eta(),genParts.at(iGen).phi(),iPFCluster.eta(),iPFCluster.phi())<0.1) dR_genScore.push_back(deltaR(genParts.at(iGen).eta(),genParts.at(iGen).phi(),iPFCluster.eta(),iPFCluster.phi())); 
                 else dR_genScore.push_back(999.);     
             }    
             pfCluster_dR_genScore[iPFCl] = dR_genScore;        
             if(std::all_of(dR_genScore.begin(),dR_genScore.end(),[](double i){return i==-999.;}) || std::all_of(dR_genScore.begin(),dR_genScore.end(),[](double i){return i==999.;})) pfCluster_dR_genScore_MatchedIndex.push_back(-1);
             else pfCluster_dR_genScore_MatchedIndex.push_back(std::min_element(dR_genScore.begin(),dR_genScore.end()) - dR_genScore.begin()); 
          } 
          if(saveCaloParticles_){ 
             for(unsigned int iCalo=0; iCalo<caloParts.size(); iCalo++){
                 caloParticle_position = calculateAndSetPositionActual(&hitsAndEnergies_CaloPart.at(iCalo), 7.4, 3.1, 1.2, 4.2, 0.89, 0.,false);
                 if (caloParticle_position == GlobalPoint(-999999., -999999., -999999.)) {
                     std::cout << "Invalid position for caloparticle, skipping event" << std::endl;
                     return;
                 }
                 std::vector<double> scores = getScores(&hitsAndEnergies_PFCluster.at(iPFCl),&hitsAndEnergies_CaloPart.at(iCalo),recHitsEB,recHitsEE); 
                 std::vector<double> scores_1MeVCut = getScores(&hitsAndEnergies_PFCluster.at(iPFCl), &hitsAndEnergies_CaloPart_1MeVCut.at(iCalo), recHitsEB, recHitsEE);    
                 std::vector<double> scores_5MeVCut = getScores(&hitsAndEnergies_PFCluster.at(iPFCl), &hitsAndEnergies_CaloPart_5MeVCut.at(iCalo), recHitsEB, recHitsEE); 
                 std::vector<double> scores_10MeVCut = getScores(&hitsAndEnergies_PFCluster.at(iPFCl), &hitsAndEnergies_CaloPart_10MeVCut.at(iCalo), recHitsEB, recHitsEE); 
                 std::vector<double> scores_50MeVCut = getScores(&hitsAndEnergies_PFCluster.at(iPFCl), &hitsAndEnergies_CaloPart_50MeVCut.at(iCalo), recHitsEB,recHitsEE);  
                 std::vector<double> scores_100MeVCut = getScores(&hitsAndEnergies_PFCluster.at(iPFCl), &hitsAndEnergies_CaloPart_100MeVCut.at(iCalo), recHitsEB,recHitsEE);              
                 if(deltaR(caloParticle_position.eta(),caloParticle_position.phi(),iPFCluster.eta(),iPFCluster.phi())<0.1) dR_simScore.push_back(deltaR(caloParticle_position.eta(),caloParticle_position.phi(),iPFCluster.eta(),iPFCluster.phi())); 
                 else dR_simScore.push_back(999.);  
 
                 if(scoreType_=="n_shared_xtals") simScore.push_back(scores[0]);  
                 if(scoreType_=="sim_fraction") simScore.push_back(scores[5]);  
                 if(scoreType_=="simScore_final_combination") simScore.push_back(scores[1]);  
                 if(scoreType_=="sim_fraction_1MeVCut") simScore.push_back(scores[10]);  
                 if(scoreType_=="sim_fraction_5MeVCut") simScore.push_back(scores[11]);  
                 if(scoreType_=="sim_fraction_10MeVCut") simScore.push_back(scores[12]);  
                 if(scoreType_=="sim_fraction_50MeVCut") simScore.push_back(scores[13]);
                 if(scoreType_=="sim_fraction_100MeVCut") simScore.push_back(scores[14]);    
                 if(scoreType_=="sim_fraction_500MeVCut") simScore.push_back(scores[15]);    
                 if(scoreType_=="sim_fraction_1GeVCut") simScore.push_back(scores[16]);      
                 if(scoreType_=="sim_rechit_diff") simScore.push_back(scores[2]); 
                 if(scoreType_=="sim_rechit_fraction") simScore.push_back(scores[3]);           
                 if(scoreType_=="global_sim_rechit_fraction") simScore.push_back(scores[4]);
                 if(scoreType_=="hgcal_caloToCluster") simScore.push_back(scores[7]);  
                 if(scoreType_=="hgcal_clusterToCalo") simScore.push_back(scores[8]);  
                 if(scoreType_=="sim_rechit_combined_fraction") simScore.push_back(scores[17]);  
                 if(scoreType_=="rechit_sim_combined_fraction") simScore.push_back(scores[18]);  
                 
                 sim_fraction_old.push_back(scores[9]);  
                 n_shared_xtals.push_back(scores[0]);  
                 sim_fraction.push_back(scores[5]);  
                 sim_fraction_1MeVCut.push_back(scores[10]);  
                 sim_fraction_5MeVCut.push_back(scores[11]);  
                 sim_fraction_10MeVCut.push_back(scores[12]);  
                 sim_fraction_50MeVCut.push_back(scores[13]);
                 sim_fraction_100MeVCut.push_back(scores[14]);    
                 sim_fraction_500MeVCut.push_back(scores[15]);    
                 sim_fraction_1GeVCut.push_back(scores[16]);      
                 sim_rechit_diff.push_back(scores[2]); 
                 sim_rechit_fraction.push_back(scores[3]);           
                 global_sim_rechit_fraction.push_back(scores[4]);
                 hgcal_caloToCluster.push_back(scores[7]);  
                 hgcal_clusterToCalo.push_back(scores[8]);
                 sim_rechit_combined_fraction.push_back(scores[17]);  
                 rechit_sim_combined_fraction.push_back(scores[18]);        
                     
             } 
             pfCluster_nXtals.push_back((iPFCluster.hitsAndFractions()).size());   
             pfCluster_dR_simScore[iPFCl] = dR_simScore;  
             pfCluster_sim_fraction_old[iPFCl] = sim_fraction_old;        
             pfCluster_simScore[iPFCl] = simScore; 
             pfCluster_n_shared_xtals[iPFCl] = n_shared_xtals;
             pfCluster_sim_fraction[iPFCl] = sim_fraction; 
             pfCluster_sim_fraction_1MeVCut[iPFCl] = sim_fraction_1MeVCut; 
             pfCluster_sim_fraction_5MeVCut[iPFCl] = sim_fraction_5MeVCut;  
             pfCluster_sim_fraction_10MeVCut[iPFCl] = sim_fraction_10MeVCut;  
             pfCluster_sim_fraction_50MeVCut[iPFCl] = sim_fraction_50MeVCut;  
             pfCluster_sim_fraction_100MeVCut[iPFCl] = sim_fraction_100MeVCut;       
             pfCluster_sim_fraction_500MeVCut[iPFCl] = sim_fraction_500MeVCut;       
             pfCluster_sim_fraction_1GeVCut[iPFCl] = sim_fraction_1GeVCut;            
             pfCluster_sim_rechit_diff[iPFCl] = sim_rechit_diff; 
             pfCluster_sim_rechit_fraction[iPFCl] = sim_rechit_fraction;           
             pfCluster_global_sim_rechit_fraction[iPFCl] = global_sim_rechit_fraction;
             pfCluster_hgcal_caloToCluster[iPFCl] = hgcal_caloToCluster; 
             pfCluster_hgcal_clusterToCalo[iPFCl] = hgcal_clusterToCalo;  
             pfCluster_sim_rechit_combined_fraction[iPFCl] = sim_rechit_combined_fraction; 
             pfCluster_rechit_sim_combined_fraction[iPFCl] = rechit_sim_combined_fraction;     
           
             pfCluster_dR_simScore_MatchedIndex.push_back(getMatchedIndex(&pfCluster_dR_simScore, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iPFCl));
             pfCluster_sim_fraction_old_MatchedIndex.push_back(getMatchedIndex(&pfCluster_sim_fraction_old, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iPFCl));
             if(!saveScores_){
                if(scoreType_!="simScore_final_combination"){ 
                   if(scoreType_=="sim_rechit_diff" || scoreType_=="sim_rechit_fraction" || scoreType_=="global_sim_rechit_fraction" || scoreType_=="hgcal_caloToCluster" || scoreType_=="hgcal_clusterToCalo" || scoreType_=="rechit_sim_combined_fraction" || scoreType_=="sim_rechit_combined_fraction") pfCluster_simScore_MatchedIndex.push_back(getMatchedIndex(&pfCluster_simScore, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iPFCl));  
                   else pfCluster_simScore_MatchedIndex.push_back(getMatchedIndex(&pfCluster_simScore, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iPFCl));  
                }else{  
                   if(iPFCluster.layer() == PFLayer::ECAL_BARREL) pfCluster_simScore_MatchedIndex.push_back(getMatchedIndex(&pfCluster_simScore, 0.04, true, std::vector<std::vector<std::vector<double>>>({pfCluster_sim_fraction_100MeVCut}), std::vector<double>({0.01}), std::vector<std::vector<std::vector<double>>>({pfCluster_sim_fraction_old, pfCluster_global_sim_rechit_fraction}), std::vector<double>({0.8,0.5}), iPFCl));
                   else if(iPFCluster.layer() == PFLayer::ECAL_ENDCAP) pfCluster_simScore_MatchedIndex.push_back(getMatchedIndex(&pfCluster_simScore, 0.04, true, std::vector<std::vector<std::vector<double>>>({pfCluster_sim_fraction_100MeVCut}), std::vector<double>({0.01}), std::vector<std::vector<std::vector<double>>>({pfCluster_sim_fraction_old, pfCluster_global_sim_rechit_fraction}), std::vector<double>({0.1,0.5}), iPFCl));                
                } 
             }else{
                pfCluster_n_shared_xtals_MatchedIndex.push_back(getMatchedIndex(&pfCluster_n_shared_xtals, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iPFCl));     
                pfCluster_sim_fraction_MatchedIndex.push_back(getMatchedIndex(&pfCluster_sim_fraction, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iPFCl));    
                pfCluster_sim_fraction_1MeVCut_MatchedIndex.push_back(getMatchedIndex(&pfCluster_sim_fraction_1MeVCut, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iPFCl));   
                pfCluster_sim_fraction_5MeVCut_MatchedIndex.push_back(getMatchedIndex(&pfCluster_sim_fraction_5MeVCut, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iPFCl));   
                pfCluster_sim_fraction_10MeVCut_MatchedIndex.push_back(getMatchedIndex(&pfCluster_sim_fraction_10MeVCut, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iPFCl));   
                pfCluster_sim_fraction_50MeVCut_MatchedIndex.push_back(getMatchedIndex(&pfCluster_sim_fraction_50MeVCut, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iPFCl));    
                pfCluster_sim_fraction_100MeVCut_MatchedIndex.push_back(getMatchedIndex(&pfCluster_sim_fraction_100MeVCut, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iPFCl));      
                pfCluster_sim_fraction_500MeVCut_MatchedIndex.push_back(getMatchedIndex(&pfCluster_sim_fraction_500MeVCut, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iPFCl));          
                pfCluster_sim_fraction_1GeVCut_MatchedIndex.push_back(getMatchedIndex(&pfCluster_sim_fraction_1GeVCut, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iPFCl));          
                pfCluster_sim_rechit_diff_MatchedIndex.push_back(getMatchedIndex(&pfCluster_sim_rechit_diff, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iPFCl));  
                pfCluster_sim_rechit_fraction_MatchedIndex.push_back(getMatchedIndex(&pfCluster_sim_rechit_fraction, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iPFCl));  
                pfCluster_global_sim_rechit_fraction_MatchedIndex.push_back(getMatchedIndex(&pfCluster_global_sim_rechit_fraction, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iPFCl)); 
                pfCluster_hgcal_caloToCluster_MatchedIndex.push_back(getMatchedIndex(&pfCluster_hgcal_caloToCluster, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iPFCl));      
                pfCluster_hgcal_clusterToCalo_MatchedIndex.push_back(getMatchedIndex(&pfCluster_hgcal_clusterToCalo, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iPFCl));  
                pfCluster_rechit_sim_combined_fraction_MatchedIndex.push_back(getMatchedIndex(&pfCluster_rechit_sim_combined_fraction, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iPFCl));  
                pfCluster_sim_rechit_combined_fraction_MatchedIndex.push_back(getMatchedIndex(&pfCluster_sim_rechit_combined_fraction, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iPFCl));          
             }   
          }    
    
          iPFCl++;        
      } 
   }
  
   //save inverse of matchings
   if(saveCaloParticles_ && savePFCluster_){ 
      fillParticleMatchedIndex(&genParticle_pfCluster_dR_genScore_MatchedIndex,&pfCluster_dR_genScore_MatchedIndex);
   } 
   if(saveCaloParticles_ && savePFCluster_){ 
      fillParticleMatchedIndex(&caloParticle_pfCluster_dR_simScore_MatchedIndex,&pfCluster_dR_simScore_MatchedIndex);
      fillParticleMatchedIndex(&caloParticle_pfCluster_sim_fraction_old_MatchedIndex,&pfCluster_sim_fraction_old_MatchedIndex);
      if(!saveScores_){
         fillParticleMatchedIndex(&caloParticle_pfCluster_simScore_MatchedIndex,&pfCluster_simScore_MatchedIndex);
      }else{
         fillParticleMatchedIndex(&caloParticle_pfCluster_n_shared_xtals_MatchedIndex,&pfCluster_n_shared_xtals_MatchedIndex);  
         fillParticleMatchedIndex(&caloParticle_pfCluster_sim_fraction_MatchedIndex,&pfCluster_sim_fraction_MatchedIndex);  
         fillParticleMatchedIndex(&caloParticle_pfCluster_sim_fraction_1MeVCut_MatchedIndex,&pfCluster_sim_fraction_1MeVCut_MatchedIndex); 
         fillParticleMatchedIndex(&caloParticle_pfCluster_sim_fraction_5MeVCut_MatchedIndex,&pfCluster_sim_fraction_5MeVCut_MatchedIndex); 
         fillParticleMatchedIndex(&caloParticle_pfCluster_sim_fraction_10MeVCut_MatchedIndex,&pfCluster_sim_fraction_10MeVCut_MatchedIndex); 
         fillParticleMatchedIndex(&caloParticle_pfCluster_sim_fraction_50MeVCut_MatchedIndex,&pfCluster_sim_fraction_50MeVCut_MatchedIndex); 
         fillParticleMatchedIndex(&caloParticle_pfCluster_sim_fraction_100MeVCut_MatchedIndex,&pfCluster_sim_fraction_100MeVCut_MatchedIndex); 
         fillParticleMatchedIndex(&caloParticle_pfCluster_sim_fraction_500MeVCut_MatchedIndex,&pfCluster_sim_fraction_500MeVCut_MatchedIndex);  
         fillParticleMatchedIndex(&caloParticle_pfCluster_sim_fraction_1GeVCut_MatchedIndex,&pfCluster_sim_fraction_1GeVCut_MatchedIndex);      
         fillParticleMatchedIndex(&caloParticle_pfCluster_sim_rechit_diff_MatchedIndex,&pfCluster_sim_rechit_diff_MatchedIndex);     
         fillParticleMatchedIndex(&caloParticle_pfCluster_sim_rechit_fraction_MatchedIndex,&pfCluster_sim_rechit_fraction_MatchedIndex); 
         fillParticleMatchedIndex(&caloParticle_pfCluster_global_sim_rechit_fraction_MatchedIndex,&pfCluster_global_sim_rechit_fraction_MatchedIndex);   
         fillParticleMatchedIndex(&caloParticle_pfCluster_hgcal_caloToCluster_MatchedIndex,&pfCluster_hgcal_caloToCluster_MatchedIndex); 
         fillParticleMatchedIndex(&caloParticle_pfCluster_hgcal_clusterToCalo_MatchedIndex,&pfCluster_hgcal_clusterToCalo_MatchedIndex);  
         fillParticleMatchedIndex(&caloParticle_pfCluster_sim_rechit_combined_fraction_MatchedIndex,&pfCluster_sim_rechit_combined_fraction_MatchedIndex);   
         fillParticleMatchedIndex(&caloParticle_pfCluster_rechit_sim_combined_fraction_MatchedIndex,&pfCluster_rechit_sim_combined_fraction_MatchedIndex);            
      }      
   } 
   
   //Save SuperClusters 
   locCov_.clear();
   full5x5_locCov_.clear();
   if(saveSuperCluster_){
      int iSC=0;
      //std::cout << "SuperClustersEB size: " << (superClusterEB.product())->size() << std::endl;
      for(const auto& iSuperCluster : *(superClusterEB.product())){  

          superCluster_energy.push_back(reduceFloat(iSuperCluster.rawEnergy(),nBits_));
          superCluster_eta.push_back(reduceFloat(iSuperCluster.eta(),nBits_));
          superCluster_phi.push_back(reduceFloat(iSuperCluster.phi(),nBits_));
          superCluster_etaWidth.push_back(reduceFloat(iSuperCluster.etaWidth(),nBits_));
          superCluster_phiWidth.push_back(reduceFloat(iSuperCluster.phiWidth(),nBits_));
          superCluster_R.push_back(reduceFloat(iSuperCluster.position().R(),nBits_));
          superCluster_nPFClusters.push_back(iSuperCluster.clusters().size());
          math::XYZPoint caloPos = iSuperCluster.seed()->position();
          EBDetId eb_id(_ebGeom->getClosestCell(GlobalPoint(caloPos.x(),caloPos.y(),caloPos.z())));  
          superCluster_ieta.push_back(eb_id.ieta());
          superCluster_iphi.push_back(eb_id.iphi());
          superCluster_iz.push_back(0);   
 
          if(saveShowerShapes_){
             reco::CaloCluster caloBC(*iSuperCluster.seed());  
             showerShapes_.clear();
             showerShapes_ = getShowerShapes(&caloBC, &(*(recHitsEB.product())), topology);  
             superCluster_e5x5.push_back(reduceFloat(showerShapes_[0],nBits_));
             superCluster_e2x2Ratio.push_back(reduceFloat(showerShapes_[1],nBits_));
             superCluster_e3x3Ratio.push_back(reduceFloat(showerShapes_[2],nBits_));
      	     superCluster_eMaxRatio.push_back(reduceFloat(showerShapes_[3],nBits_));
             superCluster_e2ndRatio.push_back(reduceFloat(showerShapes_[4],nBits_));
             superCluster_eTopRatio.push_back(reduceFloat(showerShapes_[5],nBits_));
             superCluster_eRightRatio.push_back(reduceFloat(showerShapes_[6],nBits_));
             superCluster_eBottomRatio.push_back(reduceFloat(showerShapes_[7],nBits_));
             superCluster_eLeftRatio.push_back(reduceFloat(showerShapes_[8],nBits_));
             superCluster_e2x5MaxRatio.push_back(reduceFloat(showerShapes_[9],nBits_));
             superCluster_e2x5TopRatio.push_back(reduceFloat(showerShapes_[10],nBits_));
             superCluster_e2x5RightRatio.push_back(reduceFloat(showerShapes_[11],nBits_));
             superCluster_e2x5BottomRatio.push_back(reduceFloat(showerShapes_[12],nBits_));
             superCluster_e2x5LeftRatio.push_back(reduceFloat(showerShapes_[13],nBits_));
             superCluster_swissCross.push_back(reduceFloat(showerShapes_[14],nBits_));
             superCluster_r9.push_back(reduceFloat(showerShapes_[2]*showerShapes_[0]/iSuperCluster.rawEnergy(),nBits_));
             superCluster_sigmaIetaIeta.push_back(reduceFloat(showerShapes_[16],nBits_)); 
             superCluster_sigmaIetaIphi.push_back(reduceFloat(showerShapes_[17],nBits_)); 
             superCluster_sigmaIphiIphi.push_back(reduceFloat(showerShapes_[18],nBits_)); 
             superCluster_full5x5_e5x5.push_back(reduceFloat(showerShapes_[19],nBits_));
             superCluster_full5x5_e2x2Ratio.push_back(reduceFloat(showerShapes_[20],nBits_));
             superCluster_full5x5_e3x3Ratio.push_back(reduceFloat(showerShapes_[21],nBits_));
             superCluster_full5x5_eMaxRatio.push_back(reduceFloat(showerShapes_[22],nBits_));
             superCluster_full5x5_e2ndRatio.push_back(reduceFloat(showerShapes_[23],nBits_));
             superCluster_full5x5_eTopRatio.push_back(reduceFloat(showerShapes_[24],nBits_));
             superCluster_full5x5_eRightRatio.push_back(reduceFloat(showerShapes_[25],nBits_));
             superCluster_full5x5_eBottomRatio.push_back(reduceFloat(showerShapes_[26],nBits_));
             superCluster_full5x5_eLeftRatio.push_back(reduceFloat(showerShapes_[27],nBits_));
             superCluster_full5x5_e2x5MaxRatio.push_back(reduceFloat(showerShapes_[28],nBits_));
             superCluster_full5x5_e2x5TopRatio.push_back(reduceFloat(showerShapes_[29],nBits_));
             superCluster_full5x5_e2x5RightRatio.push_back(reduceFloat(showerShapes_[30],nBits_));
             superCluster_full5x5_e2x5BottomRatio.push_back(reduceFloat(showerShapes_[31],nBits_));
             superCluster_full5x5_e2x5LeftRatio.push_back(reduceFloat(showerShapes_[32],nBits_));
             superCluster_full5x5_swissCross.push_back(reduceFloat(showerShapes_[33],nBits_));
             superCluster_full5x5_r9.push_back(reduceFloat(showerShapes_[21]*showerShapes_[19]/iSuperCluster.rawEnergy(),nBits_));
             superCluster_full5x5_sigmaIetaIeta.push_back(reduceFloat(showerShapes_[35],nBits_)); 
             superCluster_full5x5_sigmaIetaIphi.push_back(reduceFloat(showerShapes_[36],nBits_)); 
             superCluster_full5x5_sigmaIphiIphi.push_back(reduceFloat(showerShapes_[37],nBits_)); 

             HoEs_.clear();
             HoEs_ = getHoE(&iSuperCluster, towerIso1_, towerIso2_, egammaHadTower_);
             superCluster_HoEraw.push_back(reduceFloat(HoEs_[0],nBits_)); 
             superCluster_HoErawBC.push_back(reduceFloat(HoEs_[1],nBits_)); 
          } 
         
          
          if(savePFCluster_){   
             //save clusters and superClusters mutual info
             reco::CaloCluster caloSeed(*iSuperCluster.seed());  
             for(reco::CaloCluster_iterator iBC = iSuperCluster.clustersBegin(); iBC != iSuperCluster.clustersEnd(); ++iBC){
                 reco::CaloCluster caloSCluster(*(*iBC)); 
                 int iPF=0;   
                 for(const auto& iPFCluster : *(pfClusters.product())){
                     reco::CaloCluster caloPFCluster(iPFCluster);
                     if(caloPFCluster == caloSCluster) superCluster_pfClustersIndex[iSC].push_back(iPF); 
                     if(caloPFCluster == caloSCluster && caloSCluster == caloSeed) superCluster_seedIndex[iSC]=iPF;   
                     iPF++;   
                 }     
             }      
          }
          iSC++;  
      } // for superCluseterEB

      // The global SuperCluster indexing for EE has an offset = nSuperClusterEB
      iSC = (superClusterEB.product())->size();
      int iSC_tmp=-1;
      //std::cout << "SuperClustersEE size: " << (superClusterEE.product())->size() << std::endl;

      for(const auto& iSuperCluster : *(superClusterEE.product())){    

          iSC_tmp++;
        
          superCluster_energy.push_back(reduceFloat(iSuperCluster.rawEnergy(),nBits_));
          superCluster_eta.push_back(reduceFloat(iSuperCluster.eta(),nBits_));
          superCluster_phi.push_back(reduceFloat(iSuperCluster.phi(),nBits_));
          superCluster_etaWidth.push_back(reduceFloat(iSuperCluster.etaWidth(),nBits_));
          superCluster_phiWidth.push_back(reduceFloat(iSuperCluster.phiWidth(),nBits_));
          superCluster_R.push_back(reduceFloat(iSuperCluster.position().R(),nBits_)); 
          superCluster_nPFClusters.push_back(iSuperCluster.clusters().size());  
          math::XYZPoint caloPos = iSuperCluster.seed()->position(); 
          EEDetId ee_id(_eeGeom->getClosestCell(GlobalPoint(caloPos.x(),caloPos.y(),caloPos.z())));   
          superCluster_ieta.push_back(ee_id.ix());
          superCluster_iphi.push_back(ee_id.iy());
          superCluster_iz.push_back(ee_id.zside());   

          if(saveShowerShapes_){ 
             reco::CaloCluster caloBC(*iSuperCluster.seed());  
             showerShapes_.clear();
             showerShapes_ = getShowerShapes(&caloBC, &(*(recHitsEE.product())), topology);  
             superCluster_e5x5.push_back(reduceFloat(showerShapes_[0],nBits_));
             superCluster_e2x2Ratio.push_back(reduceFloat(showerShapes_[1],nBits_));
             superCluster_e3x3Ratio.push_back(reduceFloat(showerShapes_[2],nBits_));
      	     superCluster_eMaxRatio.push_back(reduceFloat(showerShapes_[3],nBits_));
             superCluster_e2ndRatio.push_back(reduceFloat(showerShapes_[4],nBits_));
             superCluster_eTopRatio.push_back(reduceFloat(showerShapes_[5],nBits_));
             superCluster_eRightRatio.push_back(reduceFloat(showerShapes_[6],nBits_));
             superCluster_eBottomRatio.push_back(reduceFloat(showerShapes_[7],nBits_));
             superCluster_eLeftRatio.push_back(reduceFloat(showerShapes_[8],nBits_));
             superCluster_e2x5MaxRatio.push_back(reduceFloat(showerShapes_[9],nBits_));
             superCluster_e2x5TopRatio.push_back(reduceFloat(showerShapes_[10],nBits_));
             superCluster_e2x5RightRatio.push_back(reduceFloat(showerShapes_[11],nBits_));
             superCluster_e2x5BottomRatio.push_back(reduceFloat(showerShapes_[12],nBits_));
             superCluster_e2x5LeftRatio.push_back(reduceFloat(showerShapes_[13],nBits_));
             superCluster_swissCross.push_back(reduceFloat(showerShapes_[14],nBits_));
             superCluster_r9.push_back(reduceFloat(showerShapes_[2]*showerShapes_[0]/iSuperCluster.rawEnergy(),nBits_));
             superCluster_sigmaIetaIeta.push_back(reduceFloat(showerShapes_[16],nBits_)); 
             superCluster_sigmaIetaIphi.push_back(reduceFloat(showerShapes_[17],nBits_)); 
             superCluster_sigmaIphiIphi.push_back(reduceFloat(showerShapes_[18],nBits_)); 
             superCluster_full5x5_e5x5.push_back(reduceFloat(showerShapes_[19],nBits_));
             superCluster_full5x5_e2x2Ratio.push_back(reduceFloat(showerShapes_[20],nBits_));
             superCluster_full5x5_e3x3Ratio.push_back(reduceFloat(showerShapes_[21],nBits_));
             superCluster_full5x5_eMaxRatio.push_back(reduceFloat(showerShapes_[22],nBits_));
             superCluster_full5x5_e2ndRatio.push_back(reduceFloat(showerShapes_[23],nBits_));
             superCluster_full5x5_eTopRatio.push_back(reduceFloat(showerShapes_[24],nBits_));
             superCluster_full5x5_eRightRatio.push_back(reduceFloat(showerShapes_[25],nBits_));
             superCluster_full5x5_eBottomRatio.push_back(reduceFloat(showerShapes_[26],nBits_));
             superCluster_full5x5_eLeftRatio.push_back(reduceFloat(showerShapes_[27],nBits_));
             superCluster_full5x5_e2x5MaxRatio.push_back(reduceFloat(showerShapes_[28],nBits_));
             superCluster_full5x5_e2x5TopRatio.push_back(reduceFloat(showerShapes_[29],nBits_));
             superCluster_full5x5_e2x5RightRatio.push_back(reduceFloat(showerShapes_[30],nBits_));
             superCluster_full5x5_e2x5BottomRatio.push_back(reduceFloat(showerShapes_[31],nBits_));
             superCluster_full5x5_e2x5LeftRatio.push_back(reduceFloat(showerShapes_[32],nBits_));
             superCluster_full5x5_swissCross.push_back(reduceFloat(showerShapes_[33],nBits_));
             superCluster_full5x5_r9.push_back(reduceFloat(showerShapes_[21]*showerShapes_[19]/iSuperCluster.rawEnergy(),nBits_));
             superCluster_full5x5_sigmaIetaIeta.push_back(reduceFloat(showerShapes_[35],nBits_)); 
             superCluster_full5x5_sigmaIetaIphi.push_back(reduceFloat(showerShapes_[36],nBits_)); 
             superCluster_full5x5_sigmaIphiIphi.push_back(reduceFloat(showerShapes_[37],nBits_)); 

             HoEs_.clear();
             HoEs_ = getHoE(&iSuperCluster, towerIso1_, towerIso2_, egammaHadTower_);
             superCluster_HoEraw.push_back(reduceFloat(HoEs_[0],nBits_)); 
             superCluster_HoErawBC.push_back(reduceFloat(HoEs_[1],nBits_)); 
          }

          
          if(iSuperCluster.preshowerClusters().isAvailable()){
              for(unsigned int iPC=0; iPC<iSuperCluster.preshowerClusters().size(); iPC++){
                  if(!iSuperCluster.preshowerClusters()[iPC].isAvailable()) { continue; } 
                  superCluster_psCluster_energy[iSC_tmp].push_back(reduceFloat(iSuperCluster.preshowerClusters()[iPC]->energy(),nBits_));
                  superCluster_psCluster_eta[iSC_tmp].push_back(reduceFloat(iSuperCluster.preshowerClusters()[iPC]->eta(),nBits_));
                  superCluster_psCluster_phi[iSC_tmp].push_back(reduceFloat(iSuperCluster.preshowerClusters()[iPC]->phi(),nBits_));   
              }
          } 

          if(savePFCluster_){   
             //save clusters and superClusters mutual info
             reco::CaloCluster caloSeed(*iSuperCluster.seed());  
             for(reco::CaloCluster_iterator iBC = iSuperCluster.clustersBegin(); iBC != iSuperCluster.clustersEnd(); ++iBC){
                 reco::CaloCluster caloSCluster(*(*iBC));  
                 int iPF=0;   
                 for(const auto& iPFCluster : *(pfClusters.product())){
                     reco::CaloCluster caloPFCluster(iPFCluster);
                     if(caloPFCluster == caloSCluster) superCluster_pfClustersIndex[iSC].push_back(iPF); 
                     if(caloPFCluster == caloSCluster && caloSCluster == caloSeed) superCluster_seedIndex[iSC]=iPF;   
                     iPF++;   
                 }     
             }      
          }
          iSC++;  
      } // end for superCLustesEE
   }
   //save pfCluster_superClustersIndex
   if(savePFCluster_ && saveSuperCluster_){
      for(unsigned int iSC=0; iSC<superCluster_pfClustersIndex.size(); iSC++)
          for(unsigned int iPF=0; iPF<superCluster_pfClustersIndex.at(iSC).size(); iPF++)
              if(superCluster_pfClustersIndex[iSC].at(iPF)>=0) pfCluster_superClustersIndex[superCluster_pfClustersIndex[iSC].at(iPF)].push_back(iSC);   
    }
   //fill tree for each event
   tree->Fill();
}

void RecoSimDumperPF::beginJob()
{

}

void RecoSimDumperPF::endJob() 
{
    

}

///------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
std::vector<float> RecoSimDumperPF::getShowerShapes(reco::CaloCluster* caloBC, const EcalRecHitCollection* recHits, const CaloTopology *topology)
{
    std::vector<float> shapes;
    shapes.resize(38); 
    locCov_.clear();
    full5x5_locCov_.clear();
    locCov_ = EcalClusterTools::localCovariances(*caloBC, recHits, topology);
    full5x5_locCov_ = noZS::EcalClusterTools::localCovariances(*caloBC, recHits, topology);
    
    float e5x5 = EcalClusterTools::e5x5(*caloBC, recHits, topology); // e5x5
    float e3x3 = EcalClusterTools::e3x3(*caloBC, recHits, topology); // e3x3
    float eMax = EcalClusterTools::eMax(*caloBC, recHits); // eMax
    float eTop = EcalClusterTools::eTop(*caloBC, recHits, topology); // eTop 
    float eRight = EcalClusterTools::eRight(*caloBC, recHits, topology); // eRight
    float eBottom = EcalClusterTools::eBottom(*caloBC, recHits, topology); // eBottom
    float eLeft = EcalClusterTools::eLeft(*caloBC, recHits, topology); // eLeft
    float e4 = eTop + eRight + eBottom + eLeft;

    shapes[0] = e5x5;
    shapes[1] = EcalClusterTools::e2x2(*caloBC, recHits, topology)/e5x5; // e2x2/e5x5
    shapes[2] = EcalClusterTools::e3x3(*caloBC, recHits, topology)/e5x5; // e3x3/e5x5
    shapes[3] = EcalClusterTools::eMax(*caloBC,  recHits)/e5x5; // eMax/e5x5
    shapes[4] = EcalClusterTools::e2nd(*caloBC, recHits)/e5x5; // e2nd/e5x5
    shapes[5] = EcalClusterTools::eTop(*caloBC, recHits, topology)/e5x5; // eTop/e5x5 
    shapes[6] = EcalClusterTools::eRight(*caloBC, recHits, topology)/e5x5; // eRight/e5x5
    shapes[7] = EcalClusterTools::eBottom(*caloBC, recHits, topology)/e5x5; // eBottom/e5x5
    shapes[8] = EcalClusterTools::eLeft(*caloBC, recHits, topology)/e5x5; // eLeft/e5x5
    shapes[9] = EcalClusterTools::e2x5Max(*caloBC, recHits, topology)/e5x5; // e2x5Max/e5x5
    shapes[10] = EcalClusterTools::e2x5Top(*caloBC, recHits, topology)/e5x5; // e2x5Top/e5x5  
    shapes[11] = EcalClusterTools::e2x5Right(*caloBC, recHits, topology)/e5x5; // e2x5Bottom/e5x5  
    shapes[12] = EcalClusterTools::e2x5Bottom(*caloBC, recHits, topology)/e5x5; // e2x5Left/e5x5  
    shapes[13] = EcalClusterTools::e2x5Left(*caloBC, recHits, topology)/e5x5; // e2x5Right/e5x5   
    shapes[14] = 1.-e4/eMax; // swissCross 
    shapes[15] = e3x3/caloBC->energy(); // r9     
    shapes[16] = sqrt(locCov_[0]); // sigmaIetaIeta 
    shapes[17] = sqrt(locCov_[1]); // sigmaIetaIphi
    shapes[18] = !edm::isFinite(locCov_[2]) ? 0. : sqrt(locCov_[2]); // sigmaIphiIphi 

    // full_5x5 variables
    float full5x5_e5x5 = noZS::EcalClusterTools::e5x5(*caloBC, recHits, topology); // e5x5
    float full5x5_e3x3 = noZS::EcalClusterTools::e3x3(*caloBC, recHits, topology); // e3x3
    float full5x5_eMax = noZS::EcalClusterTools::eMax(*caloBC, recHits); // eMax
    float full5x5_eTop = noZS::EcalClusterTools::eTop(*caloBC, recHits, topology); // eTop 
    float full5x5_eRight = noZS::EcalClusterTools::eRight(*caloBC, recHits, topology); // eRight
    float full5x5_eBottom = noZS::EcalClusterTools::eBottom(*caloBC, recHits, topology); // eBottom
    float full5x5_eLeft = noZS::EcalClusterTools::eLeft(*caloBC, recHits, topology); // eLeft
    float full5x5_e4 = full5x5_eTop + full5x5_eRight + full5x5_eBottom + full5x5_eLeft;

    shapes[19] = full5x5_e5x5;
    shapes[20] = noZS::EcalClusterTools::e2x2(*caloBC, recHits, topology)/full5x5_e5x5; // e2x2/e5x5
    shapes[21] = noZS::EcalClusterTools::e3x3(*caloBC, recHits, topology)/full5x5_e5x5; // e3x3/e5x5
    shapes[22] = noZS::EcalClusterTools::eMax(*caloBC, recHits)/full5x5_e5x5; // eMax/e5x5
    shapes[23] = noZS::EcalClusterTools::e2nd(*caloBC, recHits)/full5x5_e5x5; // e2nd/e5x5
    shapes[24] = noZS::EcalClusterTools::eTop(*caloBC, recHits, topology)/full5x5_e5x5; // eTop/e5x5 
    shapes[25] = noZS::EcalClusterTools::eRight(*caloBC, recHits, topology)/full5x5_e5x5; // eRight/e5x5
    shapes[26] = noZS::EcalClusterTools::eBottom(*caloBC, recHits, topology)/full5x5_e5x5; // eBottom/e5x5
    shapes[27] = noZS::EcalClusterTools::eLeft(*caloBC, recHits, topology)/full5x5_e5x5; // eLeft/e5x5
    shapes[28] = noZS::EcalClusterTools::e2x5Max(*caloBC, recHits, topology)/full5x5_e5x5; // e2x5Max/e5x5
    shapes[29] = noZS::EcalClusterTools::e2x5Top(*caloBC, recHits, topology)/full5x5_e5x5; // e2x5Top/e5x5  
    shapes[30] = noZS::EcalClusterTools::e2x5Right(*caloBC, recHits, topology)/full5x5_e5x5; // e2x5Bottom/e5x5  
    shapes[31] = noZS::EcalClusterTools::e2x5Bottom(*caloBC, recHits, topology)/full5x5_e5x5; // e2x5Left/e5x5  
    shapes[32] = noZS::EcalClusterTools::e2x5Left(*caloBC, recHits, topology)/full5x5_e5x5; // e2x5Right/e5x5   
    shapes[33] = 1.-full5x5_e4/full5x5_eMax; // swissCross 
    shapes[34] = full5x5_e3x3/caloBC->energy(); // r9
    shapes[35] = sqrt(full5x5_locCov_[0]); // sigmaIetaIeta        
    shapes[36] = sqrt(full5x5_locCov_[1]); // sigmaIetaIphi          
    shapes[37] = !edm::isFinite(full5x5_locCov_[2]) ? 0. : sqrt(full5x5_locCov_[2]); // sigmaIphiIphi 

    return shapes; 
}

std::vector<float> RecoSimDumperPF::getHoE(const reco::SuperCluster* iSuperCluster, EgammaTowerIsolation* towerIso1, EgammaTowerIsolation* towerIso2, const EgammaHadTower* egammaHadTower)
{
     std::vector<float> HoEs;
     HoEs.resize(2);
  
     std::vector<CaloTowerDetId> towersBehindCluster = egammaHadTower->towersOf(*iSuperCluster);
     double HoEraw1 = towerIso1->getTowerESum(iSuperCluster)/iSuperCluster->rawEnergy();
     double HoEraw2 = towerIso2->getTowerESum(iSuperCluster)/iSuperCluster->rawEnergy();        
     float HoEraw1bc = egammaHadTower->getDepth1HcalESum(towersBehindCluster)/iSuperCluster->energy();
     float HoEraw2bc = egammaHadTower->getDepth2HcalESum(towersBehindCluster)/iSuperCluster->energy(); 
     HoEs[0] = HoEraw1 + HoEraw2;
     HoEs[1] = HoEraw1bc + HoEraw2bc;
     
     return HoEs;
}

std::vector<std::pair<DetId, float> >* RecoSimDumperPF::getHitsAndEnergiesCaloPart(CaloParticle* iCaloParticle, float simHitEnergy_cut)
{
    std::vector<std::pair<DetId, float> >* HitsAndEnergies_CaloPart_tmp = new std::vector<std::pair<DetId, float> >;
    std::vector<std::pair<DetId, float> >* HitsAndEnergies_tmp = new std::vector<std::pair<DetId, float> >;
    std::map<DetId, float> HitsAndEnergies_map;
    
    const auto& simClusters = iCaloParticle->simClusters();
    for(unsigned int iSC = 0; iSC < simClusters.size() ; iSC++){
        auto simCluster = simClusters[iSC];  
        auto hits_and_energies = simCluster->hits_and_energies();
        for(unsigned int i = 0; i < hits_and_energies.size(); i++){ 
            if(hits_and_energies[i].second < simHitEnergy_cut) continue; 
            HitsAndEnergies_tmp->push_back(make_pair(DetId(hits_and_energies[i].first),hits_and_energies[i].second));  
        }  
    }

    for(unsigned int i = 0; i < HitsAndEnergies_tmp->size(); i++){  
        if (HitsAndEnergies_map.find(HitsAndEnergies_tmp->at(i).first) == HitsAndEnergies_map.end()) {
            HitsAndEnergies_map[HitsAndEnergies_tmp->at(i).first]=HitsAndEnergies_tmp->at(i).second;      
        }else{
            HitsAndEnergies_map[HitsAndEnergies_tmp->at(i).first]=HitsAndEnergies_map[HitsAndEnergies_tmp->at(i).first]+HitsAndEnergies_tmp->at(i).second; 
        }
    }

    for(auto const& hit : HitsAndEnergies_map) 
         HitsAndEnergies_CaloPart_tmp->push_back(make_pair(hit.first,hit.second));

    return HitsAndEnergies_CaloPart_tmp;
}

std::vector<std::pair<DetId, float> >* RecoSimDumperPF::getHitsAndEnergiesBC(reco::CaloCluster* iPFCluster, edm::Handle<EcalRecHitCollection> recHitsEB, edm::Handle<EcalRecHitCollection> recHitsEE)
{
    std::vector<std::pair<DetId, float> >* HitsAndEnergies_tmp = new std::vector<std::pair<DetId, float> >;
    
    const std::vector<std::pair<DetId,float> > &hitsAndFractions = iPFCluster->hitsAndFractions();
    for(unsigned int i = 0; i < hitsAndFractions.size(); i++){
        if(hitsAndFractions.at(i).first.subdetId()==EcalBarrel){
           HitsAndEnergies_tmp->push_back(make_pair(hitsAndFractions.at(i).first,hitsAndFractions.at(i).second*(*(recHitsEB.product())->find(hitsAndFractions[i].first)).energy()));
        }else if(hitsAndFractions.at(i).first.subdetId()==EcalEndcap){
           HitsAndEnergies_tmp->push_back(make_pair(hitsAndFractions.at(i).first,hitsAndFractions.at(i).second*(*(recHitsEE.product())->find(hitsAndFractions[i].first)).energy()));
        }
    }

    return HitsAndEnergies_tmp;
}


std::vector<double> RecoSimDumperPF::getScores(const std::vector<std::pair<DetId, float> >*hits_and_energies_Cluster, const std::vector<std::pair<DetId, float> > *hits_and_energies_CaloPart, edm::Handle<EcalRecHitCollection> recHitsEB, edm::Handle<EcalRecHitCollection> recHitsEE)
{
    std::vector<double> scores;
    scores.resize(19);

    double nSharedXtals=0;
    double simFraction=0.;
    double simFraction_old=0.;
    double sim_rechit_diff=0.;
    double sim_rechit_fraction=0.;     
    double global_sim_rechit_fraction=0.;  
    double hgcal_caloToCluster=0.;      
    double hgcal_clusterToCalo=0.;   
    double sim_rechit_combined_fraction=0.;
    double rechit_sim_combined_fraction=0.;    
   
    double rechits_tot_CaloPart = 0.;
    double rechits_tot_CaloPart_noEnergy = 0.;
    for(const std::pair<DetId, float>& hit_CaloPart : *hits_and_energies_CaloPart) {
        rechits_tot_CaloPart+=hit_CaloPart.second;
        rechits_tot_CaloPart_noEnergy+=1.;
    }

    double rechits_tot_Cluster = 0.;
    double rechits_tot_Cluster_noEnergy = 0.;
    for(const std::pair<DetId, float>& hit_Cluster : *hits_and_energies_Cluster) {
        rechits_tot_Cluster+=hit_Cluster.second;
        rechits_tot_Cluster_noEnergy+=1.;
    }
   
    double rechits_match_Cluster = 0.;
    double rechits_match_CaloPart = 0.;
    double rechits_match_CaloPart_noEnergy = 0.;
    for(const std::pair<DetId, float>& hit_Cluster : *hits_and_energies_Cluster){     
        for(const std::pair<DetId, float>& hit_CaloPart : *hits_and_energies_CaloPart){  
            if(hit_CaloPart.first.rawId() == hit_Cluster.first.rawId()){

               rechits_match_Cluster += hit_Cluster.second;
               rechits_match_CaloPart += hit_CaloPart.second;    
               rechits_match_CaloPart_noEnergy += 1.0;

               sim_rechit_diff += fabs(hit_CaloPart.second-hit_Cluster.second);

               double reco_ratio=0.; 
               double sim_ratio = 0.;  
               if(rechits_tot_Cluster!=0.) reco_ratio = (double)hit_Cluster.second/(double)rechits_tot_Cluster;     
               if(rechits_tot_CaloPart!=0.) sim_ratio = (double)hit_CaloPart.second/(double)rechits_tot_CaloPart; 
               sim_rechit_fraction += fabs(sim_ratio - reco_ratio);    
            }         
        }
    }

    double hgcal_caloToCluster_Num = 0.;
    double hgcal_caloToCluster_Denum = 0.;
    double hgcal_clusterToCalo_Num = 0.;
    double hgcal_clusterToCalo_Denum = 0.;
    double reco_fraction = 0.;
    double sim_fraction = 0.;
    
    for(const std::pair<DetId, float>& hit_CaloPart : *hits_and_energies_CaloPart){  

        double rechitE=0.;
        if(hit_CaloPart.first.subdetId()==EcalBarrel) rechitE = (*(recHitsEB.product())->find(hit_CaloPart.first)).energy();
        else if(hit_CaloPart.first.subdetId()==EcalEndcap) rechitE = (*(recHitsEE.product())->find(hit_CaloPart.first)).energy(); 
        sim_fraction = (double)hit_CaloPart.second/(double)rechits_tot_CaloPart;   

        reco_fraction = 0.;   
        for(const std::pair<DetId, float>& hit_Cluster : *hits_and_energies_Cluster){   
            if(hit_CaloPart.first.rawId() == hit_Cluster.first.rawId()){
               reco_fraction = (double)rechitE/(double)rechits_tot_Cluster; 
            }
        }

        hgcal_caloToCluster_Num += (reco_fraction-sim_fraction)*(reco_fraction-sim_fraction)*rechitE*rechitE; 
        hgcal_caloToCluster_Denum += sim_fraction*sim_fraction*rechitE*rechitE;      
    }

    for(const std::pair<DetId, float>& hit_Cluster : *hits_and_energies_Cluster){  

        double rechitE=0.;
        if(hit_Cluster.first.subdetId()==EcalBarrel) rechitE = (*(recHitsEB.product())->find(hit_Cluster.first)).energy();
        else if(hit_Cluster.first.subdetId()==EcalEndcap) rechitE = (*(recHitsEE.product())->find(hit_Cluster.first)).energy(); 
        reco_fraction = (double)rechitE/(double)rechits_tot_Cluster;
        hgcal_clusterToCalo_Denum += reco_fraction*reco_fraction*rechitE*rechitE;  

        for(const std::pair<DetId, float>& hit_CaloPart : *hits_and_energies_CaloPart){
            reco_fraction = 0.;
            sim_fraction = 0.; 
            if(hit_CaloPart.first.rawId() == hit_Cluster.first.rawId()){
               reco_fraction = (double)rechitE/(double)rechits_tot_Cluster; 
               sim_fraction = (double)hit_CaloPart.second/(double)rechits_tot_CaloPart;      
               hgcal_clusterToCalo_Num += (reco_fraction-sim_fraction)*(reco_fraction-sim_fraction)*rechitE*rechitE;     
            }
        }   
    }

    nSharedXtals = (int)rechits_match_CaloPart_noEnergy;
    if(nSharedXtals==0) nSharedXtals=-999;

    if(rechits_tot_CaloPart!=0.) simFraction = (double)rechits_match_CaloPart/(double)rechits_tot_CaloPart;
    else simFraction=-999.; 
    if(simFraction==0.) simFraction=-999.;

    if(rechits_match_CaloPart!=0.) simFraction_old = fabs(1.-(double)rechits_match_Cluster/(double)rechits_match_CaloPart);
    else simFraction_old=999.; 
    if(simFraction_old==1.) simFraction_old=999.;

    if(rechits_match_CaloPart_noEnergy!=0.) sim_rechit_diff = (1./(double)rechits_match_CaloPart_noEnergy)*(double)sim_rechit_diff;
    else sim_rechit_diff=999.; 

    if(sim_rechit_fraction!=0. && rechits_match_CaloPart_noEnergy!=0.) sim_rechit_fraction = (1./(double)rechits_match_CaloPart_noEnergy)*sim_rechit_fraction;
    else sim_rechit_fraction=999.;
    
    if(rechits_tot_CaloPart!=0. && rechits_tot_Cluster!=0. && rechits_match_CaloPart!=0. && rechits_match_Cluster!=0. && rechits_match_CaloPart_noEnergy!=0.) global_sim_rechit_fraction = (1./(double)rechits_match_CaloPart_noEnergy)*fabs((double)rechits_match_CaloPart/(double)rechits_tot_CaloPart - (double)rechits_match_Cluster/(double)rechits_tot_Cluster);
    else global_sim_rechit_fraction=999.;  

    if(rechits_tot_Cluster!=0.) sim_rechit_combined_fraction = (double)rechits_match_CaloPart/(double)rechits_tot_Cluster;
    else sim_rechit_combined_fraction = 999.;
    if(rechits_tot_CaloPart!=0.) rechit_sim_combined_fraction = (double)rechits_match_Cluster/(double)rechits_tot_CaloPart;
    else rechit_sim_combined_fraction = 999.;   
 
    if(hgcal_caloToCluster_Denum!=0.) hgcal_caloToCluster = (double)hgcal_caloToCluster_Num/(double)hgcal_caloToCluster_Denum; 
    else hgcal_caloToCluster = 999.;
    if(hgcal_caloToCluster_Denum!=0.) hgcal_clusterToCalo = (double)hgcal_clusterToCalo_Num/(double)hgcal_caloToCluster_Denum;  
    else hgcal_caloToCluster = 999.;
    
    scores[0] = (double)nSharedXtals;
    scores[1] = simFraction;
    scores[2] = sim_rechit_diff;
    scores[3] = sim_rechit_fraction;     
    scores[4] = global_sim_rechit_fraction;  
    if(simFraction>0.01) scores[5] = simFraction; 
    else scores[5] = -999.; 
    if(simFraction>0.03) scores[6] = simFraction; 
    else scores[6] = -999.;
    scores[7] = hgcal_caloToCluster; 
    scores[8] = hgcal_clusterToCalo; 
    scores[9] = simFraction_old;  
    if((double)rechits_match_CaloPart>0.001) scores[10] = simFraction; 
    else scores[10] = -999.; 
    if((double)rechits_match_CaloPart>0.005) scores[11] = simFraction; 
    else scores[11] = -999.;  
    if((double)rechits_match_CaloPart>0.01) scores[12] = simFraction; 
    else scores[12] = -999.;  
    if((double)rechits_match_CaloPart>0.05) scores[13] = simFraction; 
    else scores[13] = -999.;   
    if((double)rechits_match_CaloPart>0.1) scores[14] = simFraction; 
    else scores[14] = -999.;  
    if((double)rechits_match_CaloPart>0.5) scores[15] = simFraction; 
    else scores[15] = -999.;  
    if((double)rechits_match_CaloPart>1.) scores[16] = simFraction; 
    else scores[16] = -999.;  
    scores[17] = sim_rechit_combined_fraction;
    scores[18] = rechit_sim_combined_fraction;

    return scores;
}

int RecoSimDumperPF::getMatchedIndex(std::vector<std::vector<double>>* score, double selection, bool useMax, std::vector<std::vector<std::vector<double>>> scoreSelMax, std::vector<double> selectionMax, std::vector<std::vector<std::vector<double>>> scoreSelMin, std::vector<double> selectionMin, int iCl)
{
   int matchedIndex = -1; 
   if(!useMax){ 
      std::replace(score->at(iCl).begin(),score->at(iCl).end(), -999., 999.);
      if(std::all_of(score->at(iCl).begin(),score->at(iCl).end(),[](double i){return i==-999.;}) || std::all_of(score->at(iCl).begin(),score->at(iCl).end(),[](double i){return i==999.;})) matchedIndex=-1;
      else matchedIndex = std::min_element(score->at(iCl).begin(),score->at(iCl).end()) - score->at(iCl).begin();
   }else{
      std::replace(score->at(iCl).begin(),score->at(iCl).end(), 999., -999.); 
      if(std::all_of(score->at(iCl).begin(),score->at(iCl).end(),[](double i){return i==-999.;}) || std::all_of(score->at(iCl).begin(),score->at(iCl).end(),[](double i){return i==999.;})) matchedIndex=-1;
      else matchedIndex = std::max_element(score->at(iCl).begin(),score->at(iCl).end()) - score->at(iCl).begin();
   }

   if(matchedIndex==-1) return -1;
   
   //std::cout << "getMatchedIndex - Max - iPF = " << iPF << " - " << matchedIndex << " - " << scoreSelMax.size() << " - " << selectionMax.size() << std::endl;
   //std::cout << "getMatchedIndex - Min - iPF = " << iPF << " - " << matchedIndex << " - " << scoreSelMin.size() << " - " << selectionMin.size() << std::endl;
    
   bool passSelection = true;
   for(unsigned int iSelMax=0; iSelMax < scoreSelMax.size(); iSelMax++)
       if(scoreSelMax.at(iSelMax).at(iCl).at(matchedIndex) < selectionMax.at(iSelMax)) passSelection = false;
      
   for(unsigned int iSelMin=0; iSelMin < scoreSelMin.size(); iSelMin++)
       if(scoreSelMin.at(iSelMin).at(iCl).at(matchedIndex) > selectionMin.at(iSelMin)) passSelection = false;
 
   if(useMax && score->at(iCl).at(matchedIndex) > selection && passSelection) return matchedIndex;
   if(!useMax && score->at(iCl).at(matchedIndex) < selection && passSelection) return matchedIndex; 
   
   return -1; 
}

void RecoSimDumperPF::fillParticleMatchedIndex(std::vector<std::vector<int>>* particleMatchedIndex, std::vector<int>* clusterMatchedIndex)
{
   for(unsigned int iCl=0; iCl<clusterMatchedIndex->size(); iCl++)
          if(clusterMatchedIndex->at(iCl)>=0) particleMatchedIndex->at(clusterMatchedIndex->at(iCl)).push_back(iCl);
}

GlobalPoint RecoSimDumperPF::calculateAndSetPositionActual(const std::vector<std::pair<DetId, float> > *hits_and_energies_CP, double _param_T0_EB, double _param_T0_EE, double _param_T0_ES, double _param_W0, double _param_X0, double _minAllowedNorm, bool useES)
{
  double preshowerStartEta = 1.653;
  double preshowerEndEta = 2.6;
  double cl_energy_float = 0;
  double max_e = 0.0;
  double clusterT0 = 0.0;
  DetId id_max; 
 
  // find the seed and max layer
  for(const std::pair<DetId, float>& hit_CP : *hits_and_energies_CP){
    const double rh_energyf = hit_CP.second;
    cl_energy_float += rh_energyf;
    if (rh_energyf > max_e) {
      max_e = rh_energyf;
      id_max = hit_CP.first;
    }
  }

  const CaloSubdetectorGeometry* ecal_geom = nullptr;
  // get seed geometry information
  if(id_max.subdetId()==EcalBarrel){
      ecal_geom = _ebGeom;
      clusterT0 = _param_T0_EB;
  }else if(id_max.subdetId()==EcalEndcap){  
      ecal_geom = _eeGeom;
      clusterT0 = _param_T0_EE;
  }else{
      //throw cms::Exception("InvalidLayer") << "ECAL Position Calc only accepts ECAL_BARREL or ECAL_ENDCAP";
      std::cout << "WARNING: wrong layer, ECAL Position Calc only accepts ECAL_BARREL or ECAL_ENDCAP, returning invalid position" << std::endl;
  }

  if (ecal_geom == nullptr)
     return GlobalPoint(-999999., -999999., -999999.);

  auto center_cell = ecal_geom->getGeometry(id_max);
  const double ctreta = center_cell->etaPos();
  const double actreta = std::abs(ctreta);
  // need to change T0 if in ES
  if (actreta > preshowerStartEta && actreta < preshowerEndEta && useES) {
    if (ctreta > 0 && _esPlus)
      clusterT0 = _param_T0_ES;
    if (ctreta < 0 && _esMinus)
      clusterT0 = _param_T0_ES;
  }

  // floats to reproduce exactly the EGM code
  const float maxDepth = _param_X0 * (clusterT0 + log(cl_energy_float));
  const float maxToFront = center_cell->getPosition().mag();
  // calculate the position
  const double logETot_inv = -log(cl_energy_float);
  double position_norm = 0.0;
  double x(0.0), y(0.0), z(0.0);
 
  for(const std::pair<DetId, float>& hit_CP : *hits_and_energies_CP) {
    if(hit_CP.first.subdetId()!=EcalBarrel && hit_CP.first.subdetId()!=EcalEndcap) continue;
    auto cell = ecal_geom->getGeometry(hit_CP.first);
    if(!cell.get()) continue;
    double weight = 0.0;
    const double rh_energy = hit_CP.second;
    if (rh_energy > 0.0)
      weight = std::max(0.0, (_param_W0 + log(rh_energy) + logETot_inv));
    const float depth = maxDepth + maxToFront - cell->getPosition().mag();
    const GlobalPoint pos = static_cast<const TruncatedPyramid*>(cell.get())->getPosition(depth);
    
    x += weight * pos.x();
    y += weight * pos.y();
    z += weight * pos.z();

    position_norm += weight;
  }

  // FALL BACK to LINEAR WEIGHTS
  if (position_norm == 0.) {
    for(const std::pair<DetId, float>& hit_CP : *hits_and_energies_CP) {
      if(hit_CP.first.subdetId()!=EcalBarrel && hit_CP.first.subdetId()!=EcalEndcap) continue; 
      auto cell = ecal_geom->getGeometry(hit_CP.first);
      if(!cell.get()) continue;
      double weight = 0.0;
      const double rh_energy = hit_CP.second;
      if (rh_energy > 0.0)
        weight = rh_energy / cl_energy_float;

      const float depth = maxDepth + maxToFront - cell->getPosition().mag();
      const GlobalPoint pos = cell->getPosition(depth);

      x += weight * pos.x();
      y += weight * pos.y();
      z += weight * pos.z();

      position_norm += weight;
    }
  }

  if (position_norm < _minAllowedNorm) {
    edm::LogError("WeirdClusterNormalization") << "Cluster too far from seeding cell: set position to (0,0,0).";
    return GlobalPoint(0., 0., 0.);
  } else {
    const double norm_inverse = 1.0 / position_norm;
    x *= norm_inverse;
    y *= norm_inverse;
    z *= norm_inverse;
  }

  return GlobalPoint(x, y, z); 
}

float RecoSimDumperPF::reduceFloat(float val, int bits)
{
    if(!doCompression_) return val;
    else return MiniFloatConverter::reduceMantissaToNbitsRounding(val,bits);
}


///------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(RecoSimDumperPF);

