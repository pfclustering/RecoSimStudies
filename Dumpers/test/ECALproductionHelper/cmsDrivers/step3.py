# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: step3 --conditions auto:phase1_2017_realistic -n 10 --era Run2_2017 --eventcontent RECOSIM,MINIAODSIM,DQM --runUnscheduled -s RAW2DIGI,L1Reco,RECO,RECOSIM,EI,PAT,VALIDATION:@standardValidation+@miniAODValidation,DQM:@standardDQM+@ExtraHLT+@miniAODDQM --datatier GEN-SIM-RECO,MINIAODSIM,DQMIO --geometry DB:Extended --filein file:step2.root --fileout file:step3.root
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing ('analysis')
# define the defaults here, changed from command line
options.maxEvents = -1 # -1 means all events, maxEvents considers the total over files considered
# add costum parameters
options.register ('dorecofile',
                  1, # default value
                  VarParsing.multiplicity.singleton, # singleton or list
                  VarParsing.varType.int,          # string, int, or float
                  'want to produce reco file?')
options.register ('dominiaodfile',
                  0, # default value
                  VarParsing.multiplicity.singleton, # singleton or list
                  VarParsing.varType.int,          # string, int, or float
                  'want to produce miniAOD file?')
options.register ('yearGT',
                  450, # default value
                  VarParsing.multiplicity.singleton, # singleton or list
                  VarParsing.varType.int,          # string, int, or float
                  'year on which conditions of detectors other than ECAL are based')
options.register ('thrsLumi',
                  450, # default value
                  VarParsing.multiplicity.singleton, # singleton or list
                  VarParsing.varType.int,          # string, int, or float
                  'lumi on which thrs are based')
options.register ('pfrhMultbelow2p5',
                  1.0, # default value
                  VarParsing.multiplicity.singleton, # singleton or list
                  VarParsing.varType.float,          # string, int, or float
                  'multiplier of noise used for PFRH thresholds for |eta|<2.5')
options.register ('pfrhMultabove2p5',
                  1.0, # default value
                  VarParsing.multiplicity.singleton, # singleton or list
                  VarParsing.varType.float,          # string, int, or float
                  'multiplier of noise used for PFRH thresholds for |eta|>2.5')
options.register ('seedMultbelow2p5',
                  3.0, # default value
                  VarParsing.multiplicity.singleton, # singleton or list
                  VarParsing.varType.float,          # string, int, or float
                  'multiplier of noise used for seeding threshold for |eta|<2.5')
options.register ('seedMultabove2p5',
                  3.0, # default value
                  VarParsing.multiplicity.singleton, # singleton or list
                  VarParsing.varType.float,          # string, int, or float
                  'multiplier of noise used for seeding threshold for |eta|>2.5')
options.register('doRefPfrh',
                 0,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.int,
                'Use reference values for pfrechit thresholds')
options.register('doRefSeed',
                 0,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.int,
                'Use reference values for seeding thresholds')
options.register('doRingAverageEB',
                 0,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.int,
                'use ring average for pfrh and seeding thresholds for EB')
options.register('doRingAverageEE',
                 0,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.int,
                'use ring average for pfrh and seeding thresholds for EE')
options.register('doSafetyMargin',
                 0,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.int,
                'raise value of pfrh and seeding thrs of 10%')
options.register('doPU',
                 0,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.int,
                'Use PU configurations')
options.register('nThr',
                 1,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.int,
                'Number of threads')
options.register('lumi',
                 450,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.int,
                 'lumi on which ECAL conditions are based, except for PFRH&PFSeeding')
options.register ('showerSigmaMult',
                  1.0,
                  VarParsing.multiplicity.singleton, 
                  VarParsing.varType.float,      
                  'multiplier of shower sigma')
options.register ('maxSigmaDist',
                  10.0,
                  VarParsing.multiplicity.singleton, 
                  VarParsing.varType.float,      
                  'max RH-cl distance allowed in PFClustering, in units of sigma')
options.parseArguments()
print options

import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Run3_cff import Run3
if options.doPU==0:
  process = cms.Process('RECO',Run3)
else:
  from Configuration.ProcessModifiers.premix_stage2_cff import premix_stage2 
  process = cms.Process('RECO',Run3,premix_stage2)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.RecoSim_cff')
process.load('CommonTools.ParticleFlow.EITopPAG_cff')
process.load('PhysicsTools.PatAlgos.slimming.metFilterPaths_cff')
process.load('Configuration.StandardSequences.PATMC_cff')
process.load('Configuration.StandardSequences.Validation_cff')
#process.load('DQMServices.Core.DQMStoreNonLegacy_cff')
#process.load('DQMOffline.Configuration.DQMOfflineMC_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)


# Message Logger settings
#process.MessageLogger = cms.Service('MessageLogger',
#         destinations = cms.untracked.vstring(
#              'detailedInfo',
#              'critical',
#              #'cerr'
#         ),
#         critical = cms.untracked.PSet(threshold = cms.untracked.string('ERROR')),
#         detailedInfo = cms.untracked.PSet(threshold = cms.untracked.string('INFO')),
#        #cerr = cms.untracked.PSet(threshold  = cms.untracked.string('WARNING')),
#)


# Input source
process.source = cms.Source('PoolSource',
    fileNames = cms.untracked.vstring('file:step2.root'),
    secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('step3 nevts:5000'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

if options.dorecofile == 1:
   process.RECOSIMoutput = cms.OutputModule('PoolOutputModule',
       dataset = cms.untracked.PSet(
           dataTier = cms.untracked.string('GEN-SIM-RECO'),
           filterName = cms.untracked.string('')
       ),
       fileName = cms.untracked.string('file:step3.root'),
       outputCommands = process.RECOSIMEventContent.outputCommands,
       splitLevel = cms.untracked.int32(0)
   )
   process.RECOSIMoutput.outputCommands.extend(['keep *_mix_MergedCaloTruth_*',
                                               #'keep *PCaloHit*_g4SimHits_EcalHitsE*_*',
                                                'keep *_particleFlowRecHitECAL_*_*',
                                                'keep *_*towerMaker*_*_*'])

if options.dominiaodfile == 1:
   process.MINIAODSIMoutput = cms.OutputModule('PoolOutputModule',
       compressionAlgorithm = cms.untracked.string('LZMA'),
       compressionLevel = cms.untracked.int32(4),
       dataset = cms.untracked.PSet(
           dataTier = cms.untracked.string('MINIAODSIM'),
           filterName = cms.untracked.string('')
       ),
       dropMetaData = cms.untracked.string('ALL'),
       eventAutoFlushCompressedSize = cms.untracked.int32(-900),
       fastCloning = cms.untracked.bool(False),
       fileName = cms.untracked.string('file:step3.root'),
       outputCommands = process.MINIAODSIMEventContent.outputCommands,
       overrideBranchesSplitLevel = cms.untracked.VPSet(
           cms.untracked.PSet(
               branch = cms.untracked.string('patPackedCandidates_packedPFCandidates__*'),
               splitLevel = cms.untracked.int32(99)
           ), 
           cms.untracked.PSet(
               branch = cms.untracked.string('recoGenParticles_prunedGenParticles__*'),
               splitLevel = cms.untracked.int32(99)
           ), 
           cms.untracked.PSet(
               branch = cms.untracked.string('patTriggerObjectStandAlones_slimmedPatTrigger__*'),
               splitLevel = cms.untracked.int32(99)
           ), 
           cms.untracked.PSet(
               branch = cms.untracked.string('patPackedGenParticles_packedGenParticles__*'),
               splitLevel = cms.untracked.int32(99)
           ), 
           cms.untracked.PSet(
               branch = cms.untracked.string('patJets_slimmedJets__*'),
               splitLevel = cms.untracked.int32(99)
           ), 
           cms.untracked.PSet(
               branch = cms.untracked.string('recoVertexs_offlineSlimmedPrimaryVertices__*'),
               splitLevel = cms.untracked.int32(99)
           ), 
           cms.untracked.PSet(
               branch = cms.untracked.string('recoCaloClusters_reducedEgamma_reducedESClusters_*'),
               splitLevel = cms.untracked.int32(99)
           ), 
           cms.untracked.PSet(
               branch = cms.untracked.string('EcalRecHitsSorted_reducedEgamma_reducedEBRecHits_*'),
               splitLevel = cms.untracked.int32(99)
           ), 
           cms.untracked.PSet(
               branch = cms.untracked.string('EcalRecHitsSorted_reducedEgamma_reducedEERecHits_*'),
               splitLevel = cms.untracked.int32(99)
           ), 
           cms.untracked.PSet(
               branch = cms.untracked.string('recoGenJets_slimmedGenJets__*'),
               splitLevel = cms.untracked.int32(99)
           ), 
           cms.untracked.PSet(
               branch = cms.untracked.string('patJets_slimmedJetsPuppi__*'),
               splitLevel = cms.untracked.int32(99)
           ), 
           cms.untracked.PSet(
               branch = cms.untracked.string('EcalRecHitsSorted_reducedEgamma_reducedESRecHits_*'),
               splitLevel = cms.untracked.int32(99)
           )
       ),
       overrideInputFileSplitLevels = cms.untracked.bool(True),
       splitLevel = cms.untracked.int32(0)
   )

   
# Additional output definition

# Other statements
process.mix.playback = True
process.mix.digitizers = cms.PSet()
for a in process.aliases: delattr(process, a)
process.RandomNumberGeneratorService.restoreStateLabel=cms.untracked.string('randomEngineStateProducer')
from Configuration.AlCa.GlobalTag import GlobalTag
if options.yearGT == 2021:
  process.GlobalTag = GlobalTag(process.GlobalTag, '106X_mcRun3_2021_realistic_v3', '')
elif options.yearGT == 2023:
  process.GlobalTag = GlobalTag(process.GlobalTag, '106X_mcRun3_2023_realistic_v3', '')
else:
  raise RuntimeError('Global tag not set properly, check logic')

# Override ECAL tags
print 'Will override following ECAL tags'
process.GlobalTag.toGet = cms.VPSet()
from override_ECAL_tags import override_tags
for rec,tag in override_tags[options.lumi].items():
  process.GlobalTag.toGet.append( cms.PSet(record = cms.string(rec), tag = cms.string(tag) )   )
  print rec,tag
  #print process.GlobalTag.toGet[0]

# override a global tag with the conditions from external module
from CalibCalorimetry.EcalTrivialCondModules.EcalTrivialCondRetriever_cfi import *
process.myCond = EcalTrivialConditionRetriever.clone()
# prefer these conditions over the globalTag's ones
process.es_prefer = cms.ESPrefer('EcalTrivialConditionRetriever','myCond')

# choose file from which to load the thresholds
if options.doRefPfrh == 0:
  # load EB thrs file
  if options.doRingAverageEB == 0:
    EB_pfrhthr_file = cms.untracked.string('./data/noise/PFRecHitThresholds_EB_TL{l}.txt'.format(l=options.thrsLumi))
  else:
    EB_pfrhthr_file = cms.untracked.string('./data/noise/PFRecHitThresholds_EB_ringaveraged_TL{l}.txt'.format(l=options.thrsLumi))
  # load EE thrs file
  if options.doRingAverageEE == 0:
    EE_pfrhthr_file = cms.untracked.string('./data/noise/PFRecHitThresholds_EE_TL{l}.txt'.format(l=options.thrsLumi))
  else:
    EE_pfrhthr_file = cms.untracked.string('./data/noise/PFRecHitThresholds_EE_ringaveraged_TL{l}.txt'.format(l=options.thrsLumi))
else: #no ringAvg for reference thrs
  EB_pfrhthr_file = cms.untracked.string('./data/noise/PFRecHitThresholds_EB_TL{l}.txt'.format(l=options.thrsLumi))
  EE_pfrhthr_file = cms.untracked.string('./data/noise/PFRecHitThresholds_EE_TL{l}.txt'.format(l=options.thrsLumi))


### set pfrh thresholds
process.myCond.producedEcalPFRecHitThresholds = cms.untracked.bool(True)
process.myCond.PFRecHitFile =   EB_pfrhthr_file 
process.myCond.PFRecHitFileEE = EE_pfrhthr_file

if options.doRefPfrh == 0:
  if options.doSafetyMargin == 1:   
    process.myCond.EcalPFRecHitThresholdNSigmas = cms.untracked.double(options.pfrhMultbelow2p5/2.0)
    process.myCond.EcalPFRecHitThresholdNSigmasHEta = cms.untracked.double(options.pfrhMultabove2p5/3.0)
  else:
    process.myCond.EcalPFRecHitThresholdNSigmas = cms.untracked.double(options.pfrhMultbelow2p5/2.2)
    process.myCond.EcalPFRecHitThresholdNSigmasHEta = cms.untracked.double(options.pfrhMultabove2p5/3.3)
else: # use the reference values
  process.myCond.EcalPFRecHitThresholdNSigmas = cms.untracked.double(1.0)
  process.myCond.EcalPFRecHitThresholdNSigmasHEta = cms.untracked.double(1.0)

### set seeding thresholds
process.myCond.producedEcalPFSeedingThresholds = cms.untracked.bool(True)
if options.doRefSeed == 0:
  if options.doSafetyMargin == 1:
    process.myCond.EcalPFSeedingThresholdNSigmas = cms.untracked.double(options.seedMultbelow2p5/2.0) # PFRHs files are at 2sigma of the noise for |eta|<2.5
    process.myCond.EcalPFSeedingThresholdNSigmasHEta = cms.untracked.double(options.seedMultabove2p5/3.0) #                3sigma of the noise for |eta|>2.5
  else:
    process.myCond.EcalPFSeedingThresholdNSigmas = cms.untracked.double(options.seedMultbelow2p5/2.2) # PFRHs files are at 2sigma of the noise for |eta|<2.5
    process.myCond.EcalPFSeedingThresholdNSigmasHEta = cms.untracked.double(options.seedMultabove2p5/3.3) #                3sigma of the noise for |eta|>2.5
  process.myCond.PFSeedingFile =   EB_pfrhthr_file
  process.myCond.PFSeedingFileEE = EE_pfrhthr_file
else: # use the reference values
  process.myCond.EcalPFSeedingThresholdNSigmas = cms.untracked.double(1.0) 
  process.myCond.EcalPFSeedingThresholdNSigmasHEta = cms.untracked.double(1.0) 
  process.myCond.PFSeedingFile =   cms.untracked.string('./data/noise/fixed_SeedingThresholds_EB.txt')
  process.myCond.PFSeedingFileEE = cms.untracked.string('./data/noise/fixed_SeedingThresholds_EE.txt')

### explicitly set to false all other ecal conditions
process.myCond.producedEcalPedestals = cms.untracked.bool(False)
process.myCond.producedEcalWeights = cms.untracked.bool(False)
process.myCond.producedEcalGainRatios = cms.untracked.bool(False)
process.myCond.producedEcalADCToGeVConstant = cms.untracked.bool(False)
process.myCond.producedEcalMappingElectronics = cms.untracked.bool(False)
process.myCond.producedEcalTimeOffsetConstant = cms.untracked.bool(False)
process.myCond.producedEcalLinearCorrections = cms.untracked.bool(False)
process.myCond.producedEcalIntercalibConstants = cms.untracked.bool(False)
process.myCond.producedEcalIntercalibConstantsMC = cms.untracked.bool(False)
process.myCond.producedEcalIntercalibErrors = cms.untracked.bool(False)
process.myCond.producedEcalTimeCalibConstants = cms.untracked.bool(False)
process.myCond.producedEcalTimeCalibErrors = cms.untracked.bool(False)
process.myCond.producedEcalSimPulseShape = cms.untracked.bool(False)
process.myCond.producedEcalChannelStatus = cms.untracked.bool(False)
process.myCond.producedEcalDQMChannelStatus = cms.untracked.bool(False)
process.myCond.producedEcalDCSTowerStatus = cms.untracked.bool(False)
process.myCond.producedEcalDAQTowerStatus = cms.untracked.bool(False)
process.myCond.producedEcalDQMTowerStatus = cms.untracked.bool(False)
process.myCond.producedEcalTrgChannelStatus = cms.untracked.bool(False)
process.myCond.producedEcalAlignmentEB = cms.untracked.bool(False)
process.myCond.producedEcalAlignmentEE = cms.untracked.bool(False)
process.myCond.producedEcalAlignmentEE = cms.untracked.bool(False)
process.myCond.producedEcalSampleMask = cms.untracked.bool(False)
# additional booleans to excplicitly set to False
process.myCond.producedEcalTimeBiasCorrections = cms.untracked.bool(False)
process.myCond.producedEcalSamplesCorrelation = cms.untracked.bool(False)
process.myCond.producedEcalLaserCorrection = cms.untracked.bool(False)
process.myCond.producedEcalClusterLocalContCorrParameters = cms.untracked.bool(False)
process.myCond.producedEcalClusterCrackCorrParameters = cms.untracked.bool(False)
process.myCond.producedEcalClusterEnergyCorrectionParameters = cms.untracked.bool(False)
process.myCond.producedEcalClusterEnergyUncertaintyParameters = cms.untracked.bool(False)
process.myCond.producedEcalClusterEnergyCorrectionObjectSpecificParameters = cms.untracked.bool(False)

########## end override

# Path and EndPath definitions 
process.raw2digi_step = cms.Path(process.RawToDigi)
process.L1Reco_step = cms.Path(process.L1Reco)
process.reconstruction_step = cms.Path(process.reconstruction)
process.recosim_step = cms.Path(process.recosim)
process.eventinterpretaion_step = cms.Path(process.EIsequence)
if options.dorecofile == 1:
   process.RECOSIMoutput_step = cms.EndPath(process.RECOSIMoutput)
if options.dominiaodfile == 1:
   process.MINIAODSIMoutput_step = cms.EndPath(process.MINIAODSIMoutput)

if options.dorecofile == 1:   
   process.schedule = cms.Schedule(process.raw2digi_step,process.L1Reco_step,process.reconstruction_step,process.recosim_step,
                               process.eventinterpretaion_step,  
                               process.RECOSIMoutput_step)
if options.dominiaodfile == 1:
   process.schedule = cms.Schedule(process.raw2digi_step,process.L1Reco_step,process.reconstruction_step,process.recosim_step,
                               process.eventinterpretaion_step,  
                               process.MINIAODSIMoutput_step)

process.schedule.associate(process.patTask)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

#Setup FWK for multithreaded 
process.options.numberOfThreads=cms.untracked.uint32(options.nThr)
process.options.numberOfStreams=cms.untracked.uint32(0)
process.options.numberOfConcurrentLuminosityBlocks=cms.untracked.uint32(1)


# customisation of the process.
if options.doRefPfrh == 0:
  process.particleFlowClusterECALUncorrected.initialClusteringStep.thresholdsByDetector = cms.VPSet(
        cms.PSet( detector = cms.string('ECAL_BARREL'),
               gatheringThreshold = cms.double(0.0),
               gatheringThresholdPt = cms.double(0.0)
               ),
        cms.PSet( detector = cms.string('ECAL_ENDCAP'),
               gatheringThreshold = cms.double(0.0),
               gatheringThresholdPt = cms.double(0.0)
               )
  )
else: # use reference values
  process.particleFlowClusterECALUncorrected.initialClusteringStep.thresholdsByDetector = cms.VPSet(
        cms.PSet( detector = cms.string('ECAL_BARREL'),
               gatheringThreshold = cms.double(0.08),
               gatheringThresholdPt = cms.double(0.0)
               ),
        cms.PSet( detector = cms.string('ECAL_ENDCAP'),
               gatheringThreshold = cms.double(0.30),
               gatheringThresholdPt = cms.double(0.0)
               )
  )

if options.showerSigmaMult != 1:
  process.particleFlowClusterECALUncorrected.pfClusterBuilder.showerSigma = cms.double(1.5*options.showerSigmaMult)

###MG commented because needs modifications in cmssw 
###process.particleFlowClusterECALUncorrected.pfClusterBuilder.maxSigmaDist = cms.double(options.maxSigmaDist)

# Automatic addition of the customisation function from SimGeneral.MixingModule.fullMixCustomize_cff
from SimGeneral.MixingModule.fullMixCustomize_cff import setCrossingFrameOn 

#call to customisation function setCrossingFrameOn imported from SimGeneral.MixingModule.fullMixCustomize_cff
process = setCrossingFrameOn(process)

# End of customisation functions
#do not add changes to your config after this point (unless you know what you are doing)
from FWCore.ParameterSet.Utilities import convertToUnscheduled
process=convertToUnscheduled(process)

# customisation of the process.

# Automatic addition of the customisation function from PhysicsTools.PatAlgos.slimming.miniAOD_tools
from PhysicsTools.PatAlgos.slimming.miniAOD_tools import miniAOD_customizeAllMC 

#call to customisation function miniAOD_customizeAllMC imported from PhysicsTools.PatAlgos.slimming.miniAOD_tools
process = miniAOD_customizeAllMC(process)

# End of customisation functions

# Customisation from command line
process.PixelCPEGenericESProducer.IrradiationBiasCorrection = True

#Have logErrorHarvester wait for the same EDProducers to finish as those providing data for the OutputModule
from FWCore.Modules.logErrorHarvester_cff import customiseLogErrorHarvesterUsingOutputCommands
process = customiseLogErrorHarvesterUsingOutputCommands(process)

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
