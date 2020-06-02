# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: SingleGammaPt35_pythia8_cfi --conditions auto:phase1_2017_realistic -n 10 --era Run2_2017 --eventcontent FEVTDEBUG --relval 9000,50 -s GEN,SIM --datatier GEN-SIM --beamspot Realistic25ns13TeVEarly2017Collision --geometry DB:Extended --fileout file:step1.root
import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing('standard')
'''
options.register('maxEvents',
                10,
                VarParsing.multiplicity.singleton,
                VarParsing.varType.int,
                "Number of events")
options.register('etmin',
                 0,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.float,
                "ET or E min depending on doFlatEnergy")
options.register('etmax',
                 0,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.float,
                "ET or E max depending on doFlatEnergy")
options.register('rmin',
                 0,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.float,
                "Radius min")
options.register('rmax',
                 0,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.float,
                "Radius max")
options.register('zmin',
                 0,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.float,
                "Z min")
options.register('zmax',
                 0,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.float,
                "Zmax")
options.register('np',
                 1,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.int,
                "Number of Particles")
'''
options.register('nThr',
                 1,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.int,
                "Number of threads")
options.register('seedOffset',
                 1,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.int,
                 "Seed offset")
'''
options.register('doFlatEnergy',
                 0,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.int,
                 "generate flat in energy, otherwise flat in pt")
'''
options.register ("yearGT",
                  450, # default value
                  VarParsing.multiplicity.singleton, # singleton or list
                  VarParsing.varType.int,          # string, int, or float
                  "year on which conditions of detectors other than ECAL are based")
options.register('lumi',
                 450,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.int,
                 "lumi on which ECAL conditions are based, except for PFRH&PFSeeding")

options.parseArguments()
print options

#my_doFlatEnergy = cms.bool(False) if options.doFlatEnergy == 0 else cms.bool(True)

from Configuration.Eras.Era_Run3_cff import Run3
process = cms.Process('SIM',Run3)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.GeometrySimDB_cff')  #UL2016 driver: process.load('Configuration.Geometry.GeometrySimDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedRun3RoundOptics25ns13TeVLowSigmaZ_cfi') #UL2016 driver: process.load('IOMC.EventVertexGenerators.VtxSmearedRealistic25ns13TeV2016Collision_cfi') 
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')


process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)

# Input source
process.source = cms.Source("EmptySource")

process.options = cms.untracked.PSet()

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('SingleGammaPt35_pythia8_cfi nevts:10'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

# in UL2016 driver: process.FEVTDEBUGoutput
process.RAWSIMoutput = cms.OutputModule("PoolOutputModule",
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('generation_step')
    ),
    # we add compression compared to UL2016 driver
    compressionAlgorithm = cms.untracked.string('LZMA'),
    compressionLevel = cms.untracked.int32(1),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM'),
        filterName = cms.untracked.string('')
    ),
    # compression option added here as well
    eventAutoFlushCompressedSize = cms.untracked.int32(20971520),
    fileName = cms.untracked.string('file:step1.root'),
    outputCommands = process.RAWSIMEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

# Additional output definition

# Other statements
# added compared to UL2016 driver:
process.XMLFromDBSource.label = cms.string("Extended")

process.genstepfilter.triggerConditions=cms.vstring("generation_step")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '106X_upgrade2021_realistic_v5', '')

# Override ECAL tags
print 'Will override following ECAL tags'
process.GlobalTag.toGet = cms.VPSet()
from override_ECAL_tags import override_tags
for rec,tag in override_tags[options.lumi].items():
  process.GlobalTag.toGet.append( cms.PSet(record = cms.string(rec), tag = cms.string(tag) )   )
  print rec,tag
  #print process.GlobalTag.toGet[0]



process.generator = cms.EDFilter("Pythia8GeneratorFilter",
    PythiaParameters = cms.PSet(
        parameterSets = cms.vstring(
            'pythia8CommonSettings', 
            'pythia8CP5Settings', 
            'processParameters'
        ),
        processParameters = cms.vstring(
            'HardQCD:all = on', 
            'PhaseSpace:pTHatMin = 15.', 
            'PhaseSpace:pTHatMax = 3000.'
        ),
        pythia8CP5Settings = cms.vstring(
            'Tune:pp 14', 
            'Tune:ee 7', 
            'MultipartonInteractions:ecmPow=0.03344', 
            'MultipartonInteractions:bProfile=2', 
            'MultipartonInteractions:pT0Ref=1.41', 
            'MultipartonInteractions:coreRadius=0.7634', 
            'MultipartonInteractions:coreFraction=0.63', 
            'ColourReconnection:range=5.176', 
            'SigmaTotal:zeroAXB=off', 
            'SpaceShower:alphaSorder=2', 
            'SpaceShower:alphaSvalue=0.118', 
            'SigmaProcess:alphaSvalue=0.118', 
            'SigmaProcess:alphaSorder=2', 
            'MultipartonInteractions:alphaSvalue=0.118', 
            'MultipartonInteractions:alphaSorder=2', 
            'TimeShower:alphaSorder=2', 
            'TimeShower:alphaSvalue=0.118', 
            'SigmaTotal:mode = 0', 
            'SigmaTotal:sigmaEl = 21.89', 
            'SigmaTotal:sigmaTot = 100.309', 
            'PDF:pSet=LHAPDF6:NNPDF31_nnlo_as_0118'
        ),
        pythia8CommonSettings = cms.vstring(
            'Tune:preferLHAPDF = 2', 
            'Main:timesAllowErrors = 10000', 
            'Check:epTolErr = 0.01', 
            'Beams:setProductionScalesFromLHEF = off', 
            'SLHA:keepSM = on', 
            'SLHA:minMassSM = 1000.', 
            'ParticleDecays:limitTau0 = on', 
            'ParticleDecays:tau0Max = 10', 
            'ParticleDecays:allowPhotonRadiation = on'
        )
    ),
    comEnergy = cms.double(13000.0),
    crossSection = cms.untracked.double(74310000.0),
    filterEfficiency = cms.untracked.double(0.037),
    maxEventsToPrint = cms.untracked.int32(0),
    pythiaHepMCVerbosity = cms.untracked.bool(False),
    pythiaPylistVerbosity = cms.untracked.int32(1),
    reweightGen = cms.PSet(
        pTRef = cms.double(15.0),
        power = cms.double(4.5)
    )
)


# Path and EndPath definitions
process.generation_step = cms.Path(process.pgen)
process.simulation_step = cms.Path(process.psim)
process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.RAWSIMoutput_step = cms.EndPath(process.RAWSIMoutput) # in UL2016 driver: process.FEVTDEBUGoutput_step = cms.EndPath(process.FEVTDEBUGoutput)

# Schedule definition
process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step,process.simulation_step,process.endjob_step,process.RAWSIMoutput_step) #UL2016 driver: process.FEVTDEBUGoutput_step
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

# filter all path with the production filter sequence
for path in process.paths:
	getattr(process,path).insert(0, process.generator)

#Setup FWK for multithreaded
process.options.numberOfThreads=cms.untracked.uint32(options.nThr)
process.options.numberOfStreams=cms.untracked.uint32(0)
process.options.numberOfConcurrentLuminosityBlocks=cms.untracked.uint32(1)

# set a different offset seed, if you run multiple jobs 
process.RandomNumberGeneratorService.eventSeedOffset=cms.untracked.uint32(options.seedOffset)

# Customisation from command line

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
