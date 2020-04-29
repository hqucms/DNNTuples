import FWCore.ParameterSet.Config as cms

# ---------------------------------------------------------

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing('analysis')

options.outputFile = 'output.root'
options.inputFiles = '/store/mc/RunIIAutumn18MiniAOD/BulkGravitonToHHTo4Q_MX-600to6000_MH-15to250_part1_TuneCP5_13TeV-madgraph_pythia8/MINIAODSIM/multigridpack_102X_upgrade2018_realistic_v15-v1/280000/F3DA5412-2411-D548-9C93-653B4840B855.root'
options.maxEvents = -1

options.register('skipEvents', 0, VarParsing.multiplicity.singleton, VarParsing.varType.int, "skip N events")
options.register('job', 0, VarParsing.multiplicity.singleton, VarParsing.varType.int, "job number")
options.register('nJobs', 1, VarParsing.multiplicity.singleton, VarParsing.varType.int, "total jobs")
options.register('inputDataset',
                 '',
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 "Input dataset")
options.register('isTrainSample', True, VarParsing.multiplicity.singleton, VarParsing.varType.bool, "if the sample is used for training")

options.setupTags(tag='%d', ifCond='nJobs > 1', tagArg='job')
options.parseArguments()

# ---------------------------------------------------------

process = cms.Process("DNNFiller")

process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.options = cms.untracked.PSet(
    allowUnscheduled = cms.untracked.bool(True),
    wantSummary=cms.untracked.bool(False)
)

print ('Using output file ' + options.outputFile)

process.TFileService = cms.Service("TFileService",
                                   fileName=cms.string(options.outputFile))

process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(options.maxEvents))

process.source = cms.Source('PoolSource',
    fileNames=cms.untracked.vstring(options.inputFiles),
    skipEvents=cms.untracked.uint32(options.skipEvents)
)


numberOfFiles = len(process.source.fileNames)
numberOfJobs = options.nJobs
jobNumber = options.job

process.source.fileNames = process.source.fileNames[jobNumber:numberOfFiles:numberOfJobs]
if options.nJobs > 1:
    print ("running over these files:")
    print (process.source.fileNames)

# ---------------------------------------------------------

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.EventContent.EventContent_cff")
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '102X_upgrade2018_realistic_v21', '')
print 'Using global tag', process.GlobalTag.globaltag

# ---------------------------------------------------------
from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection
from RecoBTag.ONNXRuntime.pfDeepBoostedJet_cff import _pfDeepBoostedJetTagsAll as pfDeepBoostedJetTagsAll

useReclusteredJets = True
jetR = 0.8

bTagInfos = [
    'pfBoostedDoubleSVAK8TagInfos'
]

bTagDiscriminators = [
    'pfCombinedInclusiveSecondaryVertexV2BJetTags',
    'pfBoostedDoubleSecondaryVertexAK8BJetTags',
    'pfDeepDoubleBvLJetTags:probHbb',
    'pfDeepDoubleCvLJetTags:probHcc',
    'pfDeepDoubleCvBJetTags:probHcc',
    'pfMassIndependentDeepDoubleBvLJetTags:probHbb',
    'pfMassIndependentDeepDoubleCvLJetTags:probHcc',
    'pfMassIndependentDeepDoubleCvBJetTags:probHcc',
]

subjetBTagDiscriminators = ['None']

if useReclusteredJets:
    JETCorrLevels = ['L2Relative', 'L3Absolute']

    from DeepNTuples.Ntupler.jetToolbox_cff import jetToolbox
    jetToolbox(process, 'ak8', 'dummySeq', 'noOutput',
               PUMethod='Puppi', JETCorrPayload='AK8PFPuppi', JETCorrLevels=JETCorrLevels,
               Cut='pt > 170.0 && abs(rapidity()) < 2.4',
               runOnMC=True,
               addNsub=True, maxTau=3,
               addSoftDrop=True, addSoftDropSubjets=True, subJETCorrPayload='AK4PFPuppi', subJETCorrLevels=JETCorrLevels,
               bTagDiscriminators=['pfCombinedInclusiveSecondaryVertexV2BJetTags'], subjetBTagDiscriminators=subjetBTagDiscriminators)

    updateJetCollection(
       process,
       jetSource=cms.InputTag('packedPatJetsAK8PFPuppiSoftDrop'),
       rParam=jetR,
       jetCorrections=('AK8PFPuppi', cms.vstring(['L2Relative', 'L3Absolute']), 'None'),
       btagDiscriminators=bTagDiscriminators + pfDeepBoostedJetTagsAll,
       btagInfos=bTagInfos,
       postfix='AK8WithPuppiDaughters',  # needed to tell the producers that the daughters are puppi-weighted
    )
    process.updatedPatJetsTransientCorrectedAK8WithPuppiDaughters.addTagInfos = cms.bool(True)
    process.updatedPatJetsTransientCorrectedAK8WithPuppiDaughters.addBTagInfo = cms.bool(True)

    srcJets = cms.InputTag('selectedUpdatedPatJetsAK8WithPuppiDaughters')
else:
    updateJetCollection(
       process,
       jetSource=cms.InputTag('slimmedJetsAK8'),
       rParam=jetR,
       jetCorrections=('AK8PFPuppi', cms.vstring(['L2Relative', 'L3Absolute']), 'None'),
       btagDiscriminators=bTagDiscriminators + pfDeepBoostedJetTagsAll,
       btagInfos=bTagInfos,
    )
    process.updatedPatJetsTransientCorrected.addTagInfos = cms.bool(True)
    process.updatedPatJetsTransientCorrected.addBTagInfo = cms.bool(True)

    srcJets = cms.InputTag('selectedUpdatedPatJets')
# ---------------------------------------------------------
from PhysicsTools.PatAlgos.tools.helpers import getPatAlgosToolsTask, addToProcessAndTask
patTask = getPatAlgosToolsTask(process)

from RecoJets.JetProducers.ak8GenJets_cfi import ak8GenJets
process.ak8GenJetsWithNu = ak8GenJets.clone(
    src='packedGenParticles',
    rParam=cms.double(jetR),
    jetPtMin=100.0
    )
process.ak8GenJetsWithNuSoftDrop = process.ak8GenJetsWithNu.clone(
    useSoftDrop=cms.bool(True),
    zcut=cms.double(0.1),
    beta=cms.double(0.0),
    R0=cms.double(jetR),
    useExplicitGhosts=cms.bool(True)
    )
process.ak8GenJetsWithNuMatch = cms.EDProducer("GenJetMatcher",  # cut on deltaR; pick best by deltaR
    src=srcJets,  # RECO jets (any View<Jet> is ok)
    matched=cms.InputTag("ak8GenJetsWithNu"),  # GEN jets  (must be GenJetCollection)
    mcPdgId=cms.vint32(),  # n/a
    mcStatus=cms.vint32(),  # n/a
    checkCharge=cms.bool(False),  # n/a
    maxDeltaR=cms.double(jetR),  # Minimum deltaR for the match
    # maxDPtRel   = cms.double(3.0),                  # Minimum deltaPt/Pt for the match (not used in GenJetMatcher)
    resolveAmbiguities=cms.bool(True),  # Forbid two RECO objects to match to the same GEN object
    resolveByMatchQuality=cms.bool(False),  # False = just match input in order; True = pick lowest deltaR pair first
)
process.ak8GenJetsWithNuSoftDropMatch = cms.EDProducer("GenJetMatcher",  # cut on deltaR; pick best by deltaR
    src=srcJets,  # RECO jets (any View<Jet> is ok)
    matched=cms.InputTag("ak8GenJetsWithNuSoftDrop"),  # GEN jets  (must be GenJetCollection)
    mcPdgId=cms.vint32(),  # n/a
    mcStatus=cms.vint32(),  # n/a
    checkCharge=cms.bool(False),  # n/a
    maxDeltaR=cms.double(jetR),  # Minimum deltaR for the match
    # maxDPtRel   = cms.double(3.0),                  # Minimum deltaPt/Pt for the match (not used in GenJetMatcher)
    resolveAmbiguities=cms.bool(True),  # Forbid two RECO objects to match to the same GEN object
    resolveByMatchQuality=cms.bool(False),  # False = just match input in order; True = pick lowest deltaR pair first
)
process.genJetTask = cms.Task(
    process.ak8GenJetsWithNu,
    process.ak8GenJetsWithNuMatch,
    process.ak8GenJetsWithNuSoftDrop,
    process.ak8GenJetsWithNuSoftDropMatch,
)

# DeepNtuplizer
process.load("DeepNTuples.Ntupler.DeepNtuplizer_cfi")
process.deepntuplizer.jets = srcJets
process.deepntuplizer.useReclusteredJets = useReclusteredJets
process.deepntuplizer.bDiscriminators = bTagDiscriminators + pfDeepBoostedJetTagsAll

process.deepntuplizer.genJetsMatch = 'ak8GenJetsWithNuMatch'
process.deepntuplizer.genJetsSoftDropMatch = 'ak8GenJetsWithNuSoftDropMatch'

process.deepntuplizer.isQCDSample = '/QCD_' in options.inputDataset
process.deepntuplizer.isPythia = 'pythia' in options.inputDataset.lower()
process.deepntuplizer.isHerwig = 'herwig' in options.inputDataset.lower()
process.deepntuplizer.isMadGraph = 'madgraph' in options.inputDataset.lower()  # note: MG can be interfaced w/ either pythia or herwig

process.deepntuplizer.isTrainSample = options.isTrainSample
if not options.inputDataset:
    # interactive running
    process.deepntuplizer.isTrainSample = False
#==============================================================================================================================#
process.p = cms.Path(process.deepntuplizer)
process.p.associate(patTask)
process.p.associate(process.genJetTask)
