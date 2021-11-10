import FWCore.ParameterSet.Config as cms

# ---------------------------------------------------------
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing('analysis')

options.outputFile = 'output.root'
options.inputFiles = '/store/mc/RunIISummer19UL17MiniAOD/BulkGravitonToHHTo4Q_MX-600to6000_MH-15to250_part2_TuneCP5_13TeV-madgraph_pythia8/MINIAODSIM/multigridpack_106X_mc2017_realistic_v6-v1/50000/FB46C2C2-73A4-A64C-A3D7-FC47C6A48871.root'
# options.inputFiles = '/store/mc/RunIISummer19UL17MiniAOD/QCD_Pt-15to7000_TuneCP5_Flat_13TeV_pythia8/MINIAODSIM/106X_mc2017_realistic_v6_ext2-v2/60000/E0502E95-AA7E-B441-9277-498113BA458C.root'
# options.inputFiles = '/store/mc/RunIISummer19UL17MiniAOD/TTToSemiLeptonic_mtop171p5_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/106X_mc2017_realistic_v6-v2/60000/F06BD23E-48AB-1F4D-A0FB-80DD370F4868.root'
options.maxEvents = -1

options.register('skipEvents', 0, VarParsing.multiplicity.singleton, VarParsing.varType.int, "skip N events")
options.register('inputDataset',
                 '',
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 "Input dataset")
options.register('isTrainSample', True, VarParsing.multiplicity.singleton, VarParsing.varType.bool, "if the sample is used for training")

options.parseArguments()

globalTagMap = {
    'auto': 'auto:phase1_2021_realistic',
}

era = 'auto' if options.inputDataset else 'auto'
for k in globalTagMap:
    if k in options.inputDataset:
        era = k
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
# ---------------------------------------------------------
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.EventContent.EventContent_cff")
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, globalTagMap[era], '')
print('Using global tag', process.GlobalTag.globaltag)
process.TransientTrackBuilderESProducer = cms.ESProducer("TransientTrackBuilderESProducer",
    ComponentName=cms.string('TransientTrackBuilder')
)
# ---------------------------------------------------------
from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection

isPuppiJets = False
jetR = 0.4

from RecoBTag.ONNXRuntime.pfParticleNetAK4_cff import _pfParticleNetAK4JetTagsAll


bTagDiscriminators = [
    'pfDeepFlavourJetTags:probb',
    'pfDeepFlavourJetTags:probbb',
    'pfDeepFlavourJetTags:problepb',
    'pfDeepFlavourJetTags:probc',
    'pfDeepFlavourJetTags:probuds',
    'pfDeepFlavourJetTags:probg',
] + _pfParticleNetAK4JetTagsAll

updateJetCollection(
    process,
    jetSource=cms.InputTag('slimmedJets'),
    jetCorrections=('AK4PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None'),
)
srcJets = cms.InputTag('selectedUpdatedPatJets')
# ---------------------------------------------------------
from RecoJets.JetProducers.QGTagger_cfi import QGTagger
process.qgtagger = QGTagger.clone(srcJets=srcJets, srcVertexCollection="offlineSlimmedPrimaryVertices")
process.updatedJetsWithUserData = cms.EDProducer("PATJetUserDataEmbedder",
    src=srcJets,
    userFloats=cms.PSet(qgl=cms.InputTag('qgtagger:qgLikelihood')),
    )
srcJets = cms.InputTag('updatedJetsWithUserData')
process.qgTask = cms.Task(
    process.qgtagger,
    process.updatedJetsWithUserData,
)
# ---------------------------------------------------------
from PhysicsTools.PatAlgos.tools.helpers import getPatAlgosToolsTask
patTask = getPatAlgosToolsTask(process)

from RecoJets.JetProducers.ak4GenJets_cfi import ak4GenJets
process.ak4GenJetsWithNu = ak4GenJets.clone(
    src='packedGenParticles'
    )
process.ak4GenJetsWithNuMatch = cms.EDProducer("GenJetMatcher",  # cut on deltaR; pick best by deltaR
    src=srcJets,  # RECO jets (any View<Jet> is ok)
    matched=cms.InputTag("ak4GenJetsWithNu"),  # GEN jets  (must be GenJetCollection)
    mcPdgId=cms.vint32(),  # n/a
    mcStatus=cms.vint32(),  # n/a
    checkCharge=cms.bool(False),  # n/a
    maxDeltaR=cms.double(jetR),  # Minimum deltaR for the match
    # maxDPtRel   = cms.double(3.0),                  # Minimum deltaPt/Pt for the match (not used in GenJetMatcher)
    resolveAmbiguities=cms.bool(True),  # Forbid two RECO objects to match to the same GEN object
    resolveByMatchQuality=cms.bool(False),  # False = just match input in order; True = pick lowest deltaR pair first
)
process.genJetTask = cms.Task(
    process.ak4GenJetsWithNu,
    process.ak4GenJetsWithNuMatch,
)

# DeepNtuplizer
process.load("DeepNTuples.Ntupler.DeepNtuplizer_cfi")
process.deepntuplizer.jets = srcJets
process.deepntuplizer.isPuppiJets = isPuppiJets
process.deepntuplizer.bDiscriminators = bTagDiscriminators

process.deepntuplizer.genJetsMatch = 'ak4GenJetsWithNuMatch'

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
process.p.associate(process.qgTask)
