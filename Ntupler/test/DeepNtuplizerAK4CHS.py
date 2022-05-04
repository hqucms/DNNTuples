import FWCore.ParameterSet.Config as cms

# ---------------------------------------------------------
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing('analysis')

options.outputFile = 'output.root'
options.inputFiles = '/store/mc/RunIISummer20UL18MiniAOD/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v11_L1v1-v2/20000/0EE3441F-AF29-7D41-935B-3A8209BD4DD2.root'
# options.inputFiles = '/store/cmst3/group/vhcc/hc/samples/HPlusCharm_H4L_JHUGEN_3FS_M125_13TeV_amcatnlo_pythia8_MINIAOD_v3/UL2018-MINIAOD/211028_065528/0001/step6_1000.root'

#options.inputFiles = '/store/mc/RunIISummer19UL17MiniAOD/BulkGravitonToHHTo4Q_MX-600to6000_MH-15to250_part2_TuneCP5_13TeV-madgraph_pythia8/MINIAODSIM/multigridpack_106X_mc2017_realistic_v6-v1/50000/FB46C2C2-73A4-A64C-A3D7-FC47C6A48871.root'
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
    'Summer20UL17': '106X_mc2017_realistic_v9',
    'Summer20UL18': '106X_upgrade2018_realistic_v11_L1v1',
    'Summer20UL16': '',
    '106X_mc2017_realistic_v7': '106X_mc2017_realistic_v7',
}

era = None if options.inputDataset else 'Summer20UL17'

#globalTagMap = {
#    'Summer19UL17': '106X_mc2017_realistic_v6',
#    'Summer19UL18': '106X_upgrade2018_realistic_v11_L1v1',
#    'Summer19UL16': '',
#    '106X_mc2017_realistic_v7': '106X_mc2017_realistic_v7',
#}

era = None if options.inputDataset else 'Summer20UL18'  #'Summer19UL18'
for k in globalTagMap:
    if k in options.inputDataset:
        era = k
# ---------------------------------------------------------
process = cms.Process("DNNFiller")

process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.options = cms.untracked.PSet(
    allowUnscheduled = cms.untracked.bool(True),
    wantSummary=cms.untracked.bool(True)
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
# read JEC from sqlite
if era == 'Summer19UL17':
    import os
    jecTag = 'Summer19UL17_V5_MC'
    jecFile = '%s.db' % jecTag
    if not os.path.exists(jecFile):
        os.symlink('../data/'+jecFile, jecFile)
    from CondCore.CondDB.CondDB_cfi import CondDB
    CondDBJECFile = CondDB.clone(connect = cms.string( 'sqlite:%s'%jecFile ) )
    process.jec = cms.ESSource('PoolDBESSource',
        CondDBJECFile,
        toGet = cms.VPSet(
            cms.PSet(
                record = cms.string('JetCorrectionsRecord'),
                tag    = cms.string('JetCorrectorParametersCollection_%s_AK4PFchs' % jecTag),
                label  = cms.untracked.string('AK4PFchs')
            ),
            cms.PSet(
                record = cms.string('JetCorrectionsRecord'),
                tag    = cms.string('JetCorrectorParametersCollection_%s_AK4PFPuppi' % jecTag),
                label  = cms.untracked.string('AK4PFPuppi')
            ),
            # ...and so on for all jet types you need
        )
    )
    print(jecTag, process.jec.toGet)
    # Add an ESPrefer to override JEC that might be available from the global tag
    process.es_prefer_jec = cms.ESPrefer('PoolDBESSource', 'jec')
# ---------------------------------------------------------
from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection

isPuppiJets = False
jetR = 0.4  # change this (reject SVs if deltaR(sv, jet) < jetR)
pfcandR = 0.4  # (for each sV, store all pfcands w/ deltaR(sv, pf) < pfcandR; also used for matching SV to b and c hadrons)

bTagDiscriminators = [
    'pfDeepFlavourJetTags:probb',
    'pfDeepFlavourJetTags:probbb',
    'pfDeepFlavourJetTags:problepb',
    'pfDeepFlavourJetTags:probc',
    'pfDeepFlavourJetTags:probuds',
    'pfDeepFlavourJetTags:probg',
]

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


# NOTE:  ALSO NEW, from:  https://github.com/cms-sw/cmssw/blob/master/PhysicsTools/NanoAOD/python/jetMC_cff.py#L32
process.patJetPartonsNano = cms.EDProducer('HadronAndPartonSelector',
    src = cms.InputTag("generator"),
    particles = cms.InputTag("prunedGenParticles"),
    partonMode = cms.string("Auto"),
    fullChainPhysPartons = cms.bool(True)
)

process.genJetTask = cms.Task(
    process.ak4GenJetsWithNu,
    process.ak4GenJetsWithNuMatch,
    process.patJetPartonsNano,
)


# DeepNtuplizer
process.load("DeepNTuples.Ntupler.DeepNtuplizer_cfi")
process.deepntuplizer.jets = srcJets
process.deepntuplizer.isPuppiJets = isPuppiJets
process.deepntuplizer.bDiscriminators = bTagDiscriminators
# NOTE: NEW
process.deepntuplizer.jetR = jetR
process.deepntuplizer.pfcandR = pfcandR


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
