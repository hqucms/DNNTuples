import FWCore.ParameterSet.Config as cms

# ---------------------------------------------------------

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing('analysis')

options.outputFile = 'output.root'
options.inputFiles = '/store/mc/RunIISummer16MiniAODv2/GluGluHToCC_M125_13TeV_powheg_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/00DE64C2-3EE6-E711-BCB8-0242AC130002.root'
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
# process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')
# process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data', '')
process.GlobalTag = GlobalTag(process.GlobalTag, '80X_mcRun2_asymptotic_2016_TrancheIV_v8', '')
print 'Using global tag', process.GlobalTag.globaltag

# ---------------------------------------------------------
usePuppi = True
jetR = 1.5

bTagInfos = [
    'pfBoostedDoubleSVAK8TagInfos'
]

bTagDiscriminators = [
    'pfCombinedInclusiveSecondaryVertexV2BJetTags',
    'pfBoostedDoubleSecondaryVertexAK8BJetTags'
]


JETCorrLevels = ['L2Relative', 'L3Absolute']

from JMEAnalysis.JetToolbox.jetToolbox_cff import jetToolbox
jetToolbox(process, 'ak15', 'jetSequence', 'out', PUMethod='Puppi', JETCorrPayload='AK8PFPuppi', JETCorrLevels=JETCorrLevels, miniAOD=True, runOnMC=True,
           Cut='pt > 120.0 && abs(rapidity()) < 2.4', addNsub=True, maxTau=3,
           addSoftDrop=True, addSoftDropSubjets=True, subJETCorrPayload='AK4PFPuppi', subJETCorrLevels=JETCorrLevels,
           bTagDiscriminators=bTagDiscriminators, bTagInfos=bTagInfos)
srcJets = cms.InputTag('selectedPatJetsAK15PFPuppi')
srcSubjets = cms.InputTag('selectedPatJetsAK15PFPuppiSoftDropPacked')

# ---------------------------------------------------------
# genJets
from RecoJets.JetProducers.ak4GenJets_cfi import ak4GenJets
process.ak15GenJetsWithNu = ak4GenJets.clone(
    src='packedGenParticles',
    rParam=cms.double(jetR),
    jetPtMin=100.0
)
process.ak15GenJetWithNuMatch = cms.EDProducer("GenJetMatcher",  # cut on deltaR; pick best by deltaR
    src=srcJets,  # RECO jets (any View<Jet> is ok)
    matched=cms.InputTag("ak15GenJetsWithNu"),  # GEN jets  (must be GenJetCollection)
    mcPdgId=cms.vint32(),  # n/a
    mcStatus=cms.vint32(),  # n/a
    checkCharge=cms.bool(False),  # n/a
    maxDeltaR=cms.double(jetR),  # Minimum deltaR for the match
    # maxDPtRel   = cms.double(3.0),                  # Minimum deltaPt/Pt for the match (not used in GenJetMatcher)
    resolveAmbiguities=cms.bool(True),  # Forbid two RECO objects to match to the same GEN object
    resolveByMatchQuality=cms.bool(False),  # False = just match input in order; True = pick lowest deltaR pair first
)
# softdrop
process.ak15GenJetsWithNuSoftDrop = process.ak15GenJetsWithNu.clone(
    useSoftDrop=cms.bool(True),
    zcut=cms.double(0.1),
    beta=cms.double(0.0),
    R0=cms.double(jetR),
    useExplicitGhosts=cms.bool(True),
)
process.ak15GenJetWithNuSoftDropMatch = cms.EDProducer("GenJetMatcher",  # cut on deltaR; pick best by deltaR
    src=srcJets,  # RECO jets (any View<Jet> is ok)
    matched=cms.InputTag("ak15GenJetsWithNuSoftDrop"),  # GEN jets  (must be GenJetCollection)
    mcPdgId=cms.vint32(),  # n/a
    mcStatus=cms.vint32(),  # n/a
    checkCharge=cms.bool(False),  # n/a
    maxDeltaR=cms.double(jetR),  # Minimum deltaR for the match
    # maxDPtRel   = cms.double(3.0),                  # Minimum deltaPt/Pt for the match (not used in GenJetMatcher)
    resolveAmbiguities=cms.bool(True),  # Forbid two RECO objects to match to the same GEN object
    resolveByMatchQuality=cms.bool(False),  # False = just match input in order; True = pick lowest deltaR pair first
)
process.genJetSequence = cms.Sequence(process.ak15GenJetsWithNu * process.ak15GenJetWithNuMatch * process.ak15GenJetsWithNuSoftDrop * process.ak15GenJetWithNuSoftDropMatch)
# ---------------------------------------------------------

# DeepNtuplizer
process.load("DeepNTuples.Ntupler.DeepNtuplizer_cfi")
process.deepntuplizer.jets = srcJets
process.deepntuplizer.subjets = srcSubjets
process.deepntuplizer.genJetsMatch = cms.InputTag('ak15GenJetWithNuMatch')
process.deepntuplizer.genJetsSoftDropMatch = cms.InputTag('ak15GenJetWithNuSoftDropMatch')
process.deepntuplizer.usePuppi = cms.bool(usePuppi)
process.deepntuplizer.bDiscriminators = bTagDiscriminators
process.deepntuplizer.jetR = jetR
process.deepntuplizer.jetPtMin = 150

process.deepntuplizer.isQCDSample = '/QCD_' in options.inputDataset
process.deepntuplizer.isTrainSample = options.isTrainSample
if not options.inputDataset:
    # interactive running
    process.deepntuplizer.isTrainSample = False
#==============================================================================================================================#
fillElectronVars = False
fillMuonVars = False

process.deepntuplizer.fillElectronVars = fillElectronVars
process.deepntuplizer.fillMuonVars = fillMuonVars

if fillElectronVars:
    # Electron ID, following prescription in
    # https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedElectronIdentificationRun2
    # set up everything that is needed to compute electron IDs and
    # add the ValueMaps with ID decisions into the event data stream

    # Load tools and function definitions
    from PhysicsTools.SelectorUtils.tools.vid_id_tools import *

    switchOnVIDElectronIdProducer(process, DataFormat.MiniAOD)

    # Define which IDs we want to produce
    my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Summer16_80X_V1_cff']

    # Add them to the VID producer
    for idmod in my_id_modules:
        setupAllVIDIdsInModule(process, idmod, setupVIDElectronSelection)

    # Set ID tags
    process.deepntuplizer.eleVetoIds = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-veto")
    process.deepntuplizer.eleLooseIds = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-loose")
    process.deepntuplizer.eleMediumIds = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-medium")
    process.deepntuplizer.eleTightIds = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-tight")
# ---------------------------------------------------------
process.p = cms.Path(process.genJetSequence * process.deepntuplizer)
