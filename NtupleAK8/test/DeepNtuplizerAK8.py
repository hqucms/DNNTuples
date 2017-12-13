import FWCore.ParameterSet.Config as cms

# ---------------------------------------------------------

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing('analysis')

options.outputFile = 'output_bb_official.root'
#options.inputFiles = 'file:/eos/user/a/anovak/022C3683-D4AB-E611-AC4D-3417EBE70078.root'  #include file: for local files, for catalogues /store...
#options.inputFiles = 'file:/afs/cern.ch/work/l/lmastrol/public/deep-c_tagging/HIG-RunIISummer16MiniAODv2-04133.root'
#options.inputFiles = 'file:/afs/cern.ch/work/l/lmastrol/public/deep-c_tagging/BulkGravHH4C-RunIISummer16MiniAODv2.root'
#options.inputFiles = 'file:/eos/user/a/anovak/DNNtuples/CMSSW_8_0_28/src/DeepNTuples/00D74F9D-9D14-E711-A4C0-1866DA879444.root'
#options.inputFiles = 'file:/afs/cern.ch/work/l/lmastrol/public/deep-c_tagging/GluGluHToCC_official_highStat_0219759C-1DD6-E711-BCCA-02163E0144F7.root'
options.inputFiles = 'file:/afs/cern.ch/work/l/lmastrol/public/deep-c_tagging/GluGluHToBB_official_highStat_E6E6BC68-F543-E611-AB6F-003048CD716E.root'
options.maxEvents = -1

options.register('inputScript', '', VarParsing.multiplicity.singleton, VarParsing.varType.string, "input Script")
options.register('skipEvents', 0, VarParsing.multiplicity.singleton, VarParsing.varType.int, "skip N events")
options.register('job', 0, VarParsing.multiplicity.singleton, VarParsing.varType.int, "job number")
options.register('nJobs', 1, VarParsing.multiplicity.singleton, VarParsing.varType.int, "total jobs")
options.register('fjKeepFlavors', [], VarParsing.multiplicity.list, VarParsing.varType.int, "Types of fatjet to keep in this sample")
# options.register('gluonReduction', 0.0, VarParsing.multiplicity.singleton, VarParsing.varType.float, "gluon reduction")
options.register('inputDataset',
                 '',
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 "Input dataset")

options.setupTags(tag='%d', ifCond='nJobs > 1', tagArg='job')
options.parseArguments()

# ---------------------------------------------------------

process = cms.Process("DNNFiller")

process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 10       
if not options.inputScript:  # this is probably for testing
	process.MessageLogger.cerr.FwkReport.reportEvery = 10

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

if options.inputScript:
    process.load(options.inputScript)

numberOfFiles = len(process.source.fileNames)
numberOfJobs = options.nJobs
jobNumber = options.job

print jobNumber
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
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')
# process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data', '')

# ---------------------------------------------------------

bTagInfos = [
	'pfImpactParameterTagInfos',
	'pfInclusiveSecondaryVertexFinderTagInfos',
	'pfDeepCSVTagInfos',
	'pfBoostedDoubleSVAK8TagInfos'
]

bTagDiscriminators = [
	'softPFMuonBJetTags',
	'softPFElectronBJetTags',
	'pfJetBProbabilityBJetTags',
	'pfJetProbabilityBJetTags',
	'pfCombinedInclusiveSecondaryVertexV2BJetTags',
	'pfDeepCSVJetTags:probudsg',
	'pfDeepCSVJetTags:probb',
	'pfDeepCSVJetTags:probc',
	'pfDeepCSVJetTags:probbb',
	'pfDeepCSVJetTags:probcc',
	'pfBoostedDoubleSecondaryVertexAK8BJetTags'
]

jetCorrectionsAK8 = ('AK8PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute'], 'None')

from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection
updateJetCollection(
        process,
        labelName = "DeepFlavour",
        jetSource = cms.InputTag('slimmedJetsAK8'),
        jetCorrections = jetCorrectionsAK8,
        pfCandidates = cms.InputTag('packedPFCandidates'),
        pvSource = cms.InputTag("offlineSlimmedPrimaryVertices"),
        svSource = cms.InputTag('slimmedSecondaryVertices'),
        muSource = cms.InputTag('slimmedMuons'),
        elSource = cms.InputTag('slimmedElectrons'),
        btagInfos = bTagInfos,
        btagDiscriminators = bTagDiscriminators,
        explicitJTA = False
)

if hasattr(process,'updatedPatJetsTransientCorrectedDeepFlavour'):
	process.updatedPatJetsTransientCorrectedDeepFlavour.addTagInfos = cms.bool(True) 
	process.updatedPatJetsTransientCorrectedDeepFlavour.addBTagInfo = cms.bool(True)
else:
	raise ValueError('I could not find updatedPatJetsTransientCorrectedDeepFlavour to embed the tagInfos, please check the cfg')
# ---------------------------------------------------------

'''
from RecoJets.JetProducers.ak4GenJets_cfi import ak8GenJets
process.ak8GenJetsWithNu = ak4GenJets.clone(src='packedGenParticles')
 
 ## Filter out neutrinos from packed GenParticles
process.packedGenParticlesForJetsNoNu = cms.EDFilter("CandPtrSelector", src = cms.InputTag("packedGenParticles"), cut = cms.string("abs(pdgId) != 12 && abs(pdgId) != 14 && abs(pdgId) != 16"))
 ## Define GenJets
process.ak8GenJetsRecluster = ak8GenJets.clone(src='packedGenParticlesForJetsNoNu')

 
process.patGenJetMatchWithNu = cms.EDProducer("GenJetMatcher",  # cut on deltaR; pick best by deltaR           
    src         = cms.InputTag("selectedUpdatedPatJetsDeepFlavour"),      # RECO jets (any View<Jet> is ok) 
    matched     = cms.InputTag("ak8GenJetsWithNu"),        # GEN jets  (must be GenJetCollection)              
    mcPdgId     = cms.vint32(),                      # n/a   
    mcStatus    = cms.vint32(),                      # n/a   
    checkCharge = cms.bool(False),                   # n/a   
    maxDeltaR   = cms.double(0.8),                   # Minimum deltaR for the match   
    #maxDPtRel   = cms.double(3.0),                  # Minimum deltaPt/Pt for the match (not used in GenJetMatcher)                     
    resolveAmbiguities    = cms.bool(True),          # Forbid two RECO objects to match to the same GEN object 
    resolveByMatchQuality = cms.bool(False),         # False = just match input in order; True = pick lowest deltaR pair first          
)

process.patGenJetMatchRecluster = cms.EDProducer("GenJetMatcher",  # cut on deltaR; pick best by deltaR           
    src         = cms.InputTag("selectedUpdatedPatJetsDeepFlavour"),      # RECO jets (any View<Jet> is ok) 
    matched     = cms.InputTag("ak8GenJetsRecluster"),        # GEN jets  (must be GenJetCollection)              
    mcPdgId     = cms.vint32(),                      # n/a   
    mcStatus    = cms.vint32(),                      # n/a   
    checkCharge = cms.bool(False),                   # n/a   
    maxDeltaR   = cms.double(0.8),                   # Minimum deltaR for the match   
    #maxDPtRel   = cms.double(3.0),                  # Minimum deltaPt/Pt for the match (not used in GenJetMatcher)                     
    resolveAmbiguities    = cms.bool(True),          # Forbid two RECO objects to match to the same GEN object 
    resolveByMatchQuality = cms.bool(False),         # False = just match input in order; True = pick lowest deltaR pair first          
)

process.genJetSequence = cms.Sequence(process.packedGenParticlesForJetsNoNu*process.ak4GenJetsWithNu*process.ak4GenJetsRecluster*process.patGenJetMatchWithNu*process.patGenJetMatchRecluster)

# ---------------------------------------------------------

# Very Loose IVF SV collection
from PhysicsTools.PatAlgos.tools.helpers import loadWithPrefix
loadWithPrefix(process, 'RecoVertex.AdaptiveVertexFinder.inclusiveVertexing_cff', "looseIVF")
process.looseIVFinclusiveCandidateVertexFinder.primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices")
process.looseIVFinclusiveCandidateVertexFinder.tracks = cms.InputTag("packedPFCandidates")
process.looseIVFinclusiveCandidateVertexFinder.vertexMinDLen2DSig = cms.double(0.)
process.looseIVFinclusiveCandidateVertexFinder.vertexMinDLenSig = cms.double(0.)
process.looseIVFinclusiveCandidateVertexFinder.fitterSigmacut = 20

process.looseIVFcandidateVertexArbitrator.primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices")
process.looseIVFcandidateVertexArbitrator.tracks = cms.InputTag("packedPFCandidates")
process.looseIVFcandidateVertexArbitrator.secondaryVertices = cms.InputTag("looseIVFcandidateVertexMerger")
process.looseIVFcandidateVertexArbitrator.fitterSigmacut = 20
'''

# ---------------------------------------------------------

# DeepNtuplizer
process.load("DeepNTuples.NtupleAK8.DeepNtuplizerAK8_cfi")
process.deepntuplizer.jets = cms.InputTag('selectedUpdatedPatJetsDeepFlavour')
process.deepntuplizer.bDiscriminators = bTagDiscriminators 
process.deepntuplizer.bDiscriminators.append('pfCombinedMVAV2BJetTags')
process.deepntuplizer.LooseSVs = cms.InputTag("looseIVFinclusiveCandidateSecondaryVertices")

process.deepntuplizer.fjKeepFlavors = cms.untracked.vuint32(options.fjKeepFlavors)
process.deepntuplizer.isQCDSample = '/QCD_' in options.inputDataset

# process.deepntuplizer.gluonReduction = cms.double(options.gluonReduction)

# process.p = cms.Path(process.QGTagger + process.genJetSequence * process.deepntuplizer)
process.p = cms.Path(process.deepntuplizer)
