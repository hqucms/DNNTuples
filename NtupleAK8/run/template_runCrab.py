import os
from CRABClient.UserUtilities import config
config = config()

config.General.requestName = '_requestName_'
config.General.workArea = '_workArea_'
config.General.transferOutputs = True
config.General.transferLogs = False

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '_psetName_'
# config.JobType.sendExternalFolder = True

config.Data.inputDataset = '_inputDataset_'
config.JobType.pyCfgParams = _pyCfgParams_
# config.JobType.inputFiles = []

config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = _unitsPerJob_
config.Data.outLFNDirBase = '_outLFNDirBase_'
config.Data.publication = False
#config.Data.allowNonValidInputDataset = True

# uncomment this part to use CRAB to submit to FNAL_LPC
# config.Data.ignoreLocality = True
# config.Site.whitelist = ['T3_US_FNALLPC']
# config.Site.ignoreGlobalBlacklist = True

config.Site.storageSite = '_storageSite_'
