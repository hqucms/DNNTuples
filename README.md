# DNNTuplesAK8

## Setup
```
cmsrel CMSSW_8_0_30
cd CMSSW_8_0_30/src/
cmsenv
# setup JetToolBox
# https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetToolbox
git clone git@github.com:cms-jet/JetToolbox.git JMEAnalysis/JetToolbox -b jetToolbox_80X_V3
# clone this repo into "DeepNTuples" directory
git clone ssh://git@gitlab.cern.ch:7999/hqu/DNNTuplesAK8.git DeepNTuples
scram b -j8
```

## Submit jobs via CRAB

```bash
# set up CRAB env; run it after cmsenv
source /cvmfs/cms.cern.ch/crab3/crab.sh
voms-proxy-init -rfc -voms cms --valid 168:00

# create the CRAB config files
cd DeepNTuples/NtupleAK8/run
./makeCrabJobs.py -i [samples/ttbar.conf] -o [/eos/cms/store/user/$USER/DeepNtuples/output_dir] --site [T2_CH_CERN|T3_US_FNALLPC|...]
# submit jobs
./submit_[ttbar].conf
```

To check all the options of the submission script, run
```
./makeCrabJobs.py -h
```
 