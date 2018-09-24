# DNNTuplesAK8

## Setup
```
cmsrel CMSSW_10_2_5
cd CMSSW_10_2_5/src/
cmsenv

# MXNet
scram setup /cvmfs/cms.cern.ch/slc6_amd64_gcc700/cms/cmssw/CMSSW_10_3_0_pre4/config/toolbox/slc6_amd64_gcc700/tools/selected/mxnet-predict.xml

# get DeepAK8 PR
git cms-merge-topic -u hqucms:deep-boosted-jets-rebase-102X

# clone this repo into "DeepNTuples" directory
git clone ssh://git@gitlab.cern.ch:7999/hqu/DNNTuplesAK8.git DeepNTuples -b 94X
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
 