# DNNTuplesAK8

## Setup
```
cmsrel CMSSW_10_2_7
cd CMSSW_10_2_7/src/
cmsenv

# clone this repo into "DeepNTuples" directory
git clone https://github.com/hqucms/DNNTuplesAK8.git DeepNTuples -b 94X
scram b -j24
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
 