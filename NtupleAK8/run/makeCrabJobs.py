#!/usr/bin/env python
import os
import sys
import re
import argparse

parser = argparse.ArgumentParser(description='Prepare crab cfg for ntuple jobs')
parser.add_argument("-i", "--input", default="datasets.txt", help="Input configuration file with list of datasets. [Default: %(default)s]")
parser.add_argument("-s", "--submit", default="submit", help="Name of shell script to run for job submission. [Default: %(default)s]")
parser.add_argument("-n", "--numperjob", type=int, default=5, help="Number of files or events per job. Splittype determines whether splitting is by number of files (default) or by number of events. [Default: %(default)d]")
parser.add_argument("-c", "--config", default="../test/DeepNtuplizerAK8.py", help="Configuration file to be run using cmsRun to run. [Default: %(default)s]")
parser.add_argument("-o", "--outdir", default="/store/user/${USER}/ntuples", help="Output directory for ntuples. [Default: %(default)s")
parser.add_argument("-j", "--jobdir", default="jobs", help="Directory for job files. [Default: %(default)s]")
parser.add_argument("--site", default="T2_CH_CERN", help="Crab storage site. [Default: %(default)s]")
parser.add_argument("--crabcmd", default=None, help="Crab command. [Default: None]")
parser.add_argument("--crab-workarea", default='crab_projects', help="Crab work area. [Default: %(default)s]")
args = parser.parse_args()

jobname = os.path.splitext(os.path.basename(args.input))[0]
args.outdir = os.path.join(args.outdir, jobname)
args.submit = args.submit + '_' + jobname

samples = []
datasets = []
arguments = []

def parseConfig():
    with open(args.input, "r") as f :
        for line in f :
            line = line.strip()
            if not line: continue
            if '=' in line:
                arguments.append(line.replace(' ', ''))
                continue
            if line.startswith('#'):
                continue

            datasets.append(line)
            sname = line.split('/')[1]
            r = re.search(r'(_ext[0-9]+)', line.split('/')[2])
            if r:
                sname += r.groups()[0]
            samples.append(sname)

os.system("mkdir -p %s" % args.jobdir)
parseConfig()

if args.crabcmd:
    submitfile = '%s_cmd.sh' % args.submit
    print "Creating submission file: ", submitfile
    with open(submitfile, "w") as script:
        script.write("#!/bin/bash\n\n")
        for samp in samples:
            script.write("crab %s -d crab_projects/crab_%s\necho '======'\necho\n\n" % (args.crabcmd, samp))
    os.system("chmod +x %s" % submitfile)
    print "Done!"
    exit()
else:
    print "Creating crab config files: "
    crab_configs = []
    pyCfgParams = ['%s=%s']
    for isamp in range(len(samples)):
        samp = samples[isamp]
        dataset = datasets[isamp]
        if not dataset: continue
        crab_template_file = 'template_runCrab.py'
        with open(crab_template_file, 'r') as crfile:
            crabconfig = crfile.read()
        replace_dict = {'_requestName_':samp,
                        '_workArea_':args.crab_workarea,
                        '_psetName_':args.config,
                        '_inputDataset_':dataset,
                        '_pyCfgParams_':repr(arguments),
                        '_unitsPerJob_':str(args.numperjob),
                        '_outLFNDirBase_':re.sub(r'/eos/[a-z]+/store', '/store', args.outdir),
                        '_storageSite_':args.site}
        for key in replace_dict:
            crabconfig = crabconfig.replace(key, replace_dict[key])
        cfgfilename = os.path.join(args.jobdir , 'crab_submit_%s.py' % samp)
        crab_configs.append(cfgfilename)
        with open(cfgfilename, 'w') as cfgfile:
            cfgfile.write(crabconfig)
        print cfgfilename

    submitfile = '%s.sh' % args.submit
    print "Creating submission file: ", submitfile
    with open(submitfile, "w") as script:
        script.write("#!/bin/bash\n\n")
        for cfg in crab_configs:
            script.write("crab submit -c %s\necho '======'\necho\n\n" % cfg)
        script.write("\necho 'Done!'\n")

    os.system("chmod +x %s" % submitfile)
    print "Done!"
    exit()

