# RecoSimStudies

This is a customised version of https://github.com/bmarzocc/RecoSimStudies
In this repository, you find all the necessary codes for the production of samples, from the generation to the dumping.

## Installation

First installation:
```
cmsrel CMSSW_10_6_1_patch1
```
If you get an error, make sure that the remote machine on which you are working on is new enough to be compatible with the CMSSW_10_6_1_patch1 release. At the moment of writing, this release only works for machines with SL7 architecture at least, and one has typically to ask for a t3ui07 account to the PSI-T3 administrators.

```
cd CMSSW_10_6_1_patch1/src/
cmsenv
git cms-init
git checkout -b base
```
either:
```
git cms-merge-topic bmarzocc:PR_CaloParticles
git cms-merge-topic bmarzocc:PR_EcalPFSeedingThresholds
git cms-merge-topic bmarzocc:PR_ParticleGuns
```
or (if you also want to run an overlap study):
```
git cms-merge-topic mgratti:mg-PR_ParticleGuns
```

```
git clone git@github.com:pfclustering/RecoSimStudies.git
scram b -j 8
```

In case you want to interact with the Storage Element, don't forget to set up your proxy:
```    
voms-proxy-init --voms cms --valid 186:00
```

## Workflow

### Development of ```cmssw```
Developments Development of ```cmssw``` by members of pfclustering team are done via fork of cmssw.
Currently there are three topics that can be changed (see above), under bmarzocc repo.
Changes to a given topic are pushed first to own fork, and then PR is done https://github.com/bmarzocc/cmssw/tree/<TOPIC_BRANCH>

*IMPORTANT NOTE* Badder is using CMSSW_10_6_4, but the topic branches are still at CMSSW_10_6_0!

The developments should be tested in the full (meaning three topics) configuration, but only the commits relevant to a given topic should be pushed
to the relevant topic, with the following workflow:

* within a clean area, create new local branch with following convention and do developments:
```
cd CMSSW_X_Y_Z/src
git checkout -b mg-<TOPIC>
git remote add my-cmssw git@github.com:mgratti/cmssw.git 
git remote add badder-cmssw https://github.com/bmarzocc/cmssw.git
git pull badder-cmssw <TOPIC>
```
* make the relevant changes (which should already have been tested)
```
git add bla.cpp
git commit -m "bla" 
```
* push to remote branch, with same name convention:
```
git push my-cmssw mg-PR_<TOPIC>
```
* pull request to relevant topic branch under bmarzocc repo

* [OBSOLETE] once PR has been accepted, go back to `base` branch, and DELETE both local and remote branches
```
git checkout base
git branch -d mg-PR_<topic>
git push my-cmssw --delete mg-PR_<topic>
```
* [OBSOLETE] re-do all the relevant merge-topic 

More information and tricks on how to work with cmssw and github here: http://cms-sw.github.io/faq.html


### Development of ```RecoSimStudies```
Development of ```RecoSimStudies``` by members of pfclustering team happens within the pfclustering fork; 
each member has his/her own branch where to develop the new features. When development is over, he/she opens a pull request,
(ideally another member checks) and merges with master branch.

After master is in sync, developments of bmarzocc/RecoSimStudies are fetched via a pull request (from web page) with a brief comment about the changes.

## Generation
For all steps of generation until reco files 
```
cd Dumpers/test/ECALproductionHelper
```
See available options:
```
python prodHelper.py --help
```
Example commands in ```Dumpers/test/ECALproductionHelper/README.md```

To resubmit jobs that failed for step 3:
```
python resubmitHelper.py --help
```
After production is over, you can run post-production to create list of files and cmsRun command for next step
```
python postProdHelper.py --help
```

You can find examples of commands run for past productions in [this spreadsheet](https://docs.google.com/spreadsheets/d/1rJ4FppHxdD6TQpCHnnnX0QD0Y0d4aTEN6ibgTVn7Vyg/edit#gid=0)

### Dumper
```                         
cd RecoSimStudies/Dumpers/python/
```

Run the dumper on a RECO sample. Example of commands are given in Cfg_RecoSimDumper_cfg.py. For instance, use

```
cmsRun Cfg_RecoSimDumper_cfg.py outputFile=../test/outputfiles/dumpedFiles/dumped_singlePhoton_5k_EB.root inputFiles=file:../test/outputfiles/singlePhoton_5k_EB/step3.root
```
