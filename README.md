# RecoSimStudies

In this repository, you find all the necessary codes for the production of samples, from the generation to the dumping.

## Installation

    * cmsrel CMSSW_10_6_0
If you get an error, make sure that the remote machine on which you are working on is new enough to be compatible with the CMSSW_10_6_0 release. At the moment of writing, this release only works for machines with SL7 architecture at least, and one has typically to ask for a t3ui07 account to the PSI-T3 administrators.

    * cd CMSSW_10_6_0/src/
    * cmsenv
    * git cms-init
    * git cms-merge-topic bmarzocc:RecoSimStudies 
    * git clone git@github.com:pfclustering/RecoSimStudies.git
    * scram b -j 5

In case you want to interact with the Storage Element, don't forget to set up your proxy:
```    
voms-proxy-init --voms cms --valid 186:00
```


## Generation
### GEN_SIM production
```    
cd RecoSimStudies/Dumpers/test/
```

Note: the following instructions will allow you to generate single photon events in front of ECAL.

Choose the parameters you want in the "User's decision board" in launch_step1.sh. You may also want to change the output directories SERESULTDIR and TOPWORKDIR.

Finally, launch this bash script following the instructions written directly in the file. For instance, to run the file on the batch, do
```
sbatch -p wn -o logs/step1.out -e logs/step1.err -q long.q --ntasks=8 launch_step1.sh
```

Note that currently the wn partition doesn't allow to set the time limit above one day. Think of using the gpu partition instead.

### DIGI-RAW production

```                         
cd RecoSimStudies/Dumpers/test/
```

Choose the parameters you want in the "User's decision board" in launch_step2.sh. You may also want to change the output directories SERESULTDIR and TOPWORKDIR.

Finally, launch this bash script following the instructions written directly in the file. For instance, to run the file on the batch, do
```
sbatch -p wn -o logs/step2.out -e logs/step2.err -q long.q --ntasks=8 launch_step2.sh
```

### RECO production

```                         
cd RecoSimStudies/Dumpers/test/
```

Choose the parameters you want in the "User's decision board" in launch_step3.sh. You may also want to change the output directories SERESULTDIR and TOPWORKDIR.

Finally, launch this bash script following the instructions written directly in the file. For instance, to run the file on the batch, do
```
sbatch -p wn -o logs/step3.out -e logs/step3.err -q long.q --ntasks=8 launch_step3.sh
```

### Dumper
    
```                         
cd RecoSimStudies/Dumpers/python/
```

Run the dumper on a RECO sample. Example of commands are given in Cfg_RecoSimDumper_cfg.py. For instance, use

```
cmsRun Cfg_RecoSimDumper_cfg.py outputFile=../test/outputfiles/dumpedFiles/dumped_singlePhoton_5k_EB.root inputFiles=file:../test/outputfiles/singlePhoton_5k_EB/step3.root
```


