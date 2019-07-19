#!/bin/bash

###############################################################
#How to launch this script:
#if you want to run it locally: 
#
# source launch_step3.sh
#
#if you want it to run it with slurm, two possibilities:
#
#1. either run it with the wn partition (will use by default the processor t3wn38, 2.6GHz, 16cores) 
#
# sbatch -p wn -o logs/step3.out -e logs/step3.err -q long.q launch_step3.sh
#
#2. or use the gpu ressources
#
# sbatch --account=gpu_gres --partition=gpu --gres=gpu:2 --job-name=step3 -o logs/step3.out -e logs/step3.err launch_step3.sh 
#
# Add nodes: --nodes=4 (max for wn) --nodes=2 (max for gpu)
###############################################################





# Job configuration
JOBOPFILENAME="step3_RAW2DIGI_L1Reco_RECO_RECOSIM_EI_PAT_VALIDATION_DQM.py"
FILENAME="step3.root"
INFILENAME="step2.root"
#SERESULTDIR="/pnfs/psi.ch/cms/trivcat/store/user/anlyon/EcalProd/step1/test/"
SERESULTDIR="/t3home/anlyon/CMSSW_10_6_0/src/RecoSimStudies/Dumpers/test/outputfiles/singlePhoton_5k"

STARTDIR=`pwd`
TOPWORKDIR="/scratch/anlyon/"
JOBDIR="gen_"$SERESULTDIR
WORKDIR=$TOPWORKDIR/$JOBDIR
SEPREFIX="root://t3dcachedb.psi.ch:1094/"




# Job instructions
source $VO_CMS_SW_DIR/cmsset_default.sh
shopt -s expand_aliases

echo ""
echo "Going to set up cms environment"
cd $STARTDIR
cmsenv

echo ""
echo "Going to create work dir"
mkdir -p $WORKDIR


echo ""
echo "Going to create the output dir"
mkdir -p $SERESULTDIR


echo ""
echo "Going to copy cms driver"
cp $JOBOPFILENAME $WORKDIR/$JOBOPFILENAME


echo ""
echo "Going to copy input file"
#if interaction with the storage element:
#xrdcp  $SEPREFIX/$SERESULTDIR/$INFILENAME $WORKDIR/$INFILENAME
#otherwise:
cp $SERESULTDIR/$INFILENAME $WORKDIR/$INFILENAME

cd $WORKDIR

echo ""
echo "Going to run"

DATE_START=`date +%s`
cmsRun $JOBOPFILENAME 
DATE_END=`date +%s`

echo ""
echo "Finished running"

echo ""
echo "Content of current directory"
ls -al

echo ""
echo "Going to copy the output to the output directory"

#if interaction with the storage element:
#xrdcp $FILENAME $SEPREFIX/$SERESULTDIR/$FILENAME
#otherwise:
cp $FILENAME $SERESULTDIR/$FILENAME

#echo ""
#echo "Cleaning up $WORKDIR"
#rm -rf $WORKDIR

RUNTIME=$((DATE_END-DATE_START))
echo "Wallclock running time: $RUNTIME s"

cd $STARTDIR


