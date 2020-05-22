# https://docs.google.com/spreadsheets/d/1rJ4FppHxdD6TQpCHnnnX0QD0Y0d4aTEN6ibgTVn7Vyg/edit#gid=0

import sys
import os
import subprocess


def getOptions():
  from argparse import ArgumentParser

  parser = ArgumentParser(description='Production helper for Clustering studies', add_help=True)

  parser.add_argument('-v','--ver', type=str, dest='ver', help='version of production, e.g. V00_v00', default='V00_v00')
  #parser.add_argument('-r','--rel', type=str, dest='rel', help='cmssw release', default='10_6_1_patch1')
  parser.add_argument('-y', '--year', type=int, dest='year', help='year, defined conditions of CMS', default=2021, choices=[2021,2023])

  parser.add_argument('-n','--nevts', type=int, dest='nevts', help='total number of events to be generated', default=10)
  parser.add_argument('-c','--ch', type=str, dest='ch', help='channel, e.g. photon', default='photon', choices=['photon', 'gjetEM', 'QCD'])
  parser.add_argument('--etmax', type=str, dest='etmax', help='max Et (GeV)', default='100')
  parser.add_argument('--etmin', type=str, dest='etmin', help='min Et (GeV)', default='1')
  parser.add_argument('--doflatenergy', dest='doflatenergy', help='generate flat in energy, otherwise in pt', action='store_true', default=False)
  parser.add_argument('--npart', type=int, dest='npart', help='number of particles to generate per event for closeEcal configuration, specify only if you want to override the default', default=None)
  parser.add_argument('-g','--geo',type=str, dest='geo', help='detector configuration: wTk, noTk, closeEcal', default='closeEcal', choices=['wTk', 'noTk', 'closeEcal'])
  parser.add_argument('-d','--det', type=str, dest='det', help='sub-detector: EB, EEclose, EEfar or all', default='EB', choices=['EB', 'EEclose', 'EEfar', 'all'])

  parser.add_argument('--pu', type=str, dest='pu', help='PU configuration', default='noPU', choices=['noPU', 'wPU'])

  parser.add_argument('--thrslumi', type=int, dest='thrslumi', help='based on which integrated luminosity were the thrs computed', default=450, choices=[2021, 2023, 150, 180, 235, 315, 400, 450, 550])
  parser.add_argument('--pfrhmult', type=float, dest='pfrhmult', help='how many sigma of the noise to use for PFRH thresholds', default=1.)
  parser.add_argument('--pfrhmultbelow2p5', type=float, dest='pfrhmultbelow2p5', help='sigma of the noise for PFRH thresholds for |eta|<2.5', default=0.)
  parser.add_argument('--pfrhmultabove2p5', type=float, dest='pfrhmultabove2p5', help='sigma of the noise for PFRH thresholds for |eta|>2.5', default=0.)
  parser.add_argument('--seedmult', type=float, dest='seedmult', help='how many sigma of the noise to use for seeding thresholds', default=3.)
  parser.add_argument('--seedmultbelow2p5', type=float, dest='seedmultbelow2p5', help='sigma of the noise for seeding thresholds for |eta|<2.5', default=0.)
  parser.add_argument('--seedmultabove2p5', type=float, dest='seedmultabove2p5', help='sigma of the noise for seeding thresholds for |eta|>2.5', default=0.)
  parser.add_argument('--dorefpfrh', dest='dorefpfrh', help='use reference values for pfrh and gathering thresholds', action='store_true', default=False)
  parser.add_argument('--dorefseed', dest='dorefseed', help='use reference values for seeding thresholds', action='store_true', default=False)
  parser.add_argument('--doringavgEB', dest='doringavgEB', help='apply ring-averaged thresholds in EB', action='store_true', default=False)
  parser.add_argument('--doringavgEE', dest='doringavgEE', help='apply ring-averaged thresholds in EE', action='store_true', default=False)
  parser.add_argument('--dosafetymargin', dest='dosafetymargin', help='raise values on pfrh and seeding thrs of 10 percent', action='store_true', default=False)
  parser.add_argument('--showersigmamult', type=float, dest='showersigmamult', help='how large do you want the shower sigma', default=1.)
  parser.add_argument('--maxsigmadist', type=float, dest='maxsigmadist', help='max RH-cl distance in PFClustering, in units of sigma, THIS OPTION CURRENTLY DOES NOTHING', default=10.)
  parser.add_argument('--doreco', dest='doreco', help='do step 3 (reconstruction) and later stages (if any) starting from an existing production', action='store_true', default=False)
  parser.add_argument('--dorecofromeos', dest='dorecofromeos', help='do reconstruction from existing production on eos', action='store_true', default=False)
  parser.add_argument('--custominput', type=str, dest='custominput', help='full path of the input file that you want to use for the reconstruction (also SE is supported)', default=None)
  parser.add_argument('--custominputdir', type=str, dest='custominputdir', help='full path of the input directory that you want to use for the reconstruction (also SE is supported), only for multijob production', default=None)
  parser.add_argument('--doold', dest='doold', help='use old_ version of the scripts', action='store_true', default=False)

  #parser.add_argument('--dostep3only', dest='dostep3only', help='do only step 3', action='store_true', default=False)
  parser.add_argument('--doshort', dest='doshort', help='set 2 hours as wall clock time instead of 1 day', action='store_true', default=False)
  parser.add_argument('--domedium', dest='domedium', help='set 2 days as wall clock time instead of 1 day', action='store_true', default=False)
  parser.add_argument('--dolong', dest='dolong', help='set 3 days as wall clock time instead of 1 day', action='store_true', default=False)
  parser.add_argument('--docustomtime', dest='docustomtime', help='set custom time', action='store_true', default=False)
  parser.add_argument('--time', type=str, dest='time', help='requires docustomtime, allowed time for each job of each step in hours, including dumper, example 01,02,01,01,04 for 1h for step1, 2h for step2, 1h for step3, 1h for step4, 4h for dumper', default='02,02,05,01,05')
  parser.add_argument('--domultithread', dest='domultithread', help='run multithreaded', action='store_true', default=False)
  parser.add_argument('--domultijob', dest='domultijob', help='run several separate jobs', action='store_true', default=False)
  parser.add_argument('--njobs', type=int, dest='njobs', help='number of parallel jobs to submit', default=10)
  parser.add_argument('--splitfactord', type=int, dest='splitfactord', help='number of jobs for the dumper', default=1)
  parser.add_argument('--dosavehome', dest='dosavehome', help='save in home, otherwise save to SE', action='store_true', default=False)
  parser.add_argument('--doskipdumper', dest='doskipdumper', help='do not run the dumper at the end', action='store_true', default=False)
  parser.add_argument('--dodumperonly', dest='dodumperonly', help='only run the dumper', action='store_true', default=False)
  parser.add_argument('--pli', type=str, dest='pli', help='full production label of input', default=None)
  parser.add_argument('--dodefaultecaltags', dest='dodefaultecaltags', help='use default ECAL tags in GT, except for PFRH tag', action='store_true', default=False)
  
  return parser.parse_args()


if __name__ == "__main__":

  opt = getOptions()

  ##############################
  # job configurations
  #############################
  user = os.environ["USER"]
  evar = 'E' if opt.doflatenergy else 'Et'
  dumpcfg = 'Cfg_RecoSimDumper_gjets_cfg.py' if opt.ch=='gjetEM' else 'Cfg_RecoSimDumper_cfg.py'
  etRange='{}{}to{}GeV'.format(evar,opt.etmin,opt.etmax)
  if opt.pfrhmultbelow2p5 == 0:
    pfrhLabel= opt.pfrhmult if not opt.dorefpfrh else 'Ref'
  else:
    pfrhLabel = '{a}-{b}'.format(a=opt.pfrhmultbelow2p5, b=opt.pfrhmultabove2p5)
  if opt.seedmultbelow2p5 == 0:
    seedLabel= opt.seedmult if not opt.dorefseed else 'Ref'
  else:
    seedLabel = '{a}-{b}'.format(a=opt.seedmultbelow2p5, b=opt.seedmultabove2p5)
  #in case different values below/above 2.5 are not set, use the value of pfhrmult everywhere
  if opt.pfrhmultbelow2p5 == 0:
    pfrhMultbelow2p5 = opt.pfrhmult
    pfrhMultabove2p5 = opt.pfrhmult
  else:
    pfrhMultbelow2p5 = opt.pfrhmultbelow2p5
    pfrhMultabove2p5 = opt.pfrhmultabove2p5
  #in case different values below/above 2.5 are not set, use the value of seedmult everywhere
  if opt.seedmultbelow2p5 == 0:
    seedMultbelow2p5 = opt.seedmult
    seedMultabove2p5 = opt.seedmult
  else:
    seedMultbelow2p5 = opt.seedmultbelow2p5
    seedMultabove2p5 = opt.seedmultabove2p5
  thrLabelEB = 'RingEB' if opt.doringavgEB else 'XtalEB' 
  thrLabelEE = 'RingEE' if opt.doringavgEE else 'XtalEE'
  thrLabel = thrLabelEB + thrLabelEE
  safetyLabel = 'noMargin' if not opt.dosafetymargin else 'wMargin'
  if opt.ch != 'QCD':
    prodLabel='{c}_{e}_{g}_{d}_{pu}_thrsLumi{thl}_pfrh{pf}_seed{s}_{mr}_thr{thr}_shs{shs}_maxd{md}_y{y}_{v}_n{n}'.format(c=opt.ch,e=etRange,g=opt.geo,d=opt.det,pu=opt.pu,thl=opt.thrslumi,pf=pfrhLabel,s=seedLabel,mr=safetyLabel,thr=thrLabel,shs=opt.showersigmamult,md=opt.maxsigmadist,y=opt.year,v=opt.ver,n=opt.nevts)
  else:
    prodLabel='{c}_{pu}_thrsLumi{thl}_pfrh{pf}_seed{s}_thr{thr}_shs{shs}_maxd{md}_y{y}_{v}_n{n}'.format(c=opt.ch,e=etRange,g=opt.geo,d=opt.det,pu=opt.pu,thl=opt.thrslumi,pf=pfrhLabel,s=seedLabel,thr=thrLabel,shs=opt.showersigmamult,md=opt.maxsigmadist,y=opt.year,v=opt.ver,n=opt.nevts)
  prodLabelIn = prodLabel if opt.pli==None else opt.pli # prod label of input for postProduction , can be different from output production label
  dopu = 1 if opt.pu=='wPU' else 0
  # force the noise conditions to be of 2021 for reference thresholds
  if opt.dorefpfrh and opt.dorefseed:
    thrlumi = '180'
  else:
    thrlumi = opt.thrslumi
  # make so that '2021' (=TL180) and '2023' (=TL450) are supported
  if thrlumi == 2021:
    thrlumi = '180'
  elif thrlumi == 2023:
    thrlumi = '450'
  print 'lumiThrs: {a}'.format(a=thrlumi)  
  doringavgEB = 1 if opt.doringavgEB else 0
  doringavgEE = 1 if opt.doringavgEE else 0
  dosafetymargin = 1 if opt.dosafetymargin else 0
  dorefpfrh = 1 if opt.dorefpfrh else 0
  dorefseed = 1 if opt.dorefseed else 0
  if    (    opt.doringavgEB and not os.path.isfile('../../data/noise/PFRecHitThresholds_EB_ringaveraged_TL{}.txt'.format(thrlumi))) \
     or (    opt.doringavgEE and not os.path.isfile('../../data/noise/PFRecHitThresholds_EE_ringaveraged_TL{}.txt'.format(thrlumi))) \
     or (not opt.doringavgEB and not os.path.isfile('../../data/noise/PFRecHitThresholds_EB_TL{}.txt'.format(thrlumi))) \
     or (not opt.doringavgEE and not os.path.isfile('../../data/noise/PFRecHitThresholds_EE_TL{}.txt'.format(thrlumi))): 
  #   or (opt.dorefseed       and ( not os.path.isfile('fixed_SeedingThresholds_EB.txt')  or not os.path.isfile('fixed_SeedingThresholds_EE.txt') ) ) :
    raise RuntimeError('file with input thresholds not available, please check')
  doflatenergy = 1 if opt.doflatenergy else 0
  dodefaultecaltags = 1 if opt.dodefaultecaltags else 0
  nthr = 8 if opt.domultithread else 1
  if opt.domultijob and opt.njobs <= 1: raise RuntimeError('when running multiple jobs, the number of parallel jobs should be larger than 1')
  if opt.domultijob and opt.nevts % opt.njobs != 0 and not opt.dorecofromeos: raise RuntimeError('cannot split events in njobs evenly, please change njobs / nevts')
  njobs = opt.njobs if opt.domultijob else 1
  nevtsjob = opt.nevts if not opt.domultijob else opt.nevts/opt.njobs
  nevtspremixfile = 600 # current number of events in each premixed file
  npremixfiles = nevtsjob / nevtspremixfile + 1
  if opt.doreco and opt.custominput == None and opt.custominputdir == None: raise RuntimeError('you must supply the custom input, when running with doreco activated')
  if opt.doreco and not opt.dorecofromeos and not opt.domultijob and not os.path.isfile(opt.custominput): raise RuntimeError('custominput {} not found').format(opt.custominput)
  if opt.doreco and not opt.dorecofromeos and opt.domultijob and not os.path.isdir(opt.custominputdir): raise RuntimeError('custominputdir {} not found').format(opt.custominputdir)
  if opt.doreco and not opt.dorecofromeos and opt.domultijob: print 'A gentle reminder that if you are running with doreco and domultijob activated, you should use the same job splitting that was used for the original production'
  if not opt.doreco and opt.dorecofromeos: raise RuntimeError('You have to activate doreco option in order to run dorecofromeos')
  # special configurations for opt.dorecofromeos
  inseprefix = 'root://t3dcachedb.psi.ch:1094/'
  outseprefix = 'root://t3dcachedb.psi.ch:1094/'
  if opt.dorecofromeos and not opt.dodumperonly:
    inseprefix='root://eoscms.cern.ch/'
    # check directory exists
    if opt.custominputdir != None:
      command = 'xrdfs {sepx} ls {d}'.format(sepx=inseprefix, d=opt.custominputdir)
      out = subprocess.check_output(command, shell=True) # if not a dir, there will be an exception
      files=out.split('\n')[0:-1] # get the list of files
      if not files: raise RuntimeError('files not found')
      alljobids = map(lambda x: int(x.split('_job')[1].split('_step2')[0]), files) # isolate the job number
      njobs = max(alljobids)
      print 'Found njobs={} in {}'.format(njobs, opt.custominputdir)
      if len(alljobids)!=njobs: print 'Will try to run nj={nj}, but I already now that nf={nf} will fail'.format(nj=njobs,nf=njobs-len(alljobids))
      nevtsjob = -1 
    elif opt.custominput != None:
      command = 'xrdfs {sepx} stat {f}'.format(sepx=inseprefix, f=opt.custominput)
      out = subprocess.check_output(command, shell=True) # if not a file, an ERROR will be present in the out
      if 'ERROR' in out: 
        print command
        raise RuntimeError('file not found, exiting')

  if opt.docustomtime: 
    time1=opt.time.split(',')[0]
    time2=opt.time.split(',')[1]
    time3=opt.time.split(',')[2]
    time4=opt.time.split(',')[3]
    timed=opt.time.split(',')[4]
    times = [time1,time2,time3,time4,timed]

    sbatch_times = map(lambda x: '--time=0-{}:00'.format(x), times)
  else:
    if opt.domedium:
      time = '--time=1-23:59'
    elif opt.dolong:
      time = '--time=2-23:59'
    elif opt.doshort:
      time = '--time=0-02:00'
    else:
      time = '--time=1-00:00'
    sbatch_times = [time, time, time, time, '--time=0-01:00']


  ##############################
  # create production directory and logs directory within
  #############################
  prodDir = './{}'.format(prodLabel)
  command = 'mkdir -p {}'.format(prodDir)
  if not os.path.isdir(prodDir):
    os.system(command)
  #else: raise RuntimeError('directory {} already present, not going to overwrite'.format(prodDir))

  command = 'mkdir {}/logs'.format(prodDir)
  os.system(command)

  ############################
  # copy the relevant cmsDriver to prod directory
  ############################
  ## find the names first 
  if opt.ch != 'QCD': 
    step1_driverName = 'step1_{c}_{g}.py'.format(c=opt.ch,g=opt.geo)
  else:
    step1_driverName = 'step1_QCD.py'
  step2_driverName = 'step2_{pu}.py'.format(pu=opt.pu)
  step3_driverName = 'step3.py'
  step4_driverName = 'step4.py'
  drivers = [step1_driverName, step2_driverName, step3_driverName]
  if opt.doold: drivers=map(lambda x : 'old_' + x, drivers)
  target_drivers = ['step1.py', 'step2.py', 'step3.py']
  infiles  = ['', 'step1_nj{nj}.root', 'step2_nj{nj}.root']
  infiles_loc = ['', 'step1.root', 'step2.root', 'step3.root']
  outfiles = ['step1_nj{nj}.root', 'step2_nj{nj}.root', 'step3_nj{nj}.root']
  outfiles_loc = ['step1.root', 'step2.root', 'step3.root']

  if opt.dorecofromeos: 
    infiles = ['', '', 'cluster_job{nj}_step2.root']
    
  if opt.ch == 'QCD':
   drivers.append(step4_driverName)
   target_drivers.append('step4.py')
   infiles.append('step3_nj{nj}.root')
   infiles_loc.append('step3.root')
   outfiles.append('step4_nj{nj}.root')
   outfiles_loc.append('step4.root')

  ## copy them to dir
  for i,idriver in enumerate(drivers):
    if opt.doreco and i<2: continue # skip everything that is not related to step3
    if not os.path.isfile('cmsDrivers/{idr}'.format(idr=idriver)):
      raise RuntimeError('cmsDriver {idr} not found, please check naming'.format(idr=idriver))
    command = 'cp cmsDrivers/{idr} {d}/{td}'.format(idr=idriver,d=prodDir,td=target_drivers[i])
    os.system(command) 

  ## also copy the file containing the conditions
  tag_filename = 'override_ECAL_tags.py'
  command = 'cp cmsDrivers/{tf} {d}/.'.format(tf=tag_filename,d=prodDir)
  os.system(command)

  ############################
  # write the cmsRun commands for all steps
  ############################
  ## step1
  if opt.geo == 'closeEcal' and opt.ch != 'QCD':
    if opt.det == 'EB':
      rmin = 123.8
      rmax = 123.8
      zmin = -304.5
      zmax = 304.5
      npart = 10
    elif opt.det == 'EEclose':
      #rmin = 58.0 # eta=2.0
      rmin = 71.1 # eta=2.2
      #rmin = 87.4 # eta=2.4
      rmax = 171.1
      zmin = 317.0
      zmax = 317.0
      npart = 10
    elif opt.det == 'EEfar':
      rmin = 31.6
      #rmax = 58.0 # eta=2.0
      rmax = 71.1 # eta=2.2
      #rmax = 87.4 # eta=2.4
      rmax = 87.4
      zmin = 317.0
      zmax = 317.0
      npart = 10
  
    if opt.npart!=None:
      npart = opt.npart

    step1_cmsRun = 'cmsRun {jo} maxEvents={n} etmin={etmin} etmax={etmax} rmin={r1} rmax={r2} zmin={z1} zmax={z2} np={np} nThr={nt} doFlatEnergy={dfe} year={y} doDefaultECALtags={ddet}'.format(jo=target_drivers[0], n=nevtsjob, etmin=float(opt.etmin), etmax=float(opt.etmax), r1=rmin, r2=rmax, z1=zmin, z2=zmax, np=npart, nt=nthr, dfe=doflatenergy, y=opt.year, ddet=dodefaultecaltags)
    step1_cmsRun_add = 'seedOffset={nj}' # format at a later stage
  elif opt.ch == 'QCD':
    step1_cmsRun = 'cmsRun {jo} maxEvents={n} nThr={nt} year={y} doDefaultECALtags={ddet}'.format(jo=target_drivers[0], n=nevtsjob, nt=nthr, y=opt.year, ddet=dodefaultecaltags)
    step1_cmsRun_add = 'seedOffset={nj}' # format at a later stage
  elif opt.doreco:
    step1_cmsRun = 'dummy'
    step1_cmsRun_add = 'dummy'
  else:
    raise RuntimeError('this option is not currently supported')
  ## step2  
  step2_cmsRun = 'cmsRun {jo} nThr={nt} year={y} doDefaultECALtags={ddet}'.format(jo=target_drivers[1], nt=nthr, y=opt.year, ddet=dodefaultecaltags)
  step2_cmsRun_add = ('nPremixFiles={npf}'.format(npf=npremixfiles) if dopu else '') + (' randomizePremix=1' if opt.domultijob and dopu else ' ')
  ## step3
  if opt.ch != 'QCD':
    dorecofile=1
    dominiaodfile=0
  else:
    dorecofile=0
    dominiaodfile=1
  step3_cmsRun = 'cmsRun {jo} dorecofile={reco} dominiaodfile={miniaod} pfrhMultbelow2p5={pfrhmb} pfrhMultabove2p5={pfrhma} seedMultbelow2p5={smb} seedMultabove2p5={sma} nThr={nt} thrsLumi={thl} doRefPfrh={drpf} doRefSeed={drsd} doSafetyMargin={mr} doPU={dp} doRingAverageEB={draeb} doRingAverageEE={draee} year={y} doDefaultECALtags={ddet} showerSigmaMult={shs} maxSigmaDist={md} maxEvents={n}'.format(jo=target_drivers[2], reco=dorecofile, miniaod=dominiaodfile, pfrhmb=pfrhMultbelow2p5, pfrhma=pfrhMultabove2p5, smb=seedMultbelow2p5, sma=seedMultabove2p5, nt=nthr, thl=thrlumi, drpf=dorefpfrh, drsd=dorefseed, mr=dosafetymargin, dp=dopu, draeb=doringavgEB, draee=doringavgEE, y=opt.year, ddet=dodefaultecaltags, shs=opt.showersigmamult, md=opt.maxsigmadist, n=nevtsjob)
  cmsRuns = [step1_cmsRun, step2_cmsRun, step3_cmsRun]
  cmsRuns_add = [step1_cmsRun_add, step2_cmsRun_add, '']
  if opt.ch == 'QCD':
    step4_cmsRun = 'cmsRun {jo} maxEvents={n} nThr={nt} year={y} doDefaultECALtags={ddet}'.format(jo=target_drivers[3], n=nevtsjob, nt=nthr, y=opt.year, ddet=dodefaultecaltags) 
    cmsRuns.append(step4_cmsRun)
    cmsRuns_add.append('')  
  
  
  if not opt.dodumperonly:
  ############################
  # write the launching scripts
  ############################
    for i,idriver in enumerate(drivers):

      if opt.doreco and i<2: continue # skip everything that is not related to step3

      for nj in range(0,njobs):
    
        ### configurations for the template script
        outputDir = '"/pnfs/psi.ch/cms/trivcat/store/user/"$USER"/EcalProd/"' 
        if opt.dosavehome: outputDir = '`pwd`/../' 

        mkdiroutput_command = 'xrdfs t3dcachedb03.psi.ch mkdir $SERESULTDIR'
        if opt.dosavehome: mkdiroutput_command = 'mkdir -p $SERESULTDIR'

        cpinput_command = ''
        if infiles[i]!='':
          cpinput_command = 'xrdcp $INSEPREFIX/$SERESULTDIR/{infile} $WORKDIR/{infile_loc}'.format(infile=infiles[i].format(nj=nj),infile_loc=infiles_loc[i])
          if opt.dosavehome: cpinput_command = 'cp $SERESULTDIR/{infile} $WORKDIR/{infile_loc}'.format(infile=infiles[i].format(nj=nj),infile_loc=infiles_loc[i])
          if opt.doreco: 
            custominput = opt.custominput if not opt.domultijob else opt.custominputdir + '/' + infiles[i].format(nj=nj) 
            if opt.dosavehome: cpinput_command = 'cp {ci} $WORKDIR/{infile_loc}'.format(ci=custominput,infile_loc=infiles_loc[i])
            else:              
              # in case of step4, the step3 has to be taken from the workdir and not from the custominput 
              if 'step4' not in idriver:
                cpinput_command = 'xrdcp $INSEPREFIX/{ci} $WORKDIR/{infile_loc}'.format(ci=custominput,infile_loc=infiles_loc[i])
              else: 
                cpinput_command = 'xrdcp $INSEPREFIX/$SERESULTDIR/{infile} $WORKDIR/{infile_loc}'.format(infile=infiles[i].format(nj=nj), infile_loc=infiles_loc[i])
      
              
        cpoutput_command = 'xrdcp -f {outfile_loc} $OUTSEPREFIX/$SERESULTDIR/{outfile}'.format(outfile_loc=outfiles_loc[i],outfile=outfiles[i].format(nj=nj))
        if opt.dosavehome: cpoutput_command = 'cp {outfile_loc} $SERESULTDIR/{outfile}'.format(outfile_loc=outfiles_loc[i],outfile=outfiles[i].format(nj=nj))
         
        cpaux_command = ''
        if 'step3' in idriver:
          cpaux_command = 'cp -r $CMSSW_BASE/src/RecoSimStudies/Dumpers/data $WORKDIR'

        ### define a template script  
        template = [
        '#!/bin/bash',
        '',

        #### variables
        'DIRNAME="{ind}"',
        'STARTDIR=`pwd`',
        'TOPWORKDIR="/scratch/"$USER/',
        'JOBDIR="gen_"$SLURM_JOB_ID',
        'WORKDIR=$TOPWORKDIR/$JOBDIR',
        'INSEPREFIX="{isepx}"',
        'OUTSEPREFIX="{osepx}"',
        'SERESULTDIR={od}/$DIRNAME',
        'JOBOPFILENAME="{jo}"',
        '',

        #### environment
        'source $VO_CMS_SW_DIR/cmsset_default.sh',
        'shopt -s expand_aliases',
        'echo ""',
        'echo "Going to set up cms environment"',
        'cd $STARTDIR',
        'cmsenv',
        'echo ""',
        '',

        #### workdir
        'echo "Going to create work dir"',
        'mkdir -p $WORKDIR',
        'echo "workdir: "',
        'echo $WORKDIR',
        'echo ""',
        '',

        #### outputdir
        'echo "Going to create the output dir"',
        'echo "May give an error if the directory already exists, which can be safely ignored"',
        '{mkdir}',
        'echo ""',
        '',

        #### copy driver and other aux files 
        'echo "Going to copy cms driver and aux files"',
        'cp $JOBOPFILENAME $WORKDIR/$JOBOPFILENAME',
        '{cpaux}',
        'cp {tf} $WORKDIR/.',
        'echo ""',
        '',

        #### copy input file 
        'echo "Going to copy input file if needed"',
        '{cpin}',
        'echo ""',
        '',

        #### run
        'cd $WORKDIR',
        'echo ""',
        '',
        'echo "Going to run"',
        'DATE_START=`date +%s`',
        '{cmsRun}',
        'DATE_END=`date +%s`',
        'echo ""',
        '',
        'echo "Finished running"',
        'echo "Content of current directory"',
        'ls -al',

        #### copy back output
        'echo "Going to copy the output to the output directory"',
        '{cpout}',
        '',
        'echo ""',

        #### clean and go 
        'echo "Cleaning up $WORKDIR"',
        'rm -rf $WORKDIR',
        'RUNTIME=$((DATE_END-DATE_START))',
        'echo "Wallclock running time: $RUNTIME s"',
        'cd $STARTDIR',

        ] 
        template = '\n'.join(template)
        template = template.format(ind=prodLabel,od=outputDir,jo=target_drivers[i],mkdir=mkdiroutput_command,isepx=inseprefix,osepx=outseprefix,
                                   #cpin=cpinput_command,cpout=cpoutput_command,cmsRun=cmsRuns[i]+' '+cmsRuns_add[i].format(nj=nj+nthr+1),cpaux=cpaux_command) 
                                   cpin=cpinput_command,cpout=cpoutput_command,cmsRun=cmsRuns[i]+' '+cmsRuns_add[i].format(nj=nj+1),cpaux=cpaux_command,tf=tag_filename) 

        launcherFile = '{}/launch_{}.sh'.format(prodDir,outfiles[i].format(nj=nj).split('.root')[0])
        with open(launcherFile, 'w') as f:
          f.write(template)

    #############################
    # write the template script to run the postproduction helper (one job for the full task)
    ############################
  if not opt.doskipdumper and opt.ch != 'QCD':
 
    # the template should run a python script to get the sample list
    # all done locally
    template_postprod = [
      '#!/bin/bash',
      '',
      'STARTDIR=$PWD/../',
      #### environment
      'source $VO_CMS_SW_DIR/cmsset_default.sh',
      'shopt -s expand_aliases',
      'echo ""',
      'echo "Going to set up cms environment"',
      'cd $STARTDIR',
      'cmsenv',
      'echo ""',
      '',
      ### running part
      'echo "Going to run postProdHelper.py"',
      'DATE_START=`date +%s`',
      'python postProdHelper.py --pli {pli} --plo {plo} --user {u} --splitfactord {spd}'.format(pli=prodLabelIn,plo=prodLabel, u=user, spd=opt.splitfactord),
      'echo ""',
      'DATE_END=`date +%s`',
      'echo ""',
      '',
      ### print job time
      'echo "Finished running"',
      'RUNTIME=$((DATE_END-DATE_START))',
      'echo "Wallclock running time: $RUNTIME s"',
    ]
    
    template_postprod = '\n'.join(template_postprod)

    launcherFile_postprod = '{}/launch_postprod.sh'.format(prodDir)
    with open(launcherFile_postprod, 'w') as f:
      f.write(template_postprod)

    #############################
    # write the template for the dumper
    ############################
    
    for njd in range(0,opt.splitfactord):
      template_dumper = [
        '#!/bin/bash',
        '',
        'STARTDIR=$PWD/../',
        'TOPWORKDIR="/scratch/"$USER/',
        'JOBDIR="dump_"$SLURM_JOB_ID',
        'WORKDIR=$TOPWORKDIR/$JOBDIR',
        'INSEPREFIX="{isepx}"',
        'OUTSEPREFIX="{osepx}"',
        'SERESULTDIR=$OUTSEPREFIX//pnfs/psi.ch/cms/trivcat/store/user/$USER/Dumper/',        
        '',

        #### environment
        'source $VO_CMS_SW_DIR/cmsset_default.sh',
        'shopt -s expand_aliases',
        'echo ""',
        'echo "Going to set up cms environment"',
        'cd $STARTDIR',
        'cmsenv',
        'echo ""',
        '',

        #### workdir
        'echo "Going to create work dir, will contain intermediate output"',
        'mkdir -p $WORKDIR',
         'echo "workdir: "',
        'echo $WORKDIR',
        'echo ""',
        '',

        ### running part
        'echo "Going to run the Dumper"',
        'DATE_START=`date +%s`',
        'if [ "$USER" == "anlyon" ] ; then ',
          'cmsRun ../../python/{cfg} outputFile=/work/anlyon/dumpedFiles/{pl}_njd{njd}.root inputFiles_load=../../data/samples/{pl}_njd{njd}.txt',
        'else',
          'cmsRun ../../python/{cfg} outputFile=$WORKDIR/{pl}_njd{njd}.root inputFiles_load=../../data/samples/{pl}_njd{njd}.txt',
        'fi',
        'echo ""',
        '',

        ### copy back output
        'echo "Going to copy the output to the output directory"',
        'xrdcp -f $WORKDIR/{pl}_njd{njd}.root $SERESULTDIR/{pl}_njd{njd}.root',
        'DATE_END=`date +%s`',
        'echo ""',

        #### clean and go 
        'echo "Cleaning up $WORKDIR"',
        'rm -rf $WORKDIR',
        'RUNTIME=$((DATE_END-DATE_START))',
        'echo "Wallclock running time: $RUNTIME s"',
        'cd $STARTDIR',
      ]
      
      template_dumper = '\n'.join(template_dumper)
      template_dumper = template_dumper.format(pl=prodLabel, njd=njd, isepx=inseprefix, osepx=outseprefix, cfg=dumpcfg)

      launcherFile_dumper = '{}/launch_dumper_njd{}.sh'.format(prodDir,njd)
      with open(launcherFile_dumper, 'w') as f:
        f.write(template_dumper)


  #############################
  # finally write the submitter
  ############################# 
  submitter_template = []
  postprod_dependencies = '' # 

  if not opt.dodumperonly:
    postprod_dependencies = '--dependency=afterany'
    for nj in range(0,njobs):

      sbatch_command_step1 = 'jid1_nj{nj}=$(sbatch -p wn --account=t3 -o logs/step1_nj{nj}.log -e logs/step1_nj{nj}.log --job-name=step1_{pl} {t} --ntasks={nt} --mem 2500 launch_step1_nj{nj}.sh)'.format(nj=nj,pl=prodLabel,t=sbatch_times[0],nt=nthr)

      sbatch_command_step2 = 'jid2_nj{nj}=$(sbatch -p wn --account=t3 -o logs/step2_nj{nj}.log -e logs/step2_nj{nj}.log --job-name=step2_{pl} {t} --ntasks={nt} --mem 2000 --dependency=afterany:$jid1_nj{nj} launch_step2_nj{nj}.sh)'.format(nj=nj,pl=prodLabel,t=sbatch_times[1],nt=nthr)

      sbatch_command_step3 = 'jid3_nj{nj}=$(sbatch -p wn --account=t3 -o logs/step3_nj{nj}.log -e logs/step3_nj{nj}.log --job-name=step3_{pl} {t} --ntasks={nt} --mem 4000 --dependency=afterany:$jid2_nj{nj} launch_step3_nj{nj}.sh)'.format(nj=nj,pl=prodLabel,t=sbatch_times[2],nt=nthr)
      if opt.doreco: # strip the dependency away
        sbatch_command_step3 = 'jid3_nj{nj}=$(sbatch -p wn --account=t3 -o logs/step3_nj{nj}.log -e logs/step3_nj{nj}.log --job-name=step3_{pl} {t} --ntasks={nt} --mem 4000 launch_step3_nj{nj}.sh)'.format(nj=nj,pl=prodLabel,t=sbatch_times[2],nt=nthr)
   
      if opt.ch == 'QCD':  
        sbatch_command_step4 = 'jid4_nj{nj}=$(sbatch -p wn --account=t3 -o logs/step4_nj{nj}.log -e logs/step4_nj{nj}.log --job-name=step4_{pl} {t} --ntasks={nt} --mem 4000 --dependency=afterany:$jid3_nj{nj} launch_step4_nj{nj}.sh)'.format(nj=nj,pl=prodLabel,t=sbatch_times[3],nt=nthr)

      postprod_dependencies += ':$jid3_nj{nj}'.format(nj=nj)

      if not opt.doreco:
        submitter_template.append(sbatch_command_step1)
        submitter_template.append('echo "$jid1_nj%i"' % nj)
        submitter_template.append('jid1_nj%i=${jid1_nj%i#"Submitted batch job "}' % (nj,nj))
        submitter_template.append(sbatch_command_step2)
        submitter_template.append('echo "$jid2_nj%i"' % nj)
        submitter_template.append('jid2_nj%i=${jid2_nj%i#"Submitted batch job "}' % (nj,nj))

      submitter_template.append(sbatch_command_step3)
      submitter_template.append('echo "$jid3_nj%i"' % nj)
      submitter_template.append('jid3_nj%i=${jid3_nj%i#"Submitted batch job "}' % (nj,nj))

      if opt.ch == 'QCD':
        submitter_template.append(sbatch_command_step4)
        submitter_template.append('echo "$jid4_nj%i"' % nj)
        submitter_template.append('jid4_nj%i=${jid4_nj%i#"Submitted batch job "}' % (nj,nj))
      
  # add the postproduction and dumper part
  if not opt.doskipdumper and opt.ch != 'QCD':
    # post production
    sbatch_command_postprod = 'jid_pp=$(sbatch -p wn --account=t3 -o logs/postprod.log -e logs/postprod.log --job-name=postprod_{pl} {t} --ntasks=1 {dd} launch_postprod.sh)'.format(pl=prodLabel,t='--time=0-01:00',dd=postprod_dependencies)
    submitter_template.append(sbatch_command_postprod)
    submitter_template.append('echo "$jid_pp"')
    submitter_template.append('jid_pp=${jid_pp#"Submitted batch job "}')
    dumper_dependencies = ':$jid_pp'
    # dumper
    for njd in range(0, opt.splitfactord):
      sbatch_command_dumper = 'jid_njd{njd}=$(sbatch -p wn --account=t3 -o logs/dumper_njd{njd}.log -e logs/dumper_njd{njd}.log --job-name=dumper_{pl} {t} --ntasks=1 --dependency=afterany{dd} launch_dumper_njd{njd}.sh)'.format(njd=njd,pl=prodLabel,t=sbatch_times[4],dd=dumper_dependencies)
      submitter_template.append(sbatch_command_dumper)
      submitter_template.append('echo "$jid_njd%i"' % njd)
      submitter_template.append('jid_njd%i=${jid_njd%i#"Submitted batch job "}' % (njd,njd))
   

  # finish ...
  submitter_template = '\n\n'.join(submitter_template)

  submitterFile = '{}/submit.sh'.format(prodDir)
  with open(submitterFile, 'w') as f:
    f.write(submitter_template)


