# https://docs.google.com/spreadsheets/d/1rJ4FppHxdD6TQpCHnnnX0QD0Y0d4aTEN6ibgTVn7Vyg/edit#gid=0

import sys
import os
import glob
import ROOT

def getOptions():
  from argparse import ArgumentParser

  parser = ArgumentParser(description='Production helper for Clustering studies', add_help=True)

  parser.add_argument('--pl', type=str, dest='pl', help='full production label', default='photon_bla')
  parser.add_argument('--user', type=str, dest='user', help='user that produced the sample', default='mratti', choices=['anlyon', 'mratti'])
  
  return parser.parse_args()


if __name__ == "__main__":
  
  opt = getOptions()

  # input 
  path = '/pnfs/psi.ch/cms/trivcat/store/user/{u}/EcalProd/{pl}/'.format(u=opt.user,pl=opt.pl)
  expr = 'step3_nj*root'
  # output
  prepend = 'root://t3dcachedb.psi.ch:1094/'
  samplefile = '../../data/samples/' + opt.pl + '.txt'

  files = [f for f in glob.glob(path + expr)]
  if len(files)==0: raise RuntimeError('path {} does not exist or is empty'.format(path))

  files = map(lambda x: prepend+x, files) 

  valid_files = []
  # check files are not empty
  for f in files:
    rf = ROOT.TNetXNGFile.Open(f, 'r')
    if rf and rf.GetListOfKeys().Contains('Events'):
      #print 'file usable'
      valid_files.append(f)
    else: 
      print 'file not usable, will continue'
      
  print valid_files

  #create outputdir
  os.system('mkdir {}'.format('../outputfiles/dumpedFiles'))

  # write sample file with all files per production
  with open(samplefile, 'w') as of:
    of.write('\n'.join(valid_files))

  print ''
  print 'Wrote list of files for production', opt.pl
  print ''
  
  # write cmsRun command to launch
  command = 'cmsRun python/Cfg_RecoSimDumper_cfg.py outputFile=test/outputfiles/dumpedFiles/{pl}.root inputFiles_load=data/samples/{pl}.txt'.format(pl=opt.pl)
    
  print 'Command to run for production', opt.pl
  print ''
  print command
  print ''

  
