# https://docs.google.com/spreadsheets/d/1rJ4FppHxdD6TQpCHnnnX0QD0Y0d4aTEN6ibgTVn7Vyg/edit#gid=0

import sys
import os
import glob
import ROOT

def getOptions():
  from argparse import ArgumentParser

  parser = ArgumentParser(description='Production helper for Clustering studies', add_help=True)

  parser.add_argument('--pli', type=str, dest='pli', help='full production label of input', default='photon_bla')
  parser.add_argument('--plo', type=str, dest='plo', help='full production label of output', default='photon_bla')
  parser.add_argument('--user', type=str, dest='user', help='user that produced the sample', default='mratti', choices=['anlyon', 'mratti'])
  parser.add_argument('--splitfactord', type=int, dest='splitfactord', help='number of jobs for the dumper', default=1)
  
  return parser.parse_args()


if __name__ == "__main__":
  
  opt = getOptions()
  # input 
  path = '/pnfs/psi.ch/cms/trivcat/store/user/{u}/EcalProd/{pl}/'.format(u=opt.user,pl=opt.pli)
  expr = 'step3_nj*root'
  # output
  prepend = 'root://t3dcachedb.psi.ch:1094/'
  samplefile = '../../data/samples/' + opt.plo + '.txt'

  files = [f for f in glob.glob(path + expr)]
  files.sort() # sort in alphabetical order
  if len(files)==0: raise RuntimeError('path {} does not exist or is empty'.format(path))

  files = map(lambda x: prepend+x, files) 

  valid_files = []
  # check files are not empty
  for f in files:
    #print f
    rf = ROOT.TNetXNGFile.Open(f, 'r')
    if rf and rf.GetListOfKeys().Contains('Events'):
      #print 'file usable'
      valid_files.append(f)
    else: 
      print 'file not usable, will continue'
  
  print '{} usable file / {} total files = {}'.format(len(valid_files), len(files), float(len(valid_files))/float(len(files)))

  if len(valid_files) < opt.splitfactord:
    #raise RuntimeError('Invalid splitfactor, please check')
    print 'WARNING, invalid split factor, setting it to 1'
    opt.splitfactord=1
      
  #create outputdir
  os.system('mkdir {}'.format('../outputfiles/dumpedFiles'))

  # now divide based on splitfactor
  n_files_perjob = len(valid_files) / opt.splitfactord # division between two integers is an integer
  if len(valid_files) % opt.splitfactord != 0: print('N={} files will not be run'.format(len(valid_files) % opt.splitfactord))
  groups_validfiles = [valid_files[i:i+n_files_perjob] for i in range(0,len(valid_files),n_files_perjob)]
  #print 'groups_validfiles, len=', len(groups_validfiles)
  #for i in range(0,len(groups_validfiles)):
    #print 'i=', groups_validfiles[i], '\n'
  
  samplefiles = [samplefile.replace('.txt', '_njd{}.txt'.format(i)) for i in range(0,opt.splitfactord)]
  print 'samplefiles, len=', len(samplefiles)
  print samplefiles
   
  # write sample files with files needed for production
  for njd in range(0,opt.splitfactord):
    with open(samplefiles[njd], 'w') as of:
      of.write('\n'.join(groups_validfiles[njd]))
    print ''
    print 'Wrote list of files {}'.format(samplefiles[njd])
    print ''

    command = 'cmsRun python/Cfg_RecoSimDumper_cfg.py outputFile=test/outputfiles/dumpedFiles/{pl}_njd{njd}.root inputFiles_load=data/samples/{pl}_njd{njd}.txt'.format(pl=opt.plo, njd=njd)
    
    print 'Command to run',
    print ''
    print command
    print ''

#  # write sample file with all files per production
#  with open(samplefile, 'w') as of:
#    of.write('\n'.join(valid_files))
#
#  print ''
#  print 'Wrote list of files for production', opt.pl
#  print ''
  
#  # write cmsRun command to launch
#  command = 'cmsRun python/Cfg_RecoSimDumper_cfg.py outputFile=test/outputfiles/dumpedFiles/{pl}.root inputFiles_load=data/samples/{pl}.txt'.format(pl=opt.pl)
#    
#  print 'Command to run for production', opt.pl
#  print ''
#  print command
#  print ''

  
