'''
Helper for resubmission of step3 (and dumper)
Creates a "resubmit.sh" including only the jobs that failed
'''

import sys
import os
import glob
import re
from ROOT import TNetXNGFile

def getOptions():
  from argparse import ArgumentParser

  parser = ArgumentParser(description='Production helper for Clustering studies', add_help=True)

  parser.add_argument('--custominputdir', type=str, dest='custominputdir', help='SE path of the directory where previous step was run, if different from that of production label', default=None)
  parser.add_argument('--pl', type=str, dest='pl', help='production label', default='photon_bla')
  parser.add_argument('--user', type=str, dest='user', help='user that produced the sample', default='mratti', choices=['anlyon', 'mratti'])
  
  return parser.parse_args()


if __name__ == "__main__":
  
  opt = getOptions()
  pldir = '/pnfs/psi.ch/cms/trivcat/store/user/{u}/EcalProd/{pl}/'.format(u=opt.user,pl=opt.pl)
  if opt.custominputdir == None: 
    opt.custominputdir = pldir
  if not os.path.isdir(pldir):
    raise RuntimeError('Production label directory not valid, please check, {}'.format(opt.pldir))
  if not os.path.isdir(opt.custominputdir):
    raise RuntimeError('Custom input dir not valid, please check, {}'.format(opt.custominputdir))

  custominputexpr = opt.custominputdir + '/step2_nj*root' 
  plexpr = pldir + '/step3_nj*root'
  prepend = 'root://t3dcachedb.psi.ch:1094/'
  #samplefile = '../../data/samples/' + opt.pl + '.txt'

  # which njX you would expect from step2
  nj_expected = []
  step2_files = [f for f in glob.glob(custominputexpr)]
  step2_files = map(lambda x: prepend+x, step2_files)  
  for f in step2_files:
    rf = TNetXNGFile.Open(f, 'r')
    if rf and rf.GetListOfKeys().Contains('Events'):
      #print 'valid step2 file, ', f.split('nj')[1].split('.root')[0]
      nj_expected.append(f.split('nj')[1].split('.root')[0])

  # which njY you would expect from step3
  nj_actual = []
  step3_files = [f for f in glob.glob(plexpr)]
  step3_files = map(lambda x: prepend+x, step3_files)
  for f in step3_files:
    rf = TNetXNGFile.Open(f, 'r')
    if rf and rf.GetListOfKeys().Contains('Events'):
      nj_actual.append(f.split('nj')[1].split('.root')[0])

  set_expected = set(nj_expected)
  set_actual   = set(nj_actual)
  if not set_actual.issubset(set_expected): raise RuntimeError('More valid files at previous step, logic error')
  set_missing = set_expected.difference(set_actual)

  # now open the submission file and create another one where only the missing jobs are to be resubmitted 
  resubmit_file = open('./{}/resubmit.sh'.format(opt.pl), 'w')
  submit_file = open('./{}/submit.sh'.format(opt.pl), 'r')

  for line in submit_file:

    if 'jid3' in line:  # only consider resubmission of step3

      result = re.search('nj(.*)\=\$',line)
      if result and result.group(1) in set_missing:
        if 'sbatch' in line:
          start=line.split(' --ntasks=1')[0].replace('cn-test', 't3') + ' --ntasks=1 '
          end=' launch'+line.split('launch')[1]
          resubmit_file.write(start+end+'\n')
        else: 
          resubmit_file.write(line+'\n')
      result = re.search('nj(.*)\"', line)
      if result and result.group(1) in set_missing:
        resubmit_file.write(line+'\n')

  # also write the dumper part
  sub1 = 'jid_d=$(sbatch -p wn --account=t3 -o logs/dumper.log -e logs/dumper.log --job-name=dumper_photon_E1.0to200GeV_closeEcal_EEfar_wPU_pfrhRef_seedRef_thrXtalEBXtalEE_y2023_T2_v6_t0_n30000 --time=0-03:59 --ntasks=1 --dependency=afterany'
  sub2 = ''.join(map(lambda x: ':$jid3_nj{}'.format(x), list(set_missing)))
  sub3 =  ' launch_dumper.sh)'

  resubmit_file.write(sub1+sub2+sub3+'\n')
  resubmit_file.write('echo "$jid_d"\n')
  resubmit_file.write('jid_d=${jid_d#"Submitted batch job "}\n')

  # close and goodbye
  resubmit_file.close()
  submit_file.close()

  if os.path.isfile('./{}/resubmit.sh'.format(opt.pl)):

    print ''
    print 'Wrote resubmission script'
    print './{}/resubmit.sh'.format(opt.pl)
    print ''
  
