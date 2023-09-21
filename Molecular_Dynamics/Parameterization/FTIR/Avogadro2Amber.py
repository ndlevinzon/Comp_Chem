#!/usr/bin/python
import sys, os, csv, subprocess
#fh = open('decaC.txt', 'rU')
fh = open('names.dsRNA.txt', 'rU')
for line in fh:
  molecule='%s' %(line.split()[0].strip("").split('.')[0])
  os.system('sed "s/\ DA/\ A\ /g" %s.orig.pdb > %s_.pdb'%(molecule,molecule))
  os.system('sed -i "s/\ DT/\ U\ /g" %s_.pdb'%(molecule))
  os.system('sed -i "s/\ DC/\ C\ /g" %s_.pdb'%(molecule))
  os.system('sed -i "s/\ DG/\ G\ /g" %s_.pdb'%(molecule))
  myoutfile = '%s.pdb' %molecule
  myinfile = '%s_.pdb' %molecule
  outf=open(myoutfile,'w')
  fh2= open(myinfile, 'rU') 
  entry2=''
  for line in fh2:
    entry0 = line.split()[0]
    if (len(line.split())>2):
      entry2 = line.split()[2]
##    if (entry0!='CONECT'):
## MUST MAUNALLY DELETE C5 from uracil
    if (entry0!='CONECT' and entry0!='MASTER' and entry2!='C5M' and entry2!='OXT' and entry2!='P' and entry2!='O1P' and entry2!='O2P' and entry2!='H51' and entry2!='H52' and entry2!='H53' and entry2!='HTER' and entry2!='HCAP' and entry2!='H2\'2'):
      outf.write(line)
  outf.close()
  fh2.close()
fh.close()
