import os
from glob import glob
from os.path import join
import numpy as np
import shutil

def gen_charge(config):
    file_list = glob(join('charge', '*.xyz'))
    for file in file_list:
        base, ext = os.path.splitext(file)
        with open(file, 'r') as f:
            crd = f.read().split('\n')[2:]

        with open('{}.gjf'.format(base), 'w') as f:
            f.write('''%chk=HF_{i}.chk
%mem={mem}MB
%NProcShared={n_proc}
# {resp}

RESP

{charge} {multiplicity}
{crd}

--link1--
%oldchk=HF_{i}.chk
%chk=Vac_{i}.chk
%mem={mem}MB
%NProcShared={n_proc}
# Geom=AllCheck Guess=Read {resp2_vac}

--link1--
%oldchk=Vac_{i}.chk
%chk=Sol_{i}.chk
%mem={mem}MB
%NProcShared={n_proc}
# Geom=AllCheck Guess=Read {resp2_sol}
'''.format(mem=config['charge']['memory'],
           n_proc=config['charge']['n_proc'],
           charge=config['molecule']['charge'],
           multiplicity=config['molecule']['multiplicity'],
           resp=config['charge']['resp'],
           resp2_vac=config['charge']['resp2_vac'],
           resp2_sol=config['charge']['resp2_sol'],
           i=os.path.basename(base).split('_')[-1],
           crd='\n'.join(crd)))




