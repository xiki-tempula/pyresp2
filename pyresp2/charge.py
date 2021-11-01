import os
from glob import glob
from os.path import join
from pkg_resources import resource_filename
import MDAnalysis as mda
import numpy as np

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

    with open(resource_filename(__name__, 'data/RESP.in'), 'r') as f:
        Multiwfn_in = f.read()

    with open('charge/Multiwfn_resp_HF.in', 'w') as f:
        f.write(Multiwfn_in.replace('conflist.txt', 'conflist_HF.txt'))

    with open('charge/Multiwfn_resp_Vac.in', 'w') as f:
        f.write(Multiwfn_in.replace('conflist.txt', 'conflist_Vac.txt'))

    with open('charge/Multiwfn_resp_Sol.in', 'w') as f:
        f.write(Multiwfn_in.replace('conflist.txt', 'conflist_Sol.txt'))

def _find_RESP(text):
    with open(text, 'r') as f:
        lines = f.read().split('\n')
    charge_list = []
    start = False
    for line in lines:
        if 'Center       Charge' in line:
            start = True
            continue
        if 'Sum of charges' in line:
            return np.array(charge_list)
        if start:
            charge_list.append(float(line.split()[-1]))

def read_charge(mol, multiwfn_log):
    u = mda.Universe(mol)
    charges = _find_RESP(multiwfn_log)
    assert len(charges) == len(u.atoms)
    u.atoms.charges = charges
    return u

def read_all_charges(mol, HF='charge/Multiwfn_resp_HF.log',
                     Vac='charge/Multiwfn_resp_Vac.log',
                     Sol='charge/Multiwfn_resp_Sol.log'):
    HF_mol = read_charge(mol, HF)
    HF_mol.atoms.write('charge/HF.mol2')
    Vac_mol = read_charge(mol, Vac)
    Sol_mol = read_charge(mol, Sol)
    u = mda.Universe(mol)
    u.atoms.charges = (Vac_mol.atoms.charges + Sol_mol.atoms.charges) / 2
    u.atoms.write('charge/RESP2.mol2')






