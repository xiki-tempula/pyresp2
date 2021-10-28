import os
import shutil
import subprocess
import MDAnalysis
import numpy as np
import shutil

def _check_charge(mol, charge):
    # Check the total charge of the molecule.
    # If the charge is provided in the PDB file, xtb will prioritise the PDB
    # charge.
    u = MDAnalysis.Universe(mol)
    try:
        total_charge = sum(u.atoms.charges)
        if not np.isclose(total_charge, 0, atol=0.01) and \
                not np.isclose(total_charge, charge, atol=0.01):
            raise ValueError('The total charge of the input coordinates ({}) '
                             'is different from the total charge provided '
                             '({}).'.format(total_charge, charge))
    except MDAnalysis.exceptions.NoDataError:
        pass

def generate_conformers(mol, config):
    charge = config['molecule']['charge']
    base, ext = os.path.splitext(mol)
    if not ext in ['.xyz', '.pdb']:
        raise ValueError('Only pdb and xyz files are supported. Current '
                         'extension {}.'.format(ext))
    # Sanity check
    _check_charge(mol, charge)

    os.makedirs('crest', exist_ok=True)
    # Call xtb to do an initial round of optimisation
    with open('pyresp2_0_xtb_opt.log', 'a+') as f:
        cmd = '{xtb} {mol} --opt --alpb h2o --parallel {n_proc} ' \
              '--namespace crest/opt --chrg {c}'.format(
            xtb=config['bin_path']['xtb'], mol=mol,
            n_proc=config['global']['n_proc'], c=config['molecule']['charge'])
        f.write('''The initial xtb optimisation input is:
==============================================================================
{}
==============================================================================
The initial xtb optimisation output is:
'''.format(cmd))
    with open('pyresp2_0_xtb_opt.log', 'a+') as f:
        subprocess.call(cmd, shell=True, stdout=f, stderr=f)

    shutil.copy('crest/opt.xtbopt.pdb', 'crest/crest_in.pdb')
    with open('pyresp2_1_crest.log', 'a+') as f:
        cmd = '{crest} {mol} -aplb h2o -T {n_proc} -chrg {c} ' \
              '{additional}'.format(
            crest=config['bin_path']['crest'], mol='crest/crest_in.pdb',
            n_proc=config['global']['n_proc'], c=config['molecule']['charge'],
            additional=config['crest']['additional'])
        f.write('''The crest conformational search input is:
    ==============================================================================
    {}
    ==============================================================================
    The crest conformational search output is:
    '''.format(cmd))

    with open('pyresp2_1_crest.log', 'a+') as f:
        subprocess.call(cmd, shell=True, stdout=f, stderr=f)

    shutil.copy('crest_conformers.xyz', 'crest')

