import os
from subprocess import call
import MDAnalysis
import numpy as np

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
    call('{xtb} {mol} --opt --aplb h2o --parallel {n_proc} --namespace crest/opt'.format(
        xtb=config['bin_path']['xtb'], mol=mol,
        n_proc=config['global']['n_proc']), shell=True)
