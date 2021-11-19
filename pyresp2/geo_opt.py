import os
from glob import glob
from os.path import join
import shutil
import MDAnalysis as mda
from MDAnalysis.analysis.align import alignto

def gen_opt(multi_xyz, config):
    os.makedirs('opt', exist_ok=True)
    u = mda.Universe(multi_xyz)
    for i, ts in enumerate(u.trajectory[:500]):
        pos = []
        for atom in u.atoms:
            pos.append(
                '{} {:.3f} {:.3f} {:.3f}'.format(atom.type, *atom.position))

        with open('opt/opt_{}.inp'.format(i), 'w') as f:
            f.write('''! opt {method}
%maxcore  {mem}
%pal nprocs   {n_proc} end
%cpcm
smd true
SMDsolvent "water"
end
* xyz   {charge}   {multiplicity}
{pos}
*
'''.format(pos='\n'.join(pos),
           method=config['opt']['method'],
           n_proc=config['opt']['n_proc'],
           mem=int(config.getfloat('opt', 'memory') /
                   config.getfloat('opt', 'n_proc') * 0.75),
           charge=config['molecule']['charge'],
           multiplicity=config['molecule']['multiplicity']))

def read_opt(config):
    with open('pyresp2_2_opt.log', 'w') as f:
        f.write('Start examining the optimised geometry.')
    energy_text = 'FINAL SINGLE POINT ENERGY'
    min_energy = 0
    file_list = glob(join('opt', '*' + config['opt']['out_extension']))
    with open('pyresp2_2_opt.log', 'a+') as f:
        f.write('Files detected: \n{}'.format('\n'.join(file_list)))
    energy_list = []
    for file in file_list:
        with open(file, 'r') as f:
            txt = f.read()
        if energy_text in txt:
            energy = 0
            lines = txt.split('\n')
            for line in lines:
                if energy_text in line:
                    energy = float(line.split()[-1])
                    if energy < min_energy:
                        min_energy = energy
            with open('pyresp2_2_opt.log', 'a+') as f:
                f.write('Read file: {} Energy: {} Min: {}\n'.format(file,
                                                                    energy,
                                                                    min_energy))
            energy_list.append((file, energy))
        else:
            with open('pyresp2_2_opt.log', 'a+') as f:
                f.write('Error reading file: {}\n'.format(file))
    with open('pyresp2_2_opt.log', 'a+') as f:
        f.write('A total of {} structures and the energy read in.\n'.format(
            len(energy_list)))
        f.write('Mininum energy {} hatree with energy cap of {} '
                'kcal/mol\n'.format(min_energy, config['opt']['energy_cap']))
        f.write('Find the distinct conformer over the RMSD cap of {} '
                'A.\n'.format(config['opt']['rmsd_cap']))

    os.makedirs('energy', exist_ok=True)
    energy_list = sorted(energy_list, key=lambda x: x[1])
    universe_list = []
    count = 0
    for file, energy in energy_list:
        if energy < min_energy + config.getfloat('opt', 'energy_cap') / 627.5:
            # hatree to kcal/mol conversion
            with open('pyresp2_2_opt.log', 'a+') as f:
                f.write('Found conformer within the energy cap: {}\n'.format(
                    file))
            base, ext = os.path.splitext(file)
            u = mda.Universe('{}.xyz'.format(base))
            for ref in universe_list:
                old_rmsd, new_rmsd = alignto(u, ref, select='not name H*')
                if new_rmsd < config.getfloat('opt', 'rmsd_cap'):
                    break
            else:
                with open('pyresp2_2_opt.log', 'a+') as f:
                    f.write(
                        'Conformer larger than RMSD cap copy to energy/opt_{}.xyz\n'.format(count))
                universe_list.append(u)
                shutil.copy2('{}.xyz'.format(base),
                             'energy/sp_{}.xyz'.format(count))
                count += 1
        if count > config.getint('opt', 'max_conf'):
            with open('pyresp2_2_opt.log', 'a+') as f:
                f.write(
                    'Maxinum number of conformer ({}) reached. '
                    'aborting.\n'.format(
                        config.getint('opt', 'max_conf')))
            break
