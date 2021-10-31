import os
from glob import glob
from os.path import join
import numpy as np
import shutil

def _FE2occup(energy_list, T=298.15):
    energy_list = np.array(energy_list)
    energy_list = energy_list - np.min(energy_list)
    Qi = np.exp(-energy_list/(T*0.0019858995))
    return Qi/np.sum(Qi)

def gen_sp(config):
    file_list = glob(join('energy', '*.xyz'))
    for file in file_list:
        base, ext = os.path.splitext(file)
        with open('{}.inp'.format(base), 'w') as f:
            f.write('''! {method}
%maxcore  {mem}
%pal nprocs   {n_proc} end
* xyzfile {charge}   {multiplicity} {file}
'''.format(file=os.path.basename(file),
           method=config['sp']['method'],
           n_proc=config['sp']['n_proc'],
           mem=int(config.getfloat('sp', 'memory') /
                   config.getfloat('sp', 'n_proc') * 0.75),
           charge=config['molecule']['charge'],
           multiplicity=config['molecule']['multiplicity']))

def read_sp(config):
    with open('pyresp2_3_sp.log', 'w') as f:
        f.write('Start examining the single point energy.')
    energy_text = 'FINAL SINGLE POINT ENERGY'
    min_energy = 0
    file_list = glob(join('energy', '*' + config['sp']['out_extension']))
    with open('pyresp2_3_sp.log', 'a+') as f:
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
            with open('pyresp2_3_sp.log', 'a+') as f:
                f.write('Read file: {} Energy: {} Min: {}\n'.format(file,
                                                                    energy,
                                                                    min_energy))
            energy_list.append((file, energy))
        else:
            with open('pyresp2_3_sp.log', 'a+') as f:
                f.write('Error reading file: {}\n'.format(file))

    energy_list = sorted(energy_list, key=lambda x: x[1])
    file_list = [point[0] for point in energy_list]
    # hatree to kcal/mol conversion: 627.5
    energy_list = [point[1] * 627.5 for point in energy_list]
    occupancy_list = _FE2occup(energy_list)
    os.makedirs('charge', exist_ok=True)
    count = 0
    output = []
    for file, occupancy in zip(file_list, occupancy_list):
        if occupancy > config.getfloat('sp', 'occupancy_cap'):
            with open('pyresp2_3_sp.log', 'a+') as f:
                f.write('File {} Occupancy {:.2f} included.\n'.format(file,
                                                                      occupancy))
            base, ext = os.path.splitext(file)
            shutil.copy2('{}.xyz'.format(base),
                         'charge/charge_{}.xyz'.format(count))
            output.append(('charge_{}'.format(count), occupancy))
            count += 1

    with open('charge/conflist.txt', 'w') as f:
        for file, occupancy in output:
            f.write('{} {:.2f}\n'.format(file, occupancy))

    with open('charge/conflist_HF.txt', 'w') as f:
        for file, occupancy in output:
            f.write('HF_{}.chk {:.2f}\n'.format(file.split('_')[-1], occupancy))

    with open('charge/conflist_Vac.txt', 'w') as f:
        for file, occupancy in output:
            f.write('Vac_{}.chk {:.2f}\n'.format(file.split('_')[-1], occupancy))

    with open('charge/conflist_Sol.txt', 'w') as f:
        for file, occupancy in output:
            f.write('Sol_{}.chk {:.2f}\n'.format(file.split('_')[-1], occupancy))


