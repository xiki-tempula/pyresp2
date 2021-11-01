import os.path
import shutil

import numpy as np
import pytest
import configparser
from pkg_resources import resource_filename

from pyresp2.charge import gen_charge, read_charge, read_all_charges

@pytest.fixture(scope='class')
def config():
    config = configparser.ConfigParser(inline_comment_prefixes="#")
    config.read(resource_filename(__name__, '../default.ini'))
    return config

def test_gen_charge(config):
    os.makedirs('charge')
    for i in range(5):
        shutil.copy(resource_filename(__name__, 'test_data/opt_{}.xyz'.format(i)),
                    'charge/charge_{}.xyz'.format(i))
    gen_charge(config)
    for i in range(5):
        assert os.path.isfile('charge/charge_{}.gjf'.format(i))
    shutil.rmtree('charge')

def test_read_charge():
    mol = read_charge(resource_filename(__name__, 'test_data/opt.mol2'),
                      resource_filename(__name__, 'test_data/RESP.log'))
    assert np.isclose(mol.atoms.charges[0], -0.2506055214, 0.01)

def test_read_charge():
    os.makedirs('charge', exist_ok=True)
    read_all_charges(resource_filename(__name__, 'test_data/opt.mol2'),
                     HF=resource_filename(__name__, 'test_data/RESP.log'),
                     Vac=resource_filename(__name__, 'test_data/RESP.log'),
                     Sol=resource_filename(__name__, 'test_data/RESP.log'))
    assert os.path.isfile('charge/RESP2.mol2')
    assert os.path.isfile('charge/HF.mol2')
