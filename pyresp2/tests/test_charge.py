import os.path
import shutil

import pytest
import configparser
from pkg_resources import resource_filename

from pyresp2.charge import gen_charge

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
