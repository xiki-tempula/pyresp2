import os.path
import shutil

import pytest
import configparser
from pkg_resources import resource_filename

from pyresp2.sp_eng import gen_sp, read_sp

@pytest.fixture(scope='class')
def config():
    config = configparser.ConfigParser(inline_comment_prefixes="#")
    config.read(resource_filename(__name__,
                                  '../default.ini'))
    return config

def test_gen_sp(config):
    os.symlink(resource_filename(__name__, 'test_data'), 'energy')
    gen_sp(config)
    for i in range(5):
        assert os.path.isfile('energy/opt_{}.inp'.format(i))
    os.unlink('energy')

def test_read_sp(config):
    os.symlink(resource_filename(__name__, 'test_data'), 'energy')
    read_sp(config)
    for i in range(5):
        assert os.path.isfile('charge/charge_{}.xyz'.format(i))

    assert os.path.isfile('charge/conflist.txt')
    assert os.path.isfile('charge/conflist_HF.txt')
    assert os.path.isfile('charge/conflist_Sol.txt')
    assert os.path.isfile('charge/conflist_Vac.txt')
    os.unlink('energy')
