import os.path
import shutil

import pytest
import configparser
from pkg_resources import resource_filename

from pyresp2.geo_opt import gen_opt, read_opt

@pytest.fixture(scope='class')
def config():
    config = configparser.ConfigParser(inline_comment_prefixes="#")
    config.read(resource_filename(__name__,
                                  '../default.ini'))
    return config

class TestGen_opt():
    @staticmethod
    @pytest.fixture(scope='class')
    def run(config):
        gen_opt(resource_filename(__name__,
                                  'test_data/crest_conformers.xyz'), config)

    def test_files(self, run):
        for i in range(21):
            assert os.path.isfile('opt/opt_{}.inp'.format(i))
        shutil.rmtree('opt')

def test_read_opt(config):
    os.symlink(resource_filename(__name__,'test_data'), 'opt')
    read_opt(config)
    for i in range(5):
        assert os.path.isfile('energy/opt_{}.xyz'.format(i))
    os.rmdir('opt')
    shutil.rmtree('energy')
