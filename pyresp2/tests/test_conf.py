import os.path

import pytest
from pkg_resources import resource_filename
import configparser

from pyresp2.conf import generate_conformers

class TestGenerate_conformers():
    @staticmethod
    @pytest.fixture(scope='class')
    def run():
        config = configparser.ConfigParser(inline_comment_prefixes="#")
        config.read(resource_filename(__name__,
                                              '../default.ini'))
        config['crest']['additional'] = '-gff'
        generate_conformers(resource_filename(__name__,
                                              'test_data/butane.pdb'),
                            config)

    def test_opt(self, run):
        # Test if the initial xtb optimisation works.
        assert os.path.isfile('crest/opt.xtbopt.pdb')

    def test_crest(self, run):
        # Test if the crest run works.
        assert os.path.isfile('crest/crest_conformers.xyz')