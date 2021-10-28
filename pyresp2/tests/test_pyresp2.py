"""
Unit and regression test for the pyresp2 package.
"""

# Import package, test suite, and other packages as needed
import sys

import pytest

import pyresp2


def test_pyresp2_imported():
    """Sample test, will always pass so long as import statement worked."""
    assert "pyresp2" in sys.modules
