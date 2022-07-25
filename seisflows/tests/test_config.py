"""
Test the SeisFlows configuration script, which configures the compute
system and the working environment required for SF to run properly
"""
import pytest

from seisflows import config



def test_custom_import():
    """
    Test that importing based on internal modules works for various inputs
    :return:
    """
    with pytest.raises(SystemExit):
        config.custom_import()
    with pytest.raises(SystemExit):
        config.custom_import(name="NOT A VALID NAME")

    module = config.custom_import(name="optimize", module="LBFGS")
    assert(module.__name__ == "LBFGS")
    assert(module.__module__ == "seisflows.optimize.LBFGS")

    # Check one more to be safe
    module = config.custom_import(name="preprocess", module="default")
    assert(module.__name__ == "Default")
    assert(module.__module__ == "seisflows.preprocess.default")


