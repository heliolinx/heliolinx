#
# This file only exists to configure/build the C++ extension code.
# Everything else is done in pyproject.toml
#

from setuptools import setup, Extension
import pybind11
import numpy
from pybind11.setup_helpers import Pybind11Extension, build_ext

# Define the extension module
heliolinx = Pybind11Extension(
    'heliolinx.heliolinx',
    sources=['src/heliolinx.cpp', 'src/solarsyst_dyn_geo01.cpp'],
    extra_link_args=['-lgomp'],
    include_dirs=[
        pybind11.get_include(),  # Pybind11 include directory
        numpy.get_include()      # Numpy include directory
    ],
)

# Setup function
setup(ext_modules=[heliolinx])
