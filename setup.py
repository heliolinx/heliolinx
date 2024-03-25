from setuptools import setup, Extension, find_packages
import pybind11
import numpy
from pybind11.setup_helpers import Pybind11Extension, build_ext

# Define the extension module
heliolinx = Pybind11Extension(
    'heliolinx.heliolinx',
    sources=['src/heliolinx.cpp', 'src/solarsyst_dyn_geo01.cpp'],
    include_dirs=[
        pybind11.get_include(),  # Pybind11 include directory
        numpy.get_include()      # Numpy include directory
    ],
)

# Setup function
setup(
#    name='heliolinx',
#    use_scm_version={
#        "root": ".",  # Directory of the root of the project
#        "relative_to": __file__,  # Location of the file relative to the root
#        "version_scheme": "guess-next-dev",  # Automatically guess the next development version
#        "local_scheme": "dirty-tag"  # Append '.dirty' if the working directory is dirty
#    },
#    setup_requires=['setuptools_scm'],  # Needed for setuptools_scm to work
#    author='Aren Heinze, Mario Juric, Siegfried Eggl',
#    author_email='aheinze@uw.edu, mjuric@astro.washington.edu, eggls@uw.edu',
#    description='heliolinx: the fast asteroid linking code.',
#    long_description='',
    ext_modules=[heliolinx],
#    cmdclass={"build_ext": build_ext},
#    packages=find_packages(),
#    zip_safe=False,
)
