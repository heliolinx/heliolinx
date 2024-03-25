## Run this script with docker, from the root of the repository, as follows:
##
## docker run --rm -v "$PWD":/io quay.io/pypa/manylinux2014_x86_64 /io/tools/build_wheels.sh
##
## The output (the .whl file) will be in the dist/ subdirectory.
##

cd /io
rm -rf wheelhouse/ dist/ build/
/opt/python/cp311-cp311/bin/pip wheel . -w wheelhouse/
auditwheel repair wheelhouse/heliolinx-*.whl -w dist/
chown -R 1000:1000 build dist wheelhouse heliolinx.egg-info
rm -rf build wheelhouse heliolinx.egg-info
