"""Parity checks for the numpy dtypes shared between C++ and Python.

Every struct that crosses the pybind11 boundary is declared in three places:

  1. the C++ struct in  src/solarsyst_dyn_geo01.h        (the source of truth)
  2. PYBIND11_NUMPY_DTYPE(...) in  src/heliolinx.cpp
  3. an np.dtype in  python/heliolinx/solarsyst_dyn_geo.py

pybind11 exposes these arguments as py::array_t<T>, whose default ExtraFlags
include py::array::forcecast. numpy therefore casts a mismatched structured
array to the expected dtype *positionally*, without regard to field names and
without raising. A field that is present in one declaration and absent from
another shifts every later field by one position and silently corrupts the
data -- the failure mode reported in issue #2.

Nothing in the build catches that, so these tests compare the three
declarations field-for-field and in order.
"""
import os
import re

import numpy as np
import pytest

from heliolinx import heliolinx as hl
from heliolinx import solarsyst_dyn_geo as sdg

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
HEADER = os.path.join(REPO, "src", "solarsyst_dyn_geo01.h")
BINDINGS = os.path.join(REPO, "src", "heliolinx.cpp")

# Every type registered with PYBIND11_NUMPY_DTYPE in src/heliolinx.cpp.
REGISTERED_TYPES = [
    "hldet",
    "EarthState",
    "hlimage",
    "longpair",
    "tracklet",
    "hlradhyp",
    "hlclust",
    "point3d_index",
    "glint_trail",
]


def header_fields(name):
    """Data-member names of `class/struct name { ... };` in the C++ header."""
    with open(HEADER) as fp:
        src = fp.read()
    match = re.search(r'^(?:class|struct)\s+%s\s*\{(.*?)^\};' % re.escape(name),
                      src, re.M | re.S)
    assert match is not None, "no struct or class %r in %s" % (name, HEADER)
    body = re.sub(r'//[^\n]*', '', match.group(1))
    fields = []
    for line in body.split('\n'):
        line = line.strip()
        if not line or line in ("public:", "private:", "protected:"):
            continue
        if '(' in line:                      # constructor or method
            continue
        member = re.match(r'^((?:unsigned|signed|long|short|const)\s+)*'
                          r'[A-Za-z_][A-Za-z0-9_]*\s+'
                          r'([A-Za-z_][A-Za-z0-9_]*)\s*(?:\[[^\]]*\])?\s*;', line)
        if member:
            fields.append(member.group(2))
    assert fields, "parsed no data members out of %s" % name
    return fields


def macro_fields(name):
    """Field names listed in PYBIND11_NUMPY_DTYPE(name, ...)."""
    with open(BINDINGS) as fp:
        src = fp.read()
    match = re.search(r'PYBIND11_NUMPY_DTYPE\(\s*%s\s*,(.*?)\)\s*;' % re.escape(name),
                      src, re.S)
    assert match is not None, "no PYBIND11_NUMPY_DTYPE for %r in %s" % (name, BINDINGS)
    return [f.strip() for f in match.group(1).split(',') if f.strip()]


def fill_struct_fields(name):
    """Field names assigned in `void fill_struct(name & out, name const& in)`."""
    with open(BINDINGS) as fp:
        src = fp.read()
    match = re.search(r'void\s+fill_struct\(\s*%s\s*&\s*out\s*,.*?\)\s*\{(.*?)^\}'
                      % re.escape(name), src, re.M | re.S)
    assert match is not None, "no fill_struct overload for %r in %s" % (name, BINDINGS)
    body = re.sub(r'//[^\n]*', '', match.group(1))
    assigned = re.findall(r'out\.([A-Za-z_][A-Za-z0-9_]*)\s*=', body)
    assigned += re.findall(r'memcpy\(\s*out\.([A-Za-z_][A-Za-z0-9_]*)', body)
    return assigned


@pytest.mark.parametrize("name", REGISTERED_TYPES)
def test_macro_matches_header(name):
    """PYBIND11_NUMPY_DTYPE lists exactly the header's fields, in order."""
    assert macro_fields(name) == header_fields(name)


@pytest.mark.parametrize("name", REGISTERED_TYPES)
def test_registered_dtype_matches_header(name):
    """The dtype the extension module actually exposes matches the header."""
    dtype = getattr(hl, "create_" + name)(1).dtype
    assert list(dtype.names) == header_fields(name)


@pytest.mark.parametrize("name", REGISTERED_TYPES)
def test_python_dtype_matches_header(name):
    """The np.dtype in solarsyst_dyn_geo.py matches the header.

    Four registered types have no Python dtype at all; passing a hand-built
    array for one of them relies entirely on the caller getting the field
    order right. Recorded as a gap rather than asserted here.
    """
    dtype = getattr(sdg, name, None)
    if dtype is None:
        pytest.skip("solarsyst_dyn_geo.py defines no np.dtype for %s" % name)
    assert isinstance(dtype, np.dtype)
    assert list(dtype.names) == header_fields(name)


@pytest.mark.parametrize("name", REGISTERED_TYPES)
def test_fill_struct_copies_every_field(name):
    """fill_struct() assigns every data member the header declares.

    ndarray_to_vec default-initialises its destination, so a member that
    fill_struct forgets is left indeterminate rather than zeroed.
    """
    assigned = fill_struct_fields(name)
    missing = [f for f in header_fields(name) if f not in assigned]
    assert not missing, "fill_struct(%s&) never assigns %s" % (name, ", ".join(missing))


def test_load_earth_ephemerides_matches_EarthState(tmp_path):
    """The Earth loader emits the field names EarthState declares.

    load_earth_ephemerides() builds its dtype inline rather than from a
    module-level np.dtype, so the parity check above cannot see it, but its
    output goes straight to observer_vel()/observer_coords() as an
    py::array_t<EarthState>. Field names that differ are reinterpreted
    positionally and silently.
    """
    horizons = tmp_path / "earth.txt"
    horizons.write_text(
        "$$SOE\n"
        "2459000.500000000, A.D. 2020-May-31 00:00:00.0000,"
        "  4.383e-01, -9.230e-01,  1.1e-05,  1.545e-02,  7.372e-03, -3.5e-07,\n"
        "2459001.500000000, A.D. 2020-Jun-01 00:00:00.0000,"
        "  4.537e-01, -9.155e-01,  1.0e-05,  1.532e-02,  7.633e-03, -3.5e-07,\n"
        "$$EOE\n"
    )
    earth = sdg.load_earth_ephemerides(str(horizons))
    assert list(earth.dtype.names) == header_fields("EarthState")
