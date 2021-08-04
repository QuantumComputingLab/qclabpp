# QCLAB++ [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5160682.svg)](https://doi.org/10.5281/zenodo.5160682)

<p align="center"><img src="doc/doxygen/images/logo200x200.png?raw=true" /></p>

QCLAB++ is an object-oriented, fully templated C++ package for creating and
representing quantum circuits. QCLAB++ can be used for rapid prototyping and
testing of quantum algorithms, and allows for fast algorithm development and
discovery. QCLAB++ provides I/O through openQASM making it compatible with
quantum hardware.


## How to run?

The QCLAB++ package uses the CMake build system (CMake version ≥ 3.16).
The recommended way of building QCLAB++ is as follows:

1. Install

        git clone https://github.com/QuantumComputingLab/qclabpp.git

2. CMake

        cd qclabpp
        mkdir build
        cd build
        cmake ..
        make -j8

3. Run tests

        ./test/qclab_tests

4. Generate documentation

        doxygen doxygen.dox


## Developers - Lawrence Berkeley National Laboratory
- [Roel Van Beeumen](http://www.roelvanbeeumen.be/) - rvanbeeumen@lbl.gov
- [Daan Camps](http://campsd.github.io/) - dcamps@lbl.gov


## Funding
The QCLAB++ project is supported by the Laboratory Directed Research and
Development Program of Lawrence Berkeley National Laboratory under U.S.
Department of Energy Contract No. DE-AC02-05CH11231.


## About
QCLAB++ Copyright (c) 2021, The Regents of the University of California,
through Lawrence Berkeley National Laboratory (subject to receipt of
any required approvals from the U.S. Dept. of Energy). All rights reserved.

If you have questions about your rights to use or distribute this software,
please contact Berkeley Lab's Intellectual Property Office at
IPO@lbl.gov.

NOTICE.  This Software was developed under funding from the U.S. Department
of Energy and the U.S. Government consequently retains certain rights. As
such, the U.S. Government has been granted for itself and others acting on
its behalf a paid-up, nonexclusive, irrevocable, worldwide license in the
Software to reproduce, distribute copies to the public, prepare derivative
works, and perform publicly and display publicly, and to permit others to do so.
