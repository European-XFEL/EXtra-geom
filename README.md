[![Build Status](https://github.com/European-XFEL/EXtra-geom/workflows/Tests/badge.svg)](https://github.com/European-XFEL/EXtra-geom/actions?query=workflow%3ATests)
[![codecov](https://codecov.io/gh/European-XFEL/EXtra-geom/branch/master/graph/badge.svg)](https://codecov.io/gh/European-XFEL/EXtra-geom)

Python tools to work with EuXFEL detector geometry and assemble detector images.

[Documentation](https://extra-geom.readthedocs.io/en/latest/)

Installing
==========

*EXtra-geom* is available on our Anaconda installation on the Maxwell cluster:

    module load exfel exfel_anaconda3

You can also install it [from PyPI](https://pypi.org/project/EXtra-geom/)
to use in other environments with Python 3.5 or later:

    pip install extra_geom

If you get a permissions error, add the `--user` flag to that command.


Contributing
===========

Tests
-----

Tests can be run as follows:

    python3 -m pytest -v
