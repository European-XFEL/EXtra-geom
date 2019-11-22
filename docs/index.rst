European XFEL Python data tools
===============================

**EXtra-geom** is a Python library to describe the physical layout of
multi-module detectors at `European XFEL <https://www.xfel.eu/>`_, and to
assemble complete detector images.

Installation
------------

EXtra-geom is available on our Anaconda installation on the Maxwell cluster::

    module load exfel exfel_anaconda3

You can also install it `from PyPI <https://pypi.org/project/EXtra-geom/>`__
to use in other environments with Python 3.5 or later::

    pip install extra_geom

If you get a permissions error, add the ``--user`` flag to that command.

Documentation contents
----------------------

.. toctree::
   :maxdepth: 2

   geometry
   performance

.. toctree::
   :caption: Examples

   apply_geometry
   examine_geometry
   agipd_geometry
   dssc_geometry
   masks

.. toctree::
   :caption: Development
   :maxdepth: 1

   changelog

.. seealso::

   `Data Analysis at European XFEL
   <https://in.xfel.eu/readthedocs/docs/data-analysis-user-documentation/en/latest/>`_

Indices and tables
==================

* :ref:`genindex`
* :ref:`search`

