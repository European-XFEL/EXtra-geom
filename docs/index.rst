European XFEL Python data tools
===============================

**EXtra_geom** is a Python library to describe the physical layout of
multi-module detectors at `European XFEL <https://www.xfel.eu/>`_, and to
assemble complete detector images.

Installation
------------

(TODO: update this)

karabo_data is available on our Anaconda installation on the Maxwell cluster::

    module load exfel exfel_anaconda3

You can also install it `from PyPI <https://pypi.org/project/karabo-data/>`__
to use in other environments with Python 3.5 or later::

    pip install karabo_data

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

