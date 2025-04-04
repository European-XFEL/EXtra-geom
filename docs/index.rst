EXtra-geom: X-ray detector geometry for European XFEL
=====================================================

**EXtra-geom** is a Python library to describe the physical layout of
multi-module detectors at `European XFEL <https://www.xfel.eu/>`_, and to
assemble complete detector images.

Installation
------------

EXtra-geom is available in our Python environment on the Maxwell cluster::

    module load exfel exfel-python

You can also install it `from PyPI <https://pypi.org/project/EXtra-geom/>`__
to use in other environments with Python 3::

    pip install extra_geom

If you get a permissions error, add the ``--user`` flag to that command.

Documentation contents
----------------------

.. toctree::
   :maxdepth: 2

   geometry
   faq
   performance

.. toctree::
   :caption: Create/Load Examples

   agipd_geometry
   lpd_geometry
   dssc_geometry
   jungfrau_geometry
   epix_geometry
   generic_geometry
   
.. toctree::
   :caption: Use Examples
   
   apply_geometry
   convert_coords
   masks
   pyfai
   motor_based_geometry

.. toctree::
   :caption: Development
   :maxdepth: 1

   changelog

.. seealso::

   `Data Analysis at European XFEL
   <https://rtd.xfel.eu/docs/data-analysis-user-documentation/en/latest/>`_

Indices and tables
==================

* :ref:`genindex`
* :ref:`search`

