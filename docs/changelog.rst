Release Notes
=============

1.0
---

- Added support for AGIPD 'mini-half' detector (8 modules) - see
  :ref:`det-AGIPD-500K2G` (:ghpull:`26`).
- Added methods to write XFEL HDF5 geometry files and get quadrant positions
  from geometry objects (:ghpull:`24`).
- Fixed y-axis scale in metres for plotting DSSC data (:ghpull:`23`).
- Faster image assembly with less overhead (:ghpull:`16`).
- Allow parallel image assembly using a thread pool (:ghpull:`17`), which can
  speed up assembling several images to a single 3D array.

0.10
----

- Added support for pnCCD detector (:ghpull:`13`).

0.9
---

- Initial support for :ref:`det-JUNGFRAU` detectors (:ghpull:`6`).
- Fix :meth:`~.AGIPD_1MGeometry.compare` method to draw arrows the right size
  (:ghpull:`4`).
- New example showing how to construct masks: :doc:`masks` (:ghpull:`1`).
- Correct code in :meth:`.LPD_1MGeometry.from_h5_file_and_quad_positions`
  which was working only by numeric coincidence (:ghpull:`7`).

0.8
---

First separated version. No functional changes from karabo_data 0.7.

Earlier history
---------------

The code in EXtra-geom was previously released as *karabo_data*, up to version
0.7. See the `karabo_data release notes
<https://karabo-data.readthedocs.io/en/latest/changelog.html>`_ for changes
before the separation.
