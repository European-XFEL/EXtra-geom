Release Notes
=============

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
