Release Notes
=============

1.5
---

- Add method to make geometry from ASIC pairs positions for ePix100 detector and method
  to normalize ePix data (:ghpull:`97`). See :ref:`det-EPIX`.
- Make scipy an optional dependency (:ghpull:`90`).
- Add method to make DSSC-1M geometry from only quadrant positions (:ghpull:`89`). See
  :doc:`dssc_geometry` (example) and :ref:`det-DSSC-1M` (reference).
- Fix method name in docstring (:ghpull:`84`).

1.4
---

- Added support for ePix100 & ePix10k detectors (:ghpull:`73`). See
  :doc:`epix_geometry` (example) and :ref:`det-EPIX` (reference).
- :meth:`.GenericGeometry.inspect` now labels modules and tiles if there
  are more than one (:ghpull:`74`).
- Allocating output arrays to assemble integer data should be faster
  (:ghpull:`78`).
- Use ``NotImplementedError`` to make it clear that creating
  :class:`.GenericGeometry` from a ``.geom`` file is not yet supported
  (:ghpull:`77`).
- Some code reorganisation (:ghpull:`75`, :ghpull:`76`).

1.3
---

- A new :class:`.GenericGeometry` class allows describing the layout of an unknown
  detector, with the user specifying details such as pixel size and number of
  modules (:ghpull:`72`). See :doc:`generic_geometry` for an introduction.
- Fix a small discrepancy in module positions with
  :meth:`.JUNGFRAUGeometry.from_module_positions` (:ghpull:`69`).

1.2.1
-----

- Fix assembling images from an ``extra_data`` StackView object (:ghpull:`67`).

1.2
---

- JUNGFRAU geometry can now be saved to and loaded from CrystFEL format
  ``.geom`` files - see the :doc:`JUNGFRAU example <jungfrau_geometry>`,
  :meth:`.JUNGFRAUGeometry.write_crystfel_geom` and
  :meth:`.JUNGFRAUGeometry.from_crystfel_geom` (:ghpull:`49`).
- Images can now be assembled from an ``xarray.DataArray`` with a dimension
  named 'module' labelled with module numbers counting from 0 (:ghpull:`62`).
- Fix how ``coffset`` information is handled when reading CrystFEL geometry
  files (:ghpull:`60`).
- :class:`.PNCCDGeometry`, added in 0.10, is now documented and tested
  (:ghpull:`45`).
- New :doc:`faq` document (:ghpull:`51`)
- Avoid importing ``scipy.ndimage`` unnecessarily (:ghpull:`50`).

1.1.1
-----

- Fix module orientation for AGIPD 500k detector (:ghpull:`41`).

1.1
---

- New :meth:`~.AGIPD_1MGeometry.position_modules_symmetric` method to assemble
  data with the detector centre at the midpoint of the output array
  (:ghpull:`31`).
- New :meth:`~.AGIPD_1MGeometry.offset` method to move part or all of a geometry
  in 2 or 3 dimensions (:ghpull:`27`).
- New function :func:`.agipd_asic_seams` to select or mask the double-width
  pixels where AGIPD tiles touch.
- Examples in documentation rearranged and improve (:ghpull:`32`, :ghpull:`36`).
- CI moved to Github Actions (:ghpull:`34`) and integrated with Dependabot to
  control new versions of dependencies (:ghpull:`35`).

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
