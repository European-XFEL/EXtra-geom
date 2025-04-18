Release Notes
=============

1.15
----

2025-04-08

- Added support for :meth:`~.JUNGFRAUGeometry.data_coords_to_positions` on
  :class:`~.JUNGFRAUGeometry` (:ghpull:`346`).
- Fix for creating a PyFAI detector object for JUNGFRAU with
  :meth:`~.JUNGFRAUGeometry.to_pyfai_detector` (:ghpull:`334`).
- More informative error messages when input or output arrays for assembling
  images don't match the expected shape (:ghpull:`343`).
- The dependency on ``cfelpyutils`` has been replaced by ``cfel_fmt``, a fork
  to maintain the pieces we need (:ghpull:`344`).

1.14
----

2025-03-07

- Fix the direction of rotations in the :meth:`~.rotate` method to match what is
  described in the docs (right-hand convention). This is a breaking change for
  code which relied on the previous behaviour (:ghpull:`332`).

1.13
----

2025-01-28

- JUNGFRAU 4M (8 modules) geometry can now be updated based on the positions of
  the 2 motors moving the hemispheres (:ghpull:`303`), using
  :class:`~.JF4MMotors`.
- EXtra-geom now requires Python 3.9 or above (:ghpull:`309`).

1.12
----

2024-05-27

- :ref:`det-AGIPD-1M` detector geometry can now be updated using the positions
  from the 8 quadrant motors (:ghpull:`269`). A reference geometry is required
  as a starting point, and a new geometry can be created based on the changes in
  the motor positions. See :doc:`motor_based_geometry` for an example of how to
  use this.

1.11
----

2023-11-01

- Add an ``example()`` class method for each detector type to create a sample
  geometry, and helper class methods ``monolithic_geometry()`` and ``pair_geometry()``
  for Epix detectors  (:ghpull:`243`). The ``example()`` methods make it easier to look at the data when you
  have no specific information about its geometry, but it may be quite different
  from the real positions of modules in a given experiment. For some detector
  types, you need to pass the number of modules in your detector.
- Modules in :class:`.AGIPD_500K2GGeometry` are now labelled M0 - M7 by
  :meth:`~.AGIPD_500K2GGeometry.inspect` (:ghpull:`226`).

1.10
----

2023-04-21

- Support for LPD Mini detectors (:ghpull:`187`).
- Add :meth:`~.JUNGFRAUGeometry.to_pyfai_detector` for JUNGFRAU detectors
  (:ghpull:`197`)
- Fix :meth:`.AGIPD_500K2GGeometry.from_origin` with non-default units
  (:ghpull:`213`).
- Fix :meth:`~.JUNGFRAUGeometry.plot_data` with labelled data arrays
  (:ghpull:`198`).
- Fix DSSC's :meth:`~.DSSC_1MGeometry.position_modules_cartesian` and
  :meth:`~.DSSC_1MGeometry.plot_data_cartesian` to accept Xarray labelled array
  objects (:ghpull:`207`).
- Fix alias for LPD PyFAI detector class (:ghpull:`190`).


1.9
---

2022-10-28

- New methods to assemble DSSC data and convert the hexagonal pixels onto a
  square grid: :meth:`~.DSSC_1MGeometry.position_modules_cartesian` and
  :meth:`~.DSSC_1MGeometry.plot_data_cartesian` (:ghpull:`174`).

1.8
---

2022-10-18

- New method :meth:`.DSSC_1MGeometry.plot_data_hexes` to plot DSSC data, drawing
  a hexagon for each pixel. This is slower than regular plotting, but more
  accurately represents what the detector 'saw' (:ghpull:`167`).
- More useful labels, and an option to supply custom module labels, for JUNGFRAU
  geometry in :meth:`~.JUNGFRAUGeometry.inspect` (:ghpull:`177`).
- Fix assembling JUNGFRAU images from labelled array with module numbers
  starting from 1 (:ghpull:`169`).
- Fix a bug writing some geometry objects to CrystFEL format ``.geom``
  files (:ghpull:`163`)

1.7.1
-----

2022-05-19

- Fix the pattern of hexagonal pixels in the DSSC detector (:ghpull:`160`).
  Thanks to Loïc le Guyader for identifying and investigating this issue.

1.7
---

2022-03-02

- New method :meth:`~.DSSC_1MGeometry.to_pyfai_detector` for AGIPD, DSSC and LPD
  to make a PyFAI detector object (:ghpull:`139`). See :doc:`pyfai` for an example.
- New method :meth:`~.DSSC_1MGeometry.rotate` to rotate all or selected parts of
  the detector by given angles in 3D (:ghpull:`128`).
- Rename ``plot_data_fast`` to ``plot_data``, and ``position_modules_fast`` to
  ``position_modules`` (:ghpull:`143`). The old names remain as aliases.
- EXtra-geom now works with (and requires) cfelpyutils 2.x for reading
  CrystFEL format ``.geom`` files (:ghpull:`114`).

1.6
---

2021-09-22

- Store and read (with new :meth:`~.LPD_1MGeometry.from_h5_file` method)
  quadrant positions in EuXFEL HDF5 format geometry files (:ghpull:`92`).
- Read some metadata from CrystFEL format ``.geom`` files and use it as defaults
  when writing a new ``.geom`` file (:ghpull:`87`).
- Fix writing ``coffset`` (z coordinates) correctly in ``.geom`` files
  (:ghpull:`102`).
- Require cfelpyutils < 2 until we fix compatibility with the new version
  (:ghpull:`107`).

1.5
---

2021-08-30

- Add method to make geometry from ASIC pairs positions for ePix100 detector and method
  to normalize ePix data (:ghpull:`97`). See :ref:`det-EPIX`.
- Make scipy an optional dependency (:ghpull:`90`).
- Add method to make DSSC-1M geometry from only quadrant positions (:ghpull:`89`). See
  :doc:`dssc_geometry` (example) and :ref:`det-DSSC-1M` (reference).
- Fix method name in docstring (:ghpull:`84`).

1.4
---

2021-06-16

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

2021-05-20

- A new :class:`.GenericGeometry` class allows describing the layout of an unknown
  detector, with the user specifying details such as pixel size and number of
  modules (:ghpull:`72`). See :doc:`generic_geometry` for an introduction.
- Fix a small discrepancy in module positions with
  :meth:`.JUNGFRAUGeometry.from_module_positions` (:ghpull:`69`).

1.2.1
-----

2021-04-20

- Fix assembling images from an ``extra_data`` StackView object (:ghpull:`67`).

1.2
---

2021-04-16

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

2020-12-17

- Fix module orientation for AGIPD 500k detector (:ghpull:`41`).

1.1
---

2020-12-04

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

2020-10-01

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

2020-06-24

- Added support for pnCCD detector (:ghpull:`13`).

0.9
---

2020-03-24

- Initial support for :ref:`det-JUNGFRAU` detectors (:ghpull:`6`).
- Fix :meth:`~.AGIPD_1MGeometry.compare` method to draw arrows the right size
  (:ghpull:`4`).
- New example showing how to construct masks: :doc:`masks` (:ghpull:`1`).
- Correct code in :meth:`.LPD_1MGeometry.from_h5_file_and_quad_positions`
  which was working only by numeric coincidence (:ghpull:`7`).

0.8
---

2019-11-18

First separated version. No functional changes from karabo_data 0.7.

Earlier history
---------------

The code in EXtra-geom was previously released as *karabo_data*, up to version
0.7. See the `karabo_data release notes
<https://karabo-data.readthedocs.io/en/latest/changelog.html>`_ for changes
before the separation.
