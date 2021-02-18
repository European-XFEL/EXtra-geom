Frequently Asked Questions
==========================

Why is my image "flipped" when I display it?
--------------------------------------------

At EuXFEL, the convention is to show detector images looking at the front of the detector, meaning the beam coming from behind you. You can get the other view with EXtra-geom's plotting methods by passing ``frontview=False``.

If you use matplotlib's ``imshow(array)`` to display an image, it may be upside down. You can correct this by passing ``imshow(array, origin='lower')``.

What are "slow_scan" and "fast_scan"?
-------------------------------------

We talk about "slow_scan" and "fast_scan" rather than just "x" and "y" axes, because the axes of the data array in HDF5 files do not always represent the same physical directions. For example, working along the "fast_scan" axis might be increasing "y" in one AGIPD module, and decreasing "y" in another.

One way of know them it's to use the `geom.inspect()` in EXtra-geom, where you'll see diagrams where the dashed 'first row' lines are always in the "fast_scan" direction.

EXtra-geom provides a convenient method, `data_coords_to_position()<https://extra-geom.readthedocs.io/en/latest/geometry.html?highlight=data_coords_to#extra_geom.AGIPD_1MGeometry.data_coords_to_positions>`_, which converts (module_no, slow_scan, fast_scan) to physical positions.
