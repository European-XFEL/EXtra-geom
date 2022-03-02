Frequently Asked Questions
==========================

Why is my image "flipped" when I display it?
--------------------------------------------

At EuXFEL, the convention is to show detector images looking at the front of the detector, meaning the beam coming from behind you. You can get the other view with EXtra-geom's plotting methods by passing ``frontview=False``.

If you use matplotlib's ``imshow(array)`` to display an image, it may be upside down. You can correct this by passing ``imshow(array, origin='lower')``. If you also want to flip the displayed image horizontally, pass ``arr[:, ::-1]`` or ``np.fliplr(arr)``.

What are "slow-scan" and "fast-scan"?
-------------------------------------

We talk about "slow-scan" and "fast-scan" rather than just "x" and "y" axes, because the axes of the data array in HDF5 files do not always represent the same physical directions. For example, working along the fast-scan axis might be increasing y in one AGIPD module, and decreasing y in another.

Fast-scan is the inner dimension of the data, so a row of pixels along the fast-scan dimension are stored together.
One way to see how the dimensions are arranged is to use ``geom.inspect()`` in EXtra-geom. The dashed 'first row' lines in the diagram are always in the fast-scan direction.

EXtra-geom provides a method, :meth:`~.AGIPD_1MGeometry.data_coords_to_positions`, which converts (module_no, slow_scan, fast_scan) coordinates to physical positions.
