Frequently Asked Questions
==========================

1. Why is my image "flipped" when I display it?

At EuXFEL, the convention is to show detector images looking at the front of the detector, meaning the beam coming from behind you. You can get the other view with EXtra-geom's plotting methods by passing `frontview=False`.

If you use matplotlib's `imshow(array)` to display an image, it may be upside down. You can correct this by passing `imshow(array, origin='lower')`.

2. What are "slow_scan" and "fast_scan"?

They are directions, or axes, in the HDF5 files. They are not always the same for different detectors. For example:

- In LPD, "slow_scan" corresponds to the y-axis and "fast_scan" to the x-axis.
- In AGIPD, "slow_scan" is the y-axis, whereas "fast_scan" is the y-axis

One way of know them it's to use the `geom.inspect()` in EXtra-geom, where you'll see diagrams where the dashed 'first row' lines are always in the "fast_scan" direction.
