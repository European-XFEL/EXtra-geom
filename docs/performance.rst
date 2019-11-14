Performance notes
=================

These are some notes on how to load and process data efficiently.

Reduce before assembling
------------------------

Assembling detector images (see :doc:`geometry`) is relatively slow.
If your analysis involves a reduction step like summing or averaging over
a number of images, try to do this on the data from separate modules before
assembling them into images.

This also applies more generally: if a step in your processing makes the data
smaller, you want to do that step as near the start as possible.
