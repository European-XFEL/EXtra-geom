"""Integration with pyFAI

This module should be imported inside methods, so EXtra-geom works without
pyFAI installed.
"""

from pyFAI.detectors import Detector


class AGIPD1M(Detector):
    IS_CONTIGUOUS = False
    MAX_SHAPE = (16*512, 128)
    aliases = ["AGIPD 1M"]

    def __init__(self, pixel1=200e-6, pixel2=200e-6, **kwargs):
        super().__init__(pixel1, pixel2, **kwargs)
