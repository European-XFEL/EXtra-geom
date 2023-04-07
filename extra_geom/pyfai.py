"""Integration with pyFAI

This module should be imported inside methods, so EXtra-geom works without
pyFAI installed.
"""

import numpy as np
from pyFAI.detectors import Detector


class AGIPD1M(Detector):
    IS_CONTIGUOUS = False
    MAX_SHAPE = (16*512, 128)
    aliases = ["AGIPD 1M"]

    def __init__(self, pixel1=200e-6, pixel2=200e-6, **kwargs):
        super().__init__(pixel1, pixel2, **kwargs)


class AGIPD500K(Detector):
    IS_CONTIGUOUS = False
    MAX_SHAPE = (8*512, 128)
    aliases = ["AGIPD 500K"]

    def __init__(self, pixel1=200e-6, pixel2=200e-6, **kwargs):
        super().__init__(pixel1, pixel2, **kwargs)


class DSSC1M(Detector):
    IS_CONTIGUOUS = False
    MAX_SHAPE = (16*128, 512)
    aliases = ["DSSC 1M"]

    def __init__(self, pixel1=236e-6 * 1.5/np.sqrt(3), pixel2=236e-6, **kwargs):
        super().__init__(pixel1, pixel2, **kwargs)


class LPD1M(Detector):
    IS_CONTIGUOUS = False
    MAX_SHAPE = (16*256, 256)
    aliases = ["LPD 1M"]

    def __init__(self, pixel1=5e-4, pixel2=5e-4, **kwargs):
        super().__init__(pixel1, pixel2, **kwargs)


class JUNGFRAU_EuXFEL(Detector):
    IS_CONTIGUOUS = False

    def __init__(self, pixel1=7.5e-5, pixel2=7.5e-5, n_modules=None, **kwargs):
        self.MAX_SHAPE = (n_modules*512, 1024)
        super().__init__(pixel1, pixel2, **kwargs)

def EPIX10K2M(Detector):
    IS_CONTIGUOUS = FALSE
    MAX_SHAPE = (16*352, 384)

    def __init__(self, pixel1=100e-6, pixel2=100e-6, **kwargs):
        super().__init__(pixel1, pixel2, **kwargs)
