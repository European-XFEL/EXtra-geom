import numpy as np
from scipy.ndimage import gaussian_filter1d

import pyximport
pyximport.install(language_level=3)
from .cython_utils import condat_shift

# For more information, see Condat's paper
SQRT3 = np.sqrt(3)
R_HEX = np.sqrt(2 / SQRT3) * np.array([[1, 1 / 2], \
                                       [0, SQRT3 /2 ]])
R_CART = np.linalg.inv(R_HEX)

def pad4fft(tile):
    w0 = tile.shape[0] // 2
    w1 = tile.shape[1] // 2

    # Can we assume that a tile length is always a power of 2?
    width0 = (w0, w0) if tile.shape[0] % 2 == 0 else (w0, w0 + 1)
    width1 = (w1, w1) if tile.shape[1] % 2 == 0 else (w1, w1 + 1)

    return np.pad(tile, (width0, width1), mode="reflect")

def shear_new(tile, delay, axis, out=None):
    if out is None:
        out = np.empty_like(tile)

    cx = round(tile.shape[0] / 2)
    cy = round(tile.shape[1] / 2)

    # Vertical shear
    if axis == 0:
        for k in range(2 * cy):
            shift = delay * (k - cy)
            full_pixel_shift = round(shift)

            # Full pixel shifts
            out[:, k] = np.roll(tile[:, k], full_pixel_shift)

            # Sub-pixel shifts
            condat_shift(out[:, k], shift - full_pixel_shift)
            out[:, k] = gaussian_filter1d(out[:, k], 0.6)
    # Horizontal shear
    elif axis == 1:
        for k in range(2 * cx):
            shift = delay * (k - cx)
            full_pixel_shift = round(shift)

            # Full pixel shifts
            out[k, :] = np.roll(tile[k, :], full_pixel_shift)

            # Sub-pixel shifts
            condat_shift(out[k, :], shift - full_pixel_shift)
            out[k, :] = gaussian_filter1d(out[k, :], 0.6)

    return out

def hex2cart(tile):
    padded_tile = pad4fft(tile)
    cx = int(padded_tile.shape[0] / 2)
    cy = int(padded_tile.shape[1] / 2)

    # First we skew the image to a hexagonal shape
    height = tile.shape[0]
    y0 = height % 2 # Can we assume the tile lengths are always powers of 2?
    for y in range(-height // 2, height // 2 + y0 + 1):
        roll_amount = -int(np.floor((y + y0) / 2))
        offset = y - y0 + cx

        # Uncomment this line and comment out the next one to go back
        # to the original code.
        padded_tile[offset, :] = np.roll(padded_tile[offset, :], roll_amount)

    # Create shear coefficients
    delay1 = SQRT3 - np.sqrt(6 / SQRT3)
    delay2 = np.sqrt(SQRT3 / 6)
    delay3 = 2 - np.sqrt(6 / SQRT3)

    # Apply shears
    sheared = shear_new(padded_tile, delay3, axis=0)
    sheared = shear_new(sheared, delay2, axis=1)
    sheared = shear_new(sheared, delay1, axis=0)

    # Extract the rectangular box
    height = round(np.dot(R_CART, [tile.shape[0], 0])[0])
    width = round(np.dot(R_CART, [0, tile.shape[1]])[1])
    y0 = height % 2

    return sheared[-height // 2 + cx:height // 2 + cx + y0,
                   -width // 2 + cy:width // 2 + cy + 1]
