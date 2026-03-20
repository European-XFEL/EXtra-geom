import numpy as np
import pytest
import xarray as xr

pyFAI = pytest.importorskip("pyFAI")

from extra_geom import AGIPD_500K2GGeometry, LPD_1MGeometry, DSSC_1MGeometry, AGIPD_1MGeometry


def sample_data(geom, as_xarray=False):
    rng = np.random.default_rng(0)
    data = rng.poisson(2, size=(2, 2, *geom.expected_data_shape)).astype(np.float32)

    if as_xarray:
        data = xr.DataArray(
            data,
            dims=("train", "pulse", "module", "slow_scan", "fast_scan"),
            coords={
                "train": np.asarray([123456789, 123456790]),
                "pulse": np.asarray([0, 1]),
                "module": np.arange(geom.expected_data_shape[-3]),
                "slow_scan": np.arange(geom.expected_data_shape[-2]),
                "fast_scan": np.arange(geom.expected_data_shape[-1]),
            },
        )
    return data


def test_conserves_signal_and_area():
    geom = AGIPD_500K2GGeometry.example()
    data = sample_data(geom)

    out, _ = geom.position_modules_interpolate(data, resize=True)

    # Total signal should be conserved when resize=True.
    assert np.isclose(np.nansum(out), data.sum(dtype=np.float64), rtol=1e-6, atol=1e-3)


def test_oversample_preserves_total_signal():
    geom = AGIPD_500K2GGeometry.example()
    data = sample_data(geom, as_xarray=True)

    with pytest.raises(ValueError, match='positive integer'):
        geom.position_modules_interpolate(data, resize=True, oversample=0.5)

    out1, _ = geom.position_modules_interpolate(data, resize=True, oversample=1)
    out2, _ = geom.position_modules_interpolate(data, resize=True, oversample=2)

    assert np.isclose(np.nansum(out1), np.nansum(out2), rtol=1e-6, atol=1e-3)


def test_output_pixel_size():
    geom = LPD_1MGeometry.example()
    data = sample_data(geom)

    out1, _ = geom.position_modules_interpolate(data)
    out2, _ = geom.position_modules_interpolate(data, output_pixel_size=geom.pixel_size * 2)
    out3, _ = geom.position_modules_interpolate(data, output_pixel_size=(geom.pixel_size, geom.pixel_size))

    assert out3.shape == out1.shape
    assert out2.shape[:-2] == out1.shape[:-2]
    assert out2.shape[-1] * 2 == out1.shape[-1]
    assert out2.shape[-2] * 2 == out1.shape[-2]


def test_xarray():
    geom = LPD_1MGeometry.example()
    data = sample_data(geom, as_xarray=True)

    geom.position_modules_interpolate(data)

    # bad module dim name
    bad_mod_dim = data.rename({'module': 'mods'})
    with pytest.raises(ValueError, match='module'):
        geom.position_modules_interpolate(bad_mod_dim)

    # bad module coordinates
    bad_mod_coord = data.assign_coords(module=data.module + 1)
    with pytest.raises(ValueError, match='module number labels should be in the range 0-15'):
        geom.position_modules_interpolate(bad_mod_coord)

    # bad data dimensions
    coords = data.coords.copy()
    coords['slow_scan'] = np.arange(512)
    coords['fast_scan'] = np.arange(128)
    bad_data_shape = xr.DataArray(
        data.data.reshape(*data.shape[:-2], 512, 128),
        dims=data.dims,
        coords=coords,
    )
    with pytest.raises(ValueError, match='Wrong pixel dimensions for detector modules'):
        geom.position_modules_interpolate(bad_data_shape)

    # missing modules
    missing_mods = data.isel(module=np.s_[:4])
    geom.position_modules_interpolate(missing_mods)


def test_hex_pixels():
    geom = DSSC_1MGeometry.example()
    data = sample_data(geom)

    with pytest.raises(NotImplementedError, match='this geometry has 6 corners per pixel.'):
        geom.position_modules_interpolate(data)


def plot_smoke():
    geom = AGIPD_1MGeometry.example()
    data = sample_data(geom)
    geom.plot_data(data[0, 0], interpolate=True)
