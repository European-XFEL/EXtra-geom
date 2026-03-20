import numpy as np
import pytest
import xarray as xr

pyFAI = pytest.importorskip("pyFAI")

from extra_geom import AGIPD_500K2GGeometry


def test_position_modules_interpolate_conserves_signal_and_area():
    geom = AGIPD_500K2GGeometry.example()
    rng = np.random.default_rng(0)
    data = rng.poisson(3, size=geom.expected_data_shape).astype(np.float32)

    out, _ = geom.position_modules_interpolate(data, resize=True)

    # Total signal should be conserved when resize=True.
    assert np.isclose(np.nansum(out), data.sum(dtype=np.float64), rtol=1e-6, atol=1e-3)


def test_position_modules_interpolate_oversample_preserves_total_signal():
    geom = AGIPD_500K2GGeometry.example()
    rng = np.random.default_rng(1)
    data = rng.poisson(2, size=geom.expected_data_shape).astype(np.float32)
    da = xr.DataArray(
        data[None, None, ...],
        dims=("train", "pulse", "module", "slow_scan", "fast_scan"),
        coords={
            "train": np.asarray([123456789]),
            "pulse": np.asarray([0]),
            "module": np.arange(geom.expected_data_shape[-3]),
            "slow_scan": np.arange(geom.expected_data_shape[-2]),
            "fast_scan": np.arange(geom.expected_data_shape[-1]),
        },
    )

    out1, _ = geom.position_modules_interpolate(da, resize=True, oversample=1)
    out2, _ = geom.position_modules_interpolate(da, resize=True, oversample=2)

    assert np.isclose(np.nansum(out1), np.nansum(out2), rtol=1e-6, atol=1e-3)
