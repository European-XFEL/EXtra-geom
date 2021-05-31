import numpy as np


def all_tiles(geom):
    for mod in geom.modules:
        yield from mod


def tile_positions(geom):
    return np.stack([t.corner_pos for t in all_tiles(geom)])


def tile_vectors(geom):
    return np.stack([(t.ss_vec, t.fs_vec) for t in all_tiles(geom)])


def assert_geom_close(g1, g2):
    assert len(g1.modules) == len(g2.modules)
    assert {len(m) for m in g1.modules} == {len(m) for m in g2.modules}

    np.testing.assert_allclose(tile_positions(g1), tile_positions(g2))
    np.testing.assert_allclose(tile_vectors(g1), tile_vectors(g2))
