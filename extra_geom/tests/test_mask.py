import numpy as np

from extra_geom import mask

def test_delta_method():
    A = np.zeros((3, 4), dtype=bool)
    A[1:2, 1:3] = True
    B = np.copy(A)
    assert mask.delta_method(A) == [((1, 3), (1, 2))], "Simple rectangle test."
    assert np.array_equal(A, B), (
        "Make sure function does not change the original matrix.")

    C = np.zeros((5, 5), dtype=bool)
    C[0:3, 1:3] = C[1:3, 0:4] = C[2:4, 2:5] = C[4:5, 1:4] = True
    exp_res_C = [((1, 3), (0, 3)), ((0, 1), (1, 3)), ((3, 4), (1, 5)),
                 ((4, 5), (2, 4)), ((2, 3), (3, 5)), ((1, 2), (4, 5))]
    assert mask.delta_method(C) == exp_res_C, "Complex shape test."

    D = np.zeros((5, 5), dtype=bool)
    D[0:1, 0:1] = D[2:3, 2:3] = True
    exp_res_D = [((0, 1), (0, 1)), ((2, 3), (2, 3))]
    assert mask.delta_method(D) == exp_res_D, "Test on separate pixels."


def simple_mask():
    return mask.MaskRegions((2, 4, 5), [
        mask.RegionRect(None, 0, 1, 0, 2),
        mask.RegionRect(0, 2, 4, 0, 1),
        mask.RegionRect(1, 3, 4, 0, 5),
    ])

def test_roundtrip_array():
    m1 = simple_mask()
    arr = m1.to_mask_array()
    np.testing.assert_array_equal(arr, np.array([
        [[1, 1, 0, 0, 0],
         [0, 0, 0, 0, 0],
         [1, 0, 0, 0, 0],
         [1, 0, 0, 0, 0]],
        [[1, 1, 0, 0, 0],
         [0, 0, 0, 0, 0],
         [0, 0, 0, 0, 0],
         [1, 1, 1, 1, 1]],
    ], dtype=np.bool_))
    m2 = mask.MaskRegions.from_mask_array(arr)
    assert m2 == m1
