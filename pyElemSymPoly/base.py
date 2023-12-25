import os
import sys
import numpy as np
import importlib
sys.path.insert(1, os.path.join(os.path.dirname(__file__), "../"))
try:
    sys.path.append(os.path.join(os.path.dirname(__file__), "../lib"))
except:
    pass
import libpyElemSymPoly as lib

def is_power_of_2(n):
    return (n & (n - 1) == 0) and n != 0


def next_power_of_2(num):
    if is_power_of_2(num):
        return num

    next_pow = 1
    while next_pow < num:
        next_pow <<= 1

    return next_pow


def elem_sym_poly(input_array, precision=128):
    if isinstance(input_array, np.ndarray):
        array_len_ = input_array.shape[0]
        array_len = next_power_of_2(array_len_)

        pad_width = array_len - array_len_
        new_array = np.pad(input_array, (0, pad_width), mode='constant').tolist()

    elif isinstance(input_array, list):
        array_len_ = len(input_array)
        array_len = next_power_of_2(array_len_)

        pad_width = array_len - array_len_
        new_array = np.array(input_array, dtype=np.float64)
        new_array = np.pad(new_array, (0, pad_width), mode='constant').tolist()
    else:
        raise ValueError("Unsupported input array type. Must be a NumPy array or a Python list.")

    if array_len_ < 2:
        raise ValueError(f"Unsupported input array size. Must be at least 2 instead of {array_len_}")

    if not isinstance(precision, int) or (precision < 32 or precision > 513):
        raise ValueError("Unsupported precision type. Must be an Int between 32 and 512")

    polys_list = lib.compute_elem_sym_poly(new_array, array_len_, precision)

    return polys_list

