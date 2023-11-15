import ctypes
import os
import numpy as np

source_dir = os.path.join(os.path.dirname(__file__), "../")
lib_name = os.path.join(source_dir, "lib/pyElemSymPoly.so")
lib = ctypes.CDLL(lib_name)

lib = ctypes.CDLL("pyElemSymPoly.so")


compute_elem_sym_poly = lib.compute_elem_sym_poly
compute_elem_sym_poly.argtypes = [ctypes.POINTER(ctypes.c_double), ctypes.c_int, ctypes.c_int]
compute_elem_sym_poly.restype = ctypes.POINTER(ctypes.c_double)

free_memory = lib.free_poly
free_memory.argtypes = [ctypes.POINTER(ctypes.c_double)]


def is_power_of_2(n):
    return (n & (n-1) == 0) and n != 0


def next_power_of_2(num):
    if is_power_of_2(num):
        return num

    next_pow = 1
    while next_pow < num:
        next_pow <<= 1

    return next_pow


def elem_sym_poly(input_array, precision = 128):
    if isinstance(input_array, np.ndarray):
        array_len_ = input_array.shape[0]
        array_len = next_power_of_2(array_len_)

        pad_width = array_len - array_len_
        new_array = np.pad(input_array, (0, pad_width), mode='constant')

        array = np.ascontiguousarray(new_array, dtype=np.float64)
        array_ptr = array.ctypes.data_as(ctypes.POINTER(ctypes.c_double))

    elif isinstance(input_array, list):
        array_len_ = len(input_array)
        array_len = next_power_of_2(array_len_)

        pad_width = array_len - array_len_
        new_array = np.array(input_array, dtype=np.float64)
        new_array = np.pad(new_array, (0, pad_width), mode='constant')

        array = np.ascontiguousarray(new_array, dtype=np.float64)
        array_ptr = array.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    else:
        raise ValueError("Unsupported input array type. Must be a NumPy array or a Python list.")

    if not isinstance(precision, int) or (precision < 32 or precision > 513):
        raise ValueError("Unsupported precision type. Must be an Int between 32 and 512")

    polys = compute_elem_sym_poly(array_ptr, array_len, precision)

    num_polys = array_len_ + 1
    polys_list = [polys[i] for i in range(num_polys)]

    free_memory(polys)

    return polys_list


print(elem_sym_poly([1, 2, 3, 4, 5, 6, 7, 8, 9, 10]))
