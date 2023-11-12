import ctypes
import os
import sys

sys.path.append(os.path.join(os.path.dirname(__file__), "../"))

mpfr_t = ctypes.c_void_p
array_t = ctypes.POINTER(mpfr_t)
matrix_t = ctypes.POINTER(array_t)


class PolyMulState(ctypes.Structure):
    _fields_ = [
        ("N", ctypes.c_int),
        ("A_reals", array_t),
        ("A_imags", array_t),
        ("B_reals", array_t),
        ("B_imags", array_t),
        ("A_out_reals", array_t),
        ("A_out_imags", array_t),
        ("B_out_reals", array_t),
        ("B_out_imags", array_t),
        ("C_out_reals", array_t),
        ("C_out_imags", array_t),
        ("temp_C_reals", array_t),
        ("temp_C_imags", array_t),
        ("w_reals", matrix_t),
        ("w_imags", matrix_t),
        ("w_reals_inverse", matrix_t),
        ("w_imags_inverse", matrix_t)
    ]


class ElementarySymmetricState(ctypes.Structure):
    _fields_ = [
        ("N", ctypes.c_int),
        ("poly_mul_state", PolyMulState),
        ("aux_polys", matrix_t),
        ("elements", ctypes.POINTER(ctypes.c_double))
    ]


if __name__ == "__main__":
    source_dir = os.path.join(os.path.dirname(__file__), "../")

    lib_name = os.path.join(source_dir, "lib/elem_sym_poly.so")

    lib = ctypes.CDLL(lib_name)

    compute_elem_sym_poly = lib.compute_elem_sym_poly

    compute_elem_sym_poly.argtypes = [ctypes.POINTER(ctypes.c_double), ctypes.c_int]

    compute_elem_sym_poly.restype = ctypes.POINTER(ctypes.c_double)

    elements = [1,1,1,1]*64

    polys = compute_elem_sym_poly((ctypes.c_double * len(elements))(*elements),len(elements))

    num_polys = len(elements)+1
    polys_list = [polys[i] for i in range(num_polys)]

    print(polys_list)