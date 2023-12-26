import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), "../"))
import pyElemSymPoly as esp
import random

random.seed()

# ------------- [ Constants ] -------------

FLOAT_64_EPS = sys.float_info.epsilon


# ------------- [ Reference Elementary Symmetric Polynomial Implementation] -------------

def reference_esp(lst):

    ret_lst = []

    pass

# ------------- [ Input different lengths from 0 to 65536 ] -------------

input_len_null = []

input_len_dict = {}

for num in range(0,4):

    max_len = 2**num
    cur_lst = []
    for index in range(0,max_len):
        cur_lst.append(random.random())

    input_len_dict[max_len] = cur_lst

# ------------- [ Execution of Tests ] -------------

for key, value in input_len_dict.items():

    print(key)

    print(esp.elem_sym_poly(value))