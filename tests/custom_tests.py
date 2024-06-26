import sys
import os
import time
sys.path.append(os.path.join(os.path.dirname(__file__), "../"))
import pyElemSymPoly as esp
import random
from itertools import combinations
from functools import reduce
from decimal import Decimal
from decimal import getcontext

getcontext().prec = 28
random.seed()

# ------------- [ Constants ] -------------

FLOAT_64_EPS = sys.float_info.epsilon

# ------------- [ Reference Elementary Symmetric Polynomial Implementation] -------------

def reference_esp(lst):

    ret_lst = []

    for i in range(1,len(lst)+1):

        combs = list(combinations(lst, i))

        value = Decimal(0)

        for comb in combs:

            product = reduce(lambda x, y: x * y, comb)
            value += Decimal(product)
        
        ret_lst.append(float(value))
    
    return ret_lst

# ------------- [ Input different lengths from 0 to 65536 ] -------------

input_len_null = []

input_len_dict = {}

for num in range(1,4):

    max_len = 2**num
    cur_lst = []
    for index in range(0,max_len):
        cur_lst.append(random.random())

    input_len_dict[max_len] = cur_lst

# ------------- [ Execution of Tests ] -------------

for key, value in input_len_dict.items():


    lib_time_start = time.time()
    esp_values = esp.elem_sym_poly(value, precision = 128)
    lib_time_end = time.time()
    ref_time_start = time.time()
    ref_values = reference_esp(value)
    ref_time_end = time.time()

    for esp_val, ref_val in zip(esp_values, ref_values):

        assert abs(esp_val - ref_val) < 2 * FLOAT_64_EPS
