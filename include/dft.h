#ifndef RASCH_DFT_H
#define RASCH_DFT_H

#include <roots.h>
#include <vector.h>
#include <stdint.h>
#include <stdio.h>

void init_look_up_table(int N, float * restrict reals, float * restrict imags);

void init_look_up_inverse(int N, float * restrict reals, float * restrict imags);

void recursive_fft( float * restrict in_reals,  float *  restrict in_imags,  float * restrict out_reals,  float * restrict out_imags, float * restrict w_reals,float * restrict w_imags , int stride, int n);

void recursive_inverse_fft( float * restrict in_reals,  float *  restrict in_imags,  float * restrict out_reals,  float * restrict out_imags, float * restrict w_reals,float * restrict w_imags , int stride, int n);
#endif // RASCH_DFT_H