#ifndef RASCH_DFT_H
#define RASCH_DFT_H

#include <stdint.h>
#include <stdio.h>


void init_look_up_table_d(int N, double ** restrict reals, double ** restrict imags);

void init_look_up_inverse_d(int N, double ** restrict reals, double ** restrict imags);

void recursive_fft_d( double * restrict  in_reals,  double * restrict  in_imags,  double * restrict out_reals,  double * restrict out_imags, double ** restrict  w_reals, double ** restrict w_imags , int stride, int n);

void recursive_inverse_fft_d( double * restrict in_reals, double * restrict  in_imags, double * restrict out_reals, double * restrict out_imags, double ** restrict w_reals, double ** restrict w_imags , int stride, int n);

void recursive_rfft_half_zero_d(double * restrict in_reals, double * restrict  in_imags, double * restrict out_reals, double * restrict out_imags, double ** restrict w_reals, double ** restrict w_imags , int stride, int n);


#endif // RASCH_DFT_H