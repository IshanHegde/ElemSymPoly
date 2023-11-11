#include <fft.h>
#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <vecmath.h>
#include <common.h>
#include <panopticon.h> 

#define data_t mpfr_t
#define array_t data_t * restrict
#define matrix_t data_t ** restrict

#define PRECISION_FFT 86

void init_look_up_table(int N, matrix_t reals, matrix_t imags){

	int outer_array_size = (int)log2(N);

	for (int i =0; i < outer_array_size;i++){

		int inner_array_size = (int) pow(2,i);
		int aux_size = inner_array_size * 2;
		reals[i] = ALLOC(sizeof(data_t) * inner_array_size);
		imags[i] = ALLOC(sizeof(data_t) * inner_array_size);

		data_t op, aux_0, aux_1, pi;
		mpfr_inits2(PRECISION_FFT, op, aux_0, aux_1, pi, (mpfr_ptr)NULL );

		mpfr_const_pi(pi, MPFR_RNDN);
		mpfr_div_ui(aux_0, pi, aux_size, MPFR_RNDN);
		mpfr_mul_ui(aux_1, aux_0, 2, MPFR_RNDN);

		for (int j =0;j < inner_array_size;j++){

			mpfr_init2(reals[i][j], PRECISION_FFT);
			mpfr_init2(imags[i][j], PRECISION_FFT);
			mpfr_mul_ui(op, aux_1, j, MPFR_RNDN);

			mpfr_cos(reals[i][j],op, MPFR_RNDN);
			mpfr_sin(imags[i][j],op, MPFR_RNDN);

		}

		mpfr_clears(op, aux_0, aux_1, pi, (mpfr_ptr)NULL);
	}
}

void init_look_up_inverse(int N, matrix_t reals, matrix_t imags){

	int outer_array_size = (int)log2(N);

	for (int i =0; i < outer_array_size;i++){

		int inner_array_size = (int) pow(2,i);
		int aux_size = inner_array_size * 2;
		reals[i] = ALLOC(sizeof(data_t) * inner_array_size);
		imags[i] = ALLOC(sizeof(data_t) * inner_array_size);

		data_t op, aux_0, aux_1, pi;
		mpfr_inits2(PRECISION_FFT, op, aux_0, aux_1, pi, (mpfr_ptr)NULL );

		mpfr_const_pi(pi, MPFR_RNDN);
		mpfr_div_ui(aux_0, pi, aux_size, MPFR_RNDN);
		mpfr_mul_si(aux_1, aux_0, -2, MPFR_RNDN);

		for (int j =0;j < inner_array_size;j++){

			mpfr_init2(reals[i][j], PRECISION_FFT);
			mpfr_init2(imags[i][j], PRECISION_FFT);
			mpfr_mul_ui(op, aux_1, j, MPFR_RNDN);

			mpfr_cos(reals[i][j],op, MPFR_RNDN);
			mpfr_sin(imags[i][j],op, MPFR_RNDN);

		}

		mpfr_clears(op, aux_0, aux_1, pi, (mpfr_ptr)NULL);
	}

}

void free_look_up_table(int N, matrix_t reals, matrix_t imags){

    int outer_array_size = (int)log2(N);

    for (int i =0; i < outer_array_size;i++){
        free(reals[i]);
        free(imags[i]);
    }

    free(reals);
    free(imags);
}

void recursive_rfft_half_zero_safe(array_t in_reals, array_t in_imags, array_t out_reals, array_t out_imags, matrix_t w_reals, matrix_t w_imags , int stride, int n){

	if (n == 4){

		/*
		  X_0 = x_0 + x_1
		  X_1 = x_0 + ix_1
		  X_2 = x_0 - x_1
		  X_3 = x_0 - ix_1

		  re(X_1) = re(x_0) - im(x_1)
		  im(X_1) = im(x_0) + re(x_1)
		 */

		mpfr_add(out_reals[0],in_reals[0],in_reals[stride], MPFR_RNDN);
		mpfr_add(out_imags[0],in_imags[0],in_imags[stride], MPFR_RNDN);

		mpfr_sub(out_reals[1],in_reals[0],in_imags[stride], MPFR_RNDN);
		mpfr_add(out_imags[1],in_imags[0],in_reals[stride], MPFR_RNDN);

		mpfr_sub(out_reals[2],in_reals[0],in_reals[stride], MPFR_RNDN);
		mpfr_sub(out_imags[2],in_imags[0],in_imags[stride], MPFR_RNDN);

		mpfr_add(out_reals[3],in_reals[0],in_imags[stride], MPFR_RNDN);
		mpfr_sub(out_imags[3],in_imags[0],in_reals[stride], MPFR_RNDN);
        /*
		puts("------");
		for (int i = 0; i < 4; i++){
			mpfr_printf("Input: %Rf +I %Rf\n", in_reals[i*stride], in_imags[i*stride]);
		}
		for (int i = 0; i < 4; i++){
			mpfr_printf("output: %Rf +I %Rf\n", out_reals[i], out_imags[i]);
		}
		puts("------");
        */
	}else{
		recursive_rfft_half_zero_safe(in_reals, in_imags, out_reals, out_imags, w_reals, w_imags, stride << 1, n >> 1);
		recursive_rfft_half_zero_safe(in_reals + stride, in_imags + stride, out_reals + n/2, out_imags + n/2, w_reals, w_imags, stride << 1, n >> 1);

		// twiddle factor outer array index
		int aux_num = log2(n)-1;
		mpfr_t t_real, t_imag, aux_0, aux_1, aux_2, aux_3;
		mpfr_inits2(PRECISION_FFT, t_real, t_imag, aux_0, aux_1, aux_2, aux_3, (mpfr_ptr)NULL);

		for (int k =0; k < n/2; k++){

			/*
			 * w = LOAD(&w_reals[aux_num][k],&w_imags[aux_num][k]);
            y_1_k = LOAD(&out_reals[k+n/2],&out_imags[k+n/2]);
            y_0_k = LOAD(&out_reals[k],&out_imags[k]);
            t = MUL(w,y_1_k);

            STORE(&out_reals[k],&out_imags[k],ADD(y_0_k,t));
            STORE(&out_reals[k+n/2],&out_imags[k+n/2],SUB(y_0_k,t));
			 */
			mpfr_set(aux_2,out_reals[k],MPFR_RNDN);
			mpfr_set(aux_3,out_imags[k],MPFR_RNDN);

			mpfr_mul(aux_0, w_reals[aux_num][k], out_reals[k + n/2], MPFR_RNDN);
			mpfr_mul(aux_1, w_imags[aux_num][k], out_imags[k + n/2], MPFR_RNDN);
			mpfr_sub(t_real, aux_0, aux_1, MPFR_RNDN);

			mpfr_mul(aux_0, w_reals[aux_num][k], out_imags[k + n/2], MPFR_RNDN);
			mpfr_mul(aux_1, w_imags[aux_num][k], out_reals[k + n/2], MPFR_RNDN);
			mpfr_add(t_imag, aux_0, aux_1, MPFR_RNDN);


			mpfr_add(out_reals[k], out_reals[k], t_real, MPFR_RNDN);
			mpfr_add(out_imags[k], out_imags[k], t_imag, MPFR_RNDN);

			mpfr_sub(out_reals[k + n/2], aux_2, t_real, MPFR_RNDN);
			mpfr_sub(out_imags[k + n/2], aux_3, t_imag, MPFR_RNDN);






		}
		mpfr_clears(t_real, t_imag, aux_0, aux_1, aux_2, aux_3, (mpfr_ptr)NULL);

	}

}

void recursive_inverse_fft_safe(array_t in_reals, array_t in_imags, array_t out_reals, array_t out_imags, matrix_t w_reals, matrix_t w_imags , int stride, int n){

	if (n == 4){


		/*
		out_reals[0] = o0_real + o2_real + o1_real + o3_real;
		out_imags[0] = o0_imag + o2_imag + o1_imag + o3_imag;

		out_reals[1] = o0_real - o2_real + o1_imag - o3_imag;
		out_imags[1] = o0_imag - o2_imag - o1_real + o3_real;

		out_reals[2] = o0_real + o2_real - o1_real - o3_real;
		out_imags[2] = o0_imag + o2_imag - o1_imag - o3_imag;

		out_reals[3] = o0_real - o2_real - o1_imag + o3_imag;
		out_imags[3] = o0_imag - o2_imag + o1_real - o3_real;
		*/
		mpfr_add(out_reals[0], in_reals[0], in_reals[2*stride], MPFR_RNDN);
		mpfr_add(out_reals[0], out_reals[0], in_reals[stride], MPFR_RNDN);
		mpfr_add(out_reals[0], out_reals[0], in_reals[3*stride], MPFR_RNDN);

		mpfr_add(out_imags[0], in_imags[0], in_imags[2*stride], MPFR_RNDN);
		mpfr_add(out_imags[0], out_imags[0], in_imags[stride], MPFR_RNDN);
		mpfr_add(out_imags[0], out_imags[0], in_imags[3*stride], MPFR_RNDN);

		mpfr_sub(out_reals[1], in_reals[0], in_reals[2*stride], MPFR_RNDN);
		mpfr_add(out_reals[1], out_reals[1], in_imags[stride], MPFR_RNDN);
		mpfr_sub(out_reals[1], out_reals[1], in_imags[3*stride], MPFR_RNDN);

		mpfr_sub(out_imags[1], in_imags[0], in_imags[2*stride], MPFR_RNDN);
		mpfr_sub(out_imags[1], out_imags[1], in_reals[stride], MPFR_RNDN);
		mpfr_add(out_imags[1], out_imags[1], in_reals[3*stride], MPFR_RNDN);

		mpfr_add(out_reals[2], in_reals[0], in_reals[2*stride], MPFR_RNDN);
		mpfr_sub(out_reals[2], out_reals[2], in_reals[stride], MPFR_RNDN);
		mpfr_sub(out_reals[2], out_reals[2], in_reals[3*stride], MPFR_RNDN);

		mpfr_add(out_imags[2], in_imags[0], in_imags[2*stride], MPFR_RNDN);
		mpfr_sub(out_imags[2], out_imags[2], in_imags[stride], MPFR_RNDN);
		mpfr_sub(out_imags[2], out_imags[2], in_imags[3*stride], MPFR_RNDN);

		mpfr_sub(out_reals[3], in_reals[0], in_reals[2*stride], MPFR_RNDN);
		mpfr_sub(out_reals[3], out_reals[3], in_imags[stride], MPFR_RNDN);
		mpfr_add(out_reals[3], out_reals[3], in_imags[3*stride], MPFR_RNDN);

		mpfr_sub(out_imags[3], in_imags[0], in_imags[2*stride], MPFR_RNDN);
		mpfr_add(out_imags[3], out_imags[3], in_reals[stride], MPFR_RNDN);
		mpfr_sub(out_imags[3], out_imags[3], in_reals[3*stride], MPFR_RNDN);


	}else{
		recursive_inverse_fft_safe(in_reals, in_imags, out_reals, out_imags, w_reals, w_imags, stride << 1, n >> 1);
		recursive_inverse_fft_safe(in_reals + stride, in_imags + stride, out_reals + n/2, out_imags + n/2, w_reals, w_imags, stride << 1, n >> 1);

		int aux_num = log2(n)-1;
		mpfr_t t_real, t_imag, aux_0, aux_1, aux_2, aux_3;
		mpfr_inits2(PRECISION_FFT, t_real, t_imag, aux_0, aux_1, aux_2, aux_3, (mpfr_ptr)NULL);

		for (int k =0; k < n/2; k++){

			/*
			 * w = LOAD(&w_reals[aux_num][k],&w_imags[aux_num][k]);
            y_1_k = LOAD(&out_reals[k+n/2],&out_imags[k+n/2]);
            y_0_k = LOAD(&out_reals[k],&out_imags[k]);
            t = MUL(w,y_1_k);

            STORE(&out_reals[k],&out_imags[k],ADD(y_0_k,t));
            STORE(&out_reals[k+n/2],&out_imags[k+n/2],SUB(y_0_k,t));
			 */
			mpfr_set(aux_2,out_reals[k],MPFR_RNDN);
			mpfr_set(aux_3,out_imags[k],MPFR_RNDN);

			mpfr_mul(aux_0, w_reals[aux_num][k], out_reals[k + n/2], MPFR_RNDN);
			mpfr_mul(aux_1, w_imags[aux_num][k], out_imags[k + n/2], MPFR_RNDN);
			mpfr_sub(t_real, aux_0, aux_1, MPFR_RNDN);

			mpfr_mul(aux_0, w_reals[aux_num][k], out_imags[k + n/2], MPFR_RNDN);
			mpfr_mul(aux_1, w_imags[aux_num][k], out_reals[k + n/2], MPFR_RNDN);
			mpfr_add(t_imag, aux_0, aux_1, MPFR_RNDN);


			mpfr_add(out_reals[k], aux_2, t_real, MPFR_RNDN);
			mpfr_add(out_imags[k], aux_3, t_imag, MPFR_RNDN);

			mpfr_sub(out_reals[k + n/2], aux_2, t_real, MPFR_RNDN);
			mpfr_sub(out_imags[k + n/2], aux_3, t_imag, MPFR_RNDN);


		}
		mpfr_clears(t_real, t_imag, aux_0, aux_1, aux_2, aux_3, (mpfr_ptr)NULL);

	}
}

#undef data_t
#undef array_t
#undef matrix_t