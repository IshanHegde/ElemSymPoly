#include <dft.h>
#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <panopticon.h>

void print_bits(unsigned int num) {
    int num_bits = sizeof(num) * 8;  // Number of bits in the variable

    for (int i = num_bits - 1; i >= 0; i--) {
        int bit = (num >> i) & 1;
        printf("%d", bit);
    }

    printf("\n");
}




static int reverse_bits_odd(int x , int max_num_bits, int base_index)
{
    int aux_shift_num = 1;

    int m = base_index;


    
    while (aux_shift_num < max_num_bits){
        x = ((x & masks[m]) << aux_shift_num ) | ((x & masks[m+1]) >> aux_shift_num);
        m +=2;
        aux_shift_num *= 2;
    }



    return x >> max_num_bits -2;
}

static int reverse_bits_even(int x , int max_num_bits, int base_index)
{
    int aux_shift_num = 1;

    int m = base_index;


    
    while (aux_shift_num < max_num_bits){
        x = ((x & masks[m]) << aux_shift_num ) | ((x & masks[m+1]) >> aux_shift_num);
        m +=2;
        aux_shift_num *= 2;
    }



    return x ;
}


struct complex_vector * reverse_copy(struct complex_vector * vec){

    assert(vec->size % 2 == 0);

    int n = vec->size;

    int max_num_bits = ceil(log2(n-1));
    int mask_base_index_ = mask_base_index[max_num_bits];


    struct complex_vector * A = alloc_complex_vector(n);


    int i, n_half, index;
    double val;

    n_half = n / 2;


    if (max_num_bits % 2 == 1){
        for(i = 0; i < n; i++){
            index = reverse_bits_odd(i, max_num_bits, mask_base_index_);
            val = get_complex_vector_element(vec, index);
            
            set_complex_vector_element(A, i, val);
        }
    }
    else{
        for(i = 0; i < n; i++){
            index = reverse_bits_even(i, max_num_bits, mask_base_index_);
            val = get_complex_vector_element(vec, index);
            
            set_complex_vector_element(A, i, val);
        }
    }
    

    return A;

}

void dft(struct complex_vector * coeff_vec, struct complex_vector * A){

    
    int n = coeff_vec->size;

    // check if the size of the input is a power of 2
    // if ((n & (n - 1)) != 0){
    //     int new_size = 1 << (int) ceil(log2(n));
    //     coeff_vec = resize_vector(coeff_vec, new_size);
    //     out = resize_complex_vector(out, new_size);
    //     n = coeff_vec->size;
    // }

    

    //struct roots_of_unity * w = init_alloc_roots(n);

    //struct complex_vector * A = reverse_copy(coeff_vec);


    int k,s,m,j,i;

    

    double complex t,u, v, w_m, w;

    int max_iter = log2(n);

    // [ A_0, A_1, ..... A_N]
    // [ W_n^W_]
    // [ y_0, y_1, .... y_N]
    i =0;

    for (s = 1; s < max_iter+1; s++){

        m = pow(2,s);
        
        //printf("M val: %d \n",k);
        //w_m = init_alloc_roots(m);
        
        w_m = cexp(I*2*M_PI / m);
        
        w =1;
        //printf("W_m_k val: %f %f \n ",creal(w_m_k),cimag(w_m_k));

        
        for (j = 0; j < m/2; j++){

            
            WATCH("inner")
            for (k = j; k < n;k+=m){
                //WATCH("inner")
                //printf("K val: %d \n",k);
                //w_m_k = cpow(w->roots->data[k],n/m);
                
                //printf("val %d: \n", 2 * k / m);
                u = get_complex_vector_element(A,k);
                v = get_complex_vector_element(A,k+m/2);

                //printf("A val: %f %f \n ",creal(u),cimag(u));
                t = w * v;
                /*
                if (k==3){
                    printf("HERE %lf + %lfi \n", creal(w),cimag(w));

                    printf("vla %d ,,,, %d \n", k+m/2, m);

                }
                */
                set_complex_vector_element(A,k, u + t);

                set_complex_vector_element(A,k+m/2, u - t);
                
                i++;
                //STOP_WATCH("inner")
                
            }
            
            
            w = w_m *w;
            STOP_WATCH("inner")

        }
    
        //free_roots_of_unity(w_m);
    }


    // for (int i = 0;i < n; i++){
    //     set_complex_vector_element(out,i,get_complex_vector_element(A,i));
    // }

}


