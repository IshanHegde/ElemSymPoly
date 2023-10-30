static void _elem_symm_poly_comp(poly_mul_state_t poly_mul_state, matrix_t poly_matrix, int stride, int n){

    if ( n == 2){

        double poly_matrix_0_1 = poly_matrix[0][1];
        
        poly_matrix[0][0] = 1;
        poly_matrix[0][1] = poly_matrix_0_1 + poly_matrix[stride][1]; 
        poly_matrix[0][2] = poly_matrix_0_1 * poly_matrix[stride][1];
        poly_matrix[0][3] = 0;


    } else if (n == 4){

        _elem_symm_poly_comp(poly_mul_state, poly_matrix, 1, 2);
        _elem_symm_poly_comp(poly_mul_state, poly_matrix + 2, 1 , 2);

        double poly_matrix_0_1 = poly_matrix[0][1];
        double poly_matrix_0_2 = poly_matrix[0][2];
        double poly_matrix_0_3 = poly_matrix[0][3];

        poly_matrix[0][0] = 1;
        poly_matrix[0][1] = poly_matrix_0_1 + poly_matrix[stride][1];
        poly_matrix[0][2] = poly_matrix_0_1 * poly_matrix[stride][1] + poly_matrix_0_2 + poly_matrix[stride][2];
        poly_matrix[0][3] = poly_matrix_0_1 * poly_matrix[stride][2] + poly_matrix_0_2 * poly_matrix[stride][1] + poly_matrix_0_3 + poly_matrix[stride][3];
        poly_matrix[0][4] = poly_matrix_0_1 * poly_matrix[stride][3] + poly_matrix_0_2 * poly_matrix[stride][2] + poly_matrix_0_3 * poly_matrix[stride][1];
        poly_matrix[0][5] = poly_matrix_0_2 * poly_matrix[stride][3] + poly_matrix_0_3 * poly_matrix[stride][2];
        poly_matrix[0][6] = poly_matrix_0_3 * poly_matrix[stride][3];
        poly_matrix[0][7] = 0;


    } else{
        _elem_symm_poly_comp(poly_mul_state, poly_matrix, stride/2, n/2);
        _elem_symm_poly_comp(poly_mul_state, poly_matrix + n/2, stride/2 , n/2);

        for (int i =0 ;i < n; i+=2*stride){

            update_polynomial_mul_state(poly_mul_state, poly_matrix[i], poly_matrix[i + stride], n  );
            array_t result = polynomial_multiply(poly_mul_state);

            memcpy(poly_matrix[i], result, 2*n * sizeof(double));
            
        }

    }   
}