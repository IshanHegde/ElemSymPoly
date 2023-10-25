for (int k =0; k < n/2; k+=4){

    w = LOAD_4(&w_reals[aux_num][k],&w_imags[aux_num][k]);

    y_1_k = LOAD_4(&out_reals[k+n/2],&out_imags[k+n/2]);
    y_0_k = LOAD_4(&out_reals[k],&out_imags[k]);
    t = MUL_4(w,y_1_k);
    
    STORE_4(&out_reals[k],&out_imags[k],ADD_4(y_0_k,t));
    STORE_4(&out_reals[k+n/2],&out_imags[k+n/2],SUB_4(y_0_k,t));

}

=> 
n/2 / 8

for (int k =0; k < n/16; k+=4){

    w = LOAD_4(&w_reals[aux_num][k],&w_imags[aux_num][k]);

    y_1_k = LOAD_4(&out_reals[k+n/2],&out_imags[k+n/2]);
    y_0_k = LOAD_4(&out_reals[k],&out_imags[k]);
    t = MUL_4(w,y_1_k);
    
    STORE_4(&out_reals[k],&out_imags[k],ADD_4(y_0_k,t));
    STORE_4(&out_reals[k+n/2],&out_imags[k+n/2],SUB_4(y_0_k,t));

}


for (int k =n/16; k < n/8 ; k+=4){

    w = LOAD_4(&w_reals[aux_num][k],&w_imags[aux_num][k]);

    y_1_k = LOAD_4(&out_reals[k+n/2],&out_imags[k+n/2]);
    y_0_k = LOAD_4(&out_reals[k],&out_imags[k]);
    t = MUL_4(w,y_1_k);
    
    STORE_4(&out_reals[k],&out_imags[k],ADD_4(y_0_k,t));
    STORE_4(&out_reals[k+n/2],&out_imags[k+n/2],SUB_4(y_0_k,t));

}


for (int k =n/8; k < 3*n/16 ; k+=4){

    w = LOAD_4(&w_reals[aux_num][k],&w_imags[aux_num][k]);

    y_1_k = LOAD_4(&out_reals[k+n/2],&out_imags[k+n/2]);
    y_0_k = LOAD_4(&out_reals[k],&out_imags[k]);
    t = MUL_4(w,y_1_k);
    
    STORE_4(&out_reals[k],&out_imags[k],ADD_4(y_0_k,t));
    STORE_4(&out_reals[k+n/2],&out_imags[k+n/2],SUB_4(y_0_k,t));

}


for (int k =3*n/16; k < n/4 ; k+=4){

    w = LOAD_4(&w_reals[aux_num][k],&w_imags[aux_num][k]);

    y_1_k = LOAD_4(&out_reals[k+n/2],&out_imags[k+n/2]);
    y_0_k = LOAD_4(&out_reals[k],&out_imags[k]);
    t = MUL_4(w,y_1_k);
    
    STORE_4(&out_reals[k],&out_imags[k],ADD_4(y_0_k,t));
    STORE_4(&out_reals[k+n/2],&out_imags[k+n/2],SUB_4(y_0_k,t));

}