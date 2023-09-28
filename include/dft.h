#ifndef RASCH_DFT_H
#define RASCH_DFT_H

#include <roots.h>
#include <vector.h>
#include <stdint.h>
#include <stdio.h>

void print_bits(unsigned int num);


struct complex_vector * dft(struct vector * coeff_vec);

static int reverse_bits_odd( int x , int max_num_bits, int mask_base_index);

static int reverse_bits_even(int x , int max_num_bits, int base_index);

struct complex_vector * reverse_copy(struct vector * vec);

static const unsigned int masks[]= {
    0b01, 0b10, //2
    0b101, 0b010, 0b0011, 0b1100, //6
    0b0101, 0b1010, 0b0011,0b1100, // 10
    0b10101, 0b01010, 0b110011, 0b001100, 0b00001111, 0b11110000, //16
    0b010101, 0b101010, 0b110011, 0b001100, 0b00001111, 0b11110000, //22
    0b1010101, 0b0101010, 0b00110011, 0b11001100, 0b00001111, 0b11110000, //28
    0b01010101, 0b10101010, 0b00110011, 0b11001100, 0b00001111, 0b11110000, // 34
    0b101010101, 0b010101010, 0b1100110011, 0b0011001100, 0b111100001111, 0b000011110000, 0b0000000011111111, 0b1111111100000000, // 42
    0b0101010101, 0b1010101010, 0b1100110011, 0b0011001100, 0b111100001111, 0b000011110000, 0b0000000011111111, 0b1111111100000000, // 50
    0b10101010101, 0b01010101010, 0b001100110011, 0b110011001100, 0b111100001111, 0b000011110000, 0b0000000011111111, 0b1111111100000000, //58
    0b010101010101, 0b101010101010, 0b001100110011, 0b110011001100, 0b111100001111, 0b000011110000, 0b0000000011111111, 0b1111111100000000 //66
};

static const unsigned int mask_base_index[] = {
    0,
    0,
    0,
    2,
    6,
    10,
    16,
    22,
    28,
    34,
    42,
    50,
    58
};



#endif // RASCH_DFT_H