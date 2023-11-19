# ElemSymPoly
ElemSymPoly is a high-performance *C99* library with Python interface that utilizes the GNU MPFR Library for
arbitrarily precise computation of Elementary Symmetric Polynomials making use of a custom FFT implimentation. 

**The Library is currently in development stage and not production ready.**

### Installation

Recommended installation is via `pip`:

```bash
    pip install <placeholder>
```

Alternatively, you can install the library from source:

#### Requirements to build from source

- CMake
- GNU MPFR Library
- GNU GMP Library
- Python C header
- C compiler

#### Build from source

On linux, macOS and WSL:
```bash
    git clone https://github.com/IshanHegde/ElemSymPoly.git
    cd ElemSymPoly && mkdir build && cd build
    cmake ..
    make -j8
```

### General Idea

The library is based on applying a divide and conquer approach to compute the elementary symmetric polynomials.

- First, the elementary symmetric polynomials of order `N` is expressed as a product of `N` order `1` polynomials. 

- Next, using a divide and conquer approach, the product of `N` order 1 polynomials is computed by recursively.

- After a certain arbitrary threshold (currently set to `8`), the polynomials of order `8` or above are computed using FFT.

- The FFT algorithm is a custom implementation of the classic recursive Cooley-Tukey FFT algorithm. 

This algorithim has a time complexity of `O( N log^2 N )` compared to the naive approach of `O( N^2 )` and also has
arbitrary precision support (currently up to `512` decimal places due to stack overflow concerns). 

The library also has a Python wrapper for ease of use, and only relies on `NumPy`, `GNU MPFR` which in turn relies on
`GNU GMP`, Python C headers, CMake, glibc and a C compiler (only tested with GCC). 