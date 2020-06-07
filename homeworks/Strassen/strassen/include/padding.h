#ifndef __PADDING__

void naive_rectangular(float **C, float const *const *const A, float const *const *const B,
                       const size_t A_r, const size_t A_c, const size_t B_c);

void strassen_padding(float** C, float const *const *const A,
                      const size_t A_r, const size_t A_c, //dimension of A
                      float const *const *const B,
                      const size_t B_r, const size_t B_c //dimensions of B
                      );

#endif //__PADDING__
