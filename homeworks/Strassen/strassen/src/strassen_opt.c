#include"matrix.h"

void sub_matrix_blocks(float **C, float const *const *const A,
                       float const *const *const B, const size_t C_f_row,
                       const size_t C_f_col, const size_t A_f_row,
                       const size_t A_f_col, const size_t B_f_row,
                       const size_t B_f_col, const size_t n);

void sum_matrix_blocks(float **C, float const *const *const A,
                       float const *const *const B, const size_t C_f_row,
                       const size_t C_f_col, const size_t A_f_row,
                       const size_t A_f_col, const size_t B_f_row,
                       const size_t B_f_col, const size_t n);

void naive_aux(float **C, float const *const *const A,
               float const *const *const B, const size_t C_f_row,
               const size_t C_f_col, const size_t A_f_row,
               const size_t A_f_col, const size_t B_f_row,
               const size_t B_f_col, const size_t n);

void strassen_opt_aux(float **C, float const *const *const A,
                  float const *const *const B,  float **M,
                  const size_t C_f_row, const size_t C_f_col,
                  const size_t A_f_row, const size_t A_f_col,
                  const size_t B_f_row, const size_t B_f_col,
                  const size_t M_f_row, const size_t M_f_col,
                  const size_t n )
{
  if(n<=128)
  {
    naive_aux(C, A, B, C_f_row, C_f_col, A_f_row, A_f_col, B_f_row, B_f_col, n);
    return;
  }

  size_t n2 = n/2;

  //C11 = P5 + P4 - P2 + P6
  //P5 = (A11 + A22)X(B11+B22) --> M12
  //A11 + A22 in C22, B11 + B22 in C21
  sum_matrix_blocks(C, A, A, C_f_row+n2, C_f_col+n2, A_f_row, A_f_col, A_f_row+n2, A_f_col+n2, n2);
  sum_matrix_blocks(C, B, B, C_f_row+n2, C_f_col, B_f_row, B_f_col, B_f_row+n2, B_f_col+n2, n2);
  strassen_opt_aux(M, (const float *const *const)C, (const float *const *const)C, M, M_f_row, M_f_col+n2, C_f_row+n2, C_f_col+n2, C_f_row+n2, C_f_col, M_f_row+n2, M_f_col+n2, n2);
  //P6 = (A12 - A22)X(B21+B22) --> unused after this step --> in M11
  //A12-A22 in C22, B21+B22 in C12
  sub_matrix_blocks(C, A, A, C_f_row+n2, C_f_col+n2, A_f_row, A_f_col+n2, A_f_row+n2, A_f_col+n2, n2);
  sum_matrix_blocks(C, B, B, C_f_row, C_f_col+n2, B_f_row+n2, B_f_col, B_f_row+n2, B_f_col+n2, n2);
  strassen_opt_aux(M, (const float *const *const)C, (const float *const *const)C, M, M_f_row, M_f_col, C_f_row+n2, C_f_col+n2, C_f_row, C_f_col+n2, M_f_row+n2, M_f_col+n2, n2);
  //P4 = A22 X (B21 - B11) --> M21
  //B21 - B11 in C22
  sub_matrix_blocks(C, B, B, C_f_row+n2, C_f_col+n2, B_f_row+n2, B_f_col, B_f_row, B_f_col, n2);
  strassen_opt_aux(M, (const float *const *const)A, (const float *const *const)C, M, M_f_row+n2, M_f_col, A_f_row+n2, A_f_col+n2, C_f_row+n2, C_f_col+n2, M_f_row+n2, M_f_col+n2, n2);
  //P2 = (A11 + A12)XB22 --> M21
  //A11 + A12 in C22
  sum_matrix_blocks(C, A, A, C_f_row+n2, C_f_col+n2, A_f_row, A_f_col, A_f_row, A_f_col+n2, n2);
  strassen_opt_aux(C, (const float *const *const)C, (const float *const *const)B, M, C_f_row, C_f_col+n2, C_f_row+n2, C_f_col+n2, B_f_row+n2, B_f_col+n2, M_f_row+n2, M_f_col+n2, n2);
  //final computation of C11
  sum_matrix_blocks(C, (const float *const *const)M, (const float *const *const)M, C_f_row, C_f_col, M_f_row, M_f_col+n2, M_f_row+n2, M_f_col, n2);
  sub_matrix_blocks(C, (const float *const *const)C, (const float *const *const)C, C_f_row, C_f_col, C_f_row, C_f_col, C_f_row, C_f_col+n2, n2);
  sum_matrix_blocks(C, (const float *const *const)C, (const float *const *const)M, C_f_row, C_f_col, C_f_row, C_f_col, M_f_row, M_f_col, n2);

  //C12 = P1 + P2
  //P1 = A11 X (B12-B22) --> in C22
  //P2 --> already in M21 --> no more used after this
  //B12-B22 in C21
  sub_matrix_blocks(C, B, B, C_f_row+n2, C_f_col, B_f_row, B_f_col+n2, B_f_row+n2, B_f_col+n2, n2);
  strassen_opt_aux(C, (const float *const *const)A, (const float *const *const)C, M, C_f_row+n2, C_f_col+n2, A_f_row, A_f_col, C_f_row+n2, C_f_col, M_f_row+n2, M_f_col+n2, n2);
  sum_matrix_blocks(C, (const float *const *const)C, (const float *const *const)C, C_f_row, C_f_col+n2, C_f_row+n2, C_f_col+n2, C_f_row, C_f_col+n2, n2);

  //C21 = P3 + P4
  //P3 = (A21+A22)XB11 --> in M11
  //A21 + A22 in C21
  sum_matrix_blocks(C, A, A, C_f_row+n2, C_f_col, A_f_row+n2, A_f_col, A_f_row+n2, A_f_col+n2, n2);
  strassen_opt_aux(M, (const float *const *const)C, (const float *const *const)B, M, M_f_row, M_f_col, C_f_row+n2, C_f_col, B_f_row, B_f_col, M_f_row+n2, M_f_col+n2, n2);
  //P4 --> already in M21 --> no more used after this
  sum_matrix_blocks(C, (const float *const *const)M, (const float *const *const)M, C_f_row+n2, C_f_col, M_f_row, M_f_col, M_f_row+n2, M_f_col, n2);


  //C22 = P5 + P1 - P3 - P7
  //P1 --> already in C22
  //P5 --> already in M12
  //P3 --> already in M11
  //P5 + P1 in C22
  sum_matrix_blocks(C, (const float *const *const)C, (const float *const *const)M, C_f_row+n2, C_f_col+n2, C_f_row+n2, C_f_col+n2, M_f_row, M_f_col+n2, n2);
  //(P5+P1)-P3 in C22
  sub_matrix_blocks(C, (const float *const *const)C, (const float *const *const)M, C_f_row+n2, C_f_col+n2, C_f_row+n2, C_f_col+n2, M_f_row, M_f_col, n2);
  //P7 = (A11 - A21) X (B11 + B12)
  //A11-A21 in M11, B11 + B12 in M12
  sub_matrix_blocks(M, A, A, M_f_row, M_f_col, A_f_row, A_f_col, A_f_row+n2, A_f_col, n2);
  sum_matrix_blocks(M, B, B, M_f_row, M_f_col+n2, B_f_row, B_f_col, B_f_row, B_f_col+n2, n2);
  //P7 in M21 = M11 X M12
  strassen_opt_aux(M, (const float *const *const)M, (const float *const *const)M, M, M_f_row+n2, M_f_col, M_f_row, M_f_col, M_f_row, M_f_col+n2, M_f_row+n2, M_f_col+n2, n2);
  //C22 = C22 - M21
  sub_matrix_blocks(C, (const float *const *const)C, (const float *const *const)M, C_f_row+n2, C_f_col+n2, C_f_row+n2, C_f_col+n2, M_f_row+n2, M_f_col, n2);
}

void strassen_opt(float **C, float const *const *const A, float const *const *const B, size_t n)
{
 //in this memory-optimezed version of the Strassen Algorithm, an auxiliary nxn matrix M is used
 //to store partial results (the "old" matrices P and S), as well as some unused blocks of the matrix C
 //the block (M_f_row+n2, M_f_col+n2) is passed recursively
 float **M = allocate_matrix(n, n);
 strassen_opt_aux(C, A, B, M, 0, 0, 0, 0, 0, 0, 0, 0, n);
 deallocate_matrix(M, n);
}
