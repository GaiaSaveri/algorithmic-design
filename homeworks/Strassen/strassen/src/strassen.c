#include "matrix.h"
#include "strassen.h"

/*
 * this function performs the element-wise
 * sub_matrix_blockstraction of B from A and put the resulting
 * sub_matrix_blocks-matrix in C. The parameters *_f_row and *_f_col
 * represents the first row and the first column,
 * respectively, of the sub_matrix_blocks-matrix we want to deal with.
 */
void sub_matrix_blocks(float **C, float const *const *const A,
                       float const *const * const B,
                       const size_t C_f_row, const size_t C_f_col,
                       const size_t A_f_row, const size_t A_f_col,
                       const size_t B_f_row, const size_t B_f_col,
                       const size_t n)
{
    for (size_t y = 0; y < n; y++) {
      for (size_t x = 0; x < n; x++) {
          C[y+C_f_row][x+C_f_col] =
               A[y+A_f_row][x+A_f_col] - B[y+B_f_row][x+B_f_col];
      }
    }
}

/*
 * this function performs the element-wise
 * sum_matrix_blocks of A and B and put the resulting
 * sub_matrix_blocks-matrix in C. The parameters *_f_row and *_f_col
 * represents the first row and the first column,
 * respectively, of the sub_matrix_blocks-matrix we want to deal with.
 */
void sum_matrix_blocks(float **C, float const *const *const A,
                       float const *const * const B,
                       const size_t C_f_row, const size_t C_f_col,
                       const size_t A_f_row, const size_t A_f_col,
                       const size_t B_f_row, const size_t B_f_col,
                       const size_t n)
{
    for (size_t y = 0; y < n; y++) {
      for (size_t x = 0; x < n; x++) {
          C[y+C_f_row][x+C_f_col] =
               A[y+A_f_row][x+A_f_col] + B[y+B_f_row][x+B_f_col];
      }
    }
}

/*
 * this function implements the naive algorithm
 * for matrix multiplication between sub_matrix_blocks-matrixes.
 * The result is placed in the sub_matrix_blocks-matrix C.
 * The parameters *_f_row and *_f_col
 * represents the first row and the first column,
 * respectively, of the sub_matrix_blocks-matrix we want to deal with.
 */
void naive_aux(float **C, float const *const *const A,
               float const *const * const B,
               const size_t C_f_row, const size_t C_f_col,
               const size_t A_f_row, const size_t A_f_col,
               const size_t B_f_row, const size_t B_f_col,
               const size_t n)
{
  for (size_t y = 0; y < n; y++) {
    for (size_t x = 0; x < n; x++) {
      float value = 0.0;
      for (size_t z = 0; z < n; z++) {
        value += A[y + A_f_row][z + A_f_col]*B[z + B_f_row][x + B_f_col];
      }

      C[y + C_f_row][x + C_f_col] = value;
    }
  }
}

/*
 * This function implements the Strassen's algorithm
 * for matrix multiplication between sub_matrix_blocks-matrixes.
 * The result is placed in the sub_matrix_blocks-matrix C.
 * The parameters *_f_row and *_f_col
 * represents the first row and the first column,
 * respectively, of the sub_matrix_blocks-matrix we want to deal with.
 */
void strassen_aux(float **C, float const *const *const A,
                  float const *const * const B,
                  const size_t C_f_row, const size_t C_f_col,
                  const size_t A_f_row, const size_t A_f_col,
                  const size_t B_f_row, const size_t B_f_col,
                  const size_t n)
{
    if (n <= (1<<5)) {
        naive_aux(C, A, B,
                  C_f_row, C_f_col,
                  A_f_row, A_f_col,
                  B_f_row, B_f_col,
                  n);

        return;
    }

    size_t n2 = n/2; // This is the size of the blocks

    float ***S = (float ***)malloc(sizeof(float **) * 10);
    for (size_t i = 0; i < 10; i++) {
        S[i] = allocate_matrix(n2, n2);
    }

    float ***P = (float ***)malloc(sizeof(float **) * 7);
    for (size_t i = 0; i < 7; i++) {
        P[i] = allocate_matrix(n2, n2);
    }

    // S1 = B12 - B22
    sub_matrix_blocks(S[0], B, B,
                      0, 0,
                      B_f_row, B_f_col + n2,
                      B_f_row + n2, B_f_col + n2,
                      n2);

    // P1 = A11 x S1
    strassen_aux(P[0], A, (const float* const *const) S[0],
                 0, 0,
                 A_f_row, A_f_col,
                 0, 0,
                 n2);

    // S2 = A11 + A12
    sum_matrix_blocks(S[1], A, A,
                      0, 0,
                      A_f_row, A_f_col,
                      A_f_row, A_f_col + n2,
                      n2);


    // P2 = S2 x B22
    strassen_aux(P[1], (const float* const *const) S[1], B,
                 0, 0,
                 0, 0,
                 B_f_row + n2, B_f_col + n2,
                 n2);

    // S3 = A21 + A22
    sum_matrix_blocks(S[2], A, A,
                      0, 0,
                      A_f_row + n2, A_f_col,
                      A_f_row + n2, A_f_col + n2,
                      n2);

    // P3 = S3 x B11
    strassen_aux(P[2], (const float* const *const) S[2], B,
                 0, 0,
                 0, 0,
                 B_f_row, B_f_col,
                 n2);

    // S4 = B21 - B11
    sub_matrix_blocks(S[3], B, B,
                      0, 0,
                      B_f_row + n2, B_f_col,
                      B_f_row, B_f_col,
                      n2);

    // P4 = A22 x S4
    strassen_aux(P[3], A, (const float* const *const) S[3],
                 0, 0,
                 A_f_row + n2, A_f_col + n2,
                 0, 0,
                 n2);

    // S5 = A11 + A22
    sum_matrix_blocks(S[4], A, A,
                      0, 0,
                      A_f_row, A_f_col,
                      A_f_row + n2, A_f_col + n2,
                      n2);

    // S6 = B11 + B22
    sum_matrix_blocks(S[5], B, B,
                      0, 0,
                      B_f_row, B_f_col,
                      B_f_row + n2, B_f_col + n2,
                      n2);

    // P5 = S5 x S6
    strassen_aux(P[4], (const float* const *const) S[4],
                 (const float* const *const) S[5],
                 0, 0,
                 0, 0,
                 0, 0,
                 n2);

    // S7 = A12 - A22
    sub_matrix_blocks(S[6], A, A,
                      0, 0,
                      A_f_row, A_f_col + n2,
                      A_f_row + n2, A_f_col + n2,
                      n2);

    // S8 = B21 + B22
    sum_matrix_blocks(S[7], B, B,
                      0, 0,
                      B_f_row + n2, B_f_col,
                      B_f_row + n2, B_f_col + n2,
                      n2);

    // P6 = S7 x S8
    strassen_aux(P[5], (const float* const *const) S[6],
                 (const float* const *const) S[7],
                 0, 0,
                 0, 0,
                 0, 0,
                 n2);

    // S9 = A11 - A21
    sub_matrix_blocks(S[8], A, A,
                      0, 0,
                      A_f_row, A_f_col,
                      A_f_row + n2, A_f_col,
                      n2);

    // S10 = B11 + B12
    sum_matrix_blocks(S[9], B, B,
                      0, 0,
                      B_f_row, B_f_col,
                      B_f_row, B_f_col + n2,
                      n2);

    // P7 = S9 x S10
    strassen_aux(P[6], (const float* const *const) S[8],
                 (const float* const *const) S[9],
                 0, 0,
                 0, 0,
                 0, 0,
                 n2);

    // C11 = P5 + P4 - P2 + P6
    sum_matrix_blocks(C, (const float* const *const) P[4],
                      (const float* const *const) P[3],
                      C_f_row, C_f_col,
                      0, 0,
                      0, 0,
                      n2);
    sub_matrix_blocks(C, (const float* const *const) C,
                      (const float* const *const) P[1],
                      C_f_row, C_f_col,
                      C_f_row, C_f_col,
                      0, 0,
                      n2);
    sum_matrix_blocks(C, (const float* const *const) C,
                      (const float* const *const) P[5],
                      C_f_row, C_f_col,
                      C_f_row, C_f_col,
                      0, 0,
                      n2);

    // C12 = P1 + P2
    sum_matrix_blocks(C, (const float* const *const) P[0],
                      (const float* const *const) P[1],
                      C_f_row, C_f_col+n2,
                      0, 0,
                      0, 0,
                      n2);

    // C21 = P3 + P4
    sum_matrix_blocks(C, (const float* const *const) P[2],
                      (const float* const *const) P[3],
                      C_f_row+n2, C_f_col,
                      0, 0,
                      0, 0,
                      n2);

    // C22 = P5 + P1 - P3 - P7
    sum_matrix_blocks(C, (const float* const *const) P[4],
                      (const float* const *const) P[0],
                      C_f_row+n2, C_f_col+n2,
                      0, 0,
                      0, 0,
                      n2);
    sub_matrix_blocks(C, (const float* const *const) C,
                      (const float* const *const) P[2],
                      C_f_row+n2, C_f_col+n2,
                      C_f_row+n2, C_f_col+n2,
                      0, 0,
                      n2);
    sub_matrix_blocks(C, (const float* const *const) C,
                      (const float* const *const) P[6],
                      C_f_row+n2, C_f_col+n2,
                      C_f_row+n2, C_f_col+n2,
                      0, 0,
                      n2);

    for (size_t i = 0; i < 10; i++) {
      deallocate_matrix(S[i], n2);
    }
    free(S);

    for (size_t i = 0; i < 7; i++) {
      deallocate_matrix(P[i], n2);
    }
    free(P);
}


/*
 * this functions is exclusively meant to provide an
 * easy to use API
 */
void strassen_matrix_multiplication(float **C, float const *const *const A,
                                    float const *const *const B, size_t n)
{

  strassen_aux(C, A, B,0, 0, 0, 0,0, 0, n);

}

/**
 * This function implements the Strassen algorithm for matrix multiplication,
 * the difference with the one above is that here only six auxiliary matrices
 * are allocated.
 * Parameters are exactly the same used in the previous function. 
 */
void strassen_aux_optimal(float **C,float const *const * const A, float const *const * const B,
                          const size_t C_f_row, const size_t C_f_col,
                          const size_t A_f_row, const size_t A_f_col,
                          const size_t B_f_row, const size_t B_f_col,
                          const size_t n)
{
  if (n < (1<<6)) {
    naive_aux(C,A,B, C_f_row, C_f_col,
               A_f_row, A_f_col,
               B_f_row, B_f_col,
              n);

    return;
  }
  //just for practicality
  const size_t n2=n/2;

  const size_t C1X = C_f_row;
  const size_t C2X = C_f_row + n2;
  const size_t CX1 = C_f_col;
  const size_t CX2 = C_f_col + n2;

  const size_t A1X = A_f_row;
  const size_t A2X = A_f_row + n2;
  const size_t AX1 = A_f_col;
  const size_t AX2 = A_f_col + n2;

  const size_t B1X = B_f_row;
  const size_t B2X = B_f_row + n2;
  const size_t BX1 = B_f_col;
  const size_t BX2 = B_f_col + n2;

  //allocate just 2 S instead of 10
  float ***S=(float ***)malloc(sizeof(float **)*2);

  for (int i=0; i<2; i++) {
    S[i] = allocate_matrix(n2, n2);
  }
  //allocate just 4 P instead of 7
  float ***P=(float ***)malloc(sizeof(float **)*4);
  for (int i=0; i<4; i++) {
    P[i] = allocate_matrix(n2, n2);
  }

  //P5
  // S5 = A11 + A22
  sum_matrix_blocks(S[0],A,A,0,0,A1X,AX1,A2X,AX2,n2);

  // S6 = B11 + B22
  sum_matrix_blocks(S[1],B,B,0,0,B1X,BX1,B2X,BX2,n2);

  // P5 = S5 x S6
  strassen_aux_optimal(P[0],(const float* const *const)S[0],(const float* const *const)S[1],
                      0,0, 0, 0,0, 0,n2);

  //P4
  // S4 = B21 - B11
  sub_matrix_blocks(S[1],B,B,0,0,B2X,BX1,B1X,BX1,n2);

  // P4 = A22 x S4
  strassen_aux_optimal(P[1], A, (const float* const *const)S[1], 0, 0,A2X, AX2,0, 0,n2);

  //P2
  // S2 = A11 + A12
  sum_matrix_blocks(S[1],A,A,0,0,A1X,AX1,A1X,AX2,n2);
  // P2 = S2 x B22
  strassen_aux_optimal(P[2], (const float* const *const)S[1], B, 0, 0, 0, 0, B2X, BX2,n2);

  //P6
  // S7 = A12 - A22 (there was a bug here A21 - A22)
  sub_matrix_blocks(S[0],A,A,0,0,A1X,AX2,A2X,AX2,n2);

  // S8 = B21 + B22 (there was a bug here B21 - B22)
  sum_matrix_blocks(S[1],B,B,0,0,B2X,BX1,B2X,BX2,n2);

  // P6 = S7 x S8
  strassen_aux_optimal(P[3],(const float* const *const)S[0],(const float* const *const)S[1],
                       0, 0,0, 0,0, 0,n2);

  //P6
  // C11 = (P5 + P4) - P2 + P6
  sum_matrix_blocks(C, (const float* const *const)P[0], (const float* const *const)P[1],
                    C1X, CX1,0, 0, 0, 0,n2);

  sub_matrix_blocks(C, (const float* const *const)C, (const float* const *const)P[2],
                    C1X, CX1, C1X, CX1, 0, 0,n2);

  sum_matrix_blocks(C, (const float* const *const)C, (const float* const *const)P[3],
                    C1X, CX1,C1X, CX1, 0, 0,n2);

  //P1
  // S1 = B12 - B22
  sub_matrix_blocks(S[0],B,B,0,0,B1X,BX2,B2X,BX2,n2);
  // P1 = A11 x S1
  strassen_aux_optimal(P[3], A, (const float* const *const)S[0],0, 0,A1X, AX1, 0, 0,n2);

  // C12 = P1 + P2
  sum_matrix_blocks(C, (const float* const *const)P[2], (const float* const *const)P[3],
                    C1X, CX2, 0, 0, 0, 0,n2);

  //P1
  // S3 = A21 + A22
  sum_matrix_blocks(S[0],A,A,0,0,A2X,AX1,A2X,AX2,n2);

  // P3 = S3 x B11
  strassen_aux_optimal(P[2], (const float* const *const)S[0], B, 0, 0, 0, 0, B1X, BX1, n2);

  // C21 = P3 + P4
  sum_matrix_blocks(C, (const float* const *const)P[1], (const float* const *const)P[2],
                    C2X, CX1, 0, 0, 0, 0,n2);

  //P7
  // S9 = A11 - A21
  sub_matrix_blocks(S[0],A, A, 0,0,A1X,AX1,A2X,AX1,n2);

  // S10 = B11 + B12
  sum_matrix_blocks(S[1],B,B,0,0,B1X,BX1,B1X,BX2,n2);
  // P7 = S9 x S10
  strassen_aux_optimal(P[1], (const float* const *const)S[0], (const float* const *const)S[1],
                       0, 0, 0, 0,0, 0,n2);

  // C22 = P5 + P1 - P3 - P7
  sum_matrix_blocks(C, (const float* const *const)P[0], (const float* const *const)P[3],
                    C2X, CX2, 0, 0, 0, 0,n2);

  sub_matrix_blocks(C, (const float* const *const)C, (const float* const *const)P[2],
                    C2X, CX2, C2X, CX2, 0, 0,n2);

  sub_matrix_blocks(C, (const float* const *const)C, (const float* const *const)P[1],
                    C2X, CX2, C2X, CX2, 0, 0, n2);

  for (int i=0; i<2; i++) {
    deallocate_matrix(S[i], n2);
  }

  for (int i=0; i<4; i++) {
    deallocate_matrix(P[i], n2);
  }
}



void strassen(float **C,
      float const *const * const A, float const *const * const B, const size_t n)
{
  strassen_aux_optimal(C, A, B, 0, 0,  0, 0,  0, 0, n);
}
