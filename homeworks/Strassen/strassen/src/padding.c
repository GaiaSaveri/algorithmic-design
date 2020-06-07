#include <stdio.h>
#include <stdlib.h>
#include <math.h> //pow
#include <string.h> //memcpy

#include "padding.h"
#include "matrix.h"
#include "strassen.h"

/*
 * This function performs the naive matrix multiplication between two possibly rectangular matrices.
 */
void naive_rectangular(float **C, float const *const *const A,
                                float const *const *const B,
                                const size_t A_r, //number of rows of A
                                const size_t A_c, //number of columns of A = number of rows of B
                                const size_t B_c //number of columns of B
                              )
{
  for (size_t y = 0; y < A_r; y++) {
    for (size_t x = 0; x < B_c; x++) {
      float value = 0.0;
      for (size_t z = 0; z < A_c; z++) {
        value += A[y][z]*B[z][x];
      }

      C[y][x] = value;
    }
  }
}

/*
 * This function finds the smallest power of two bigger than the argument r.
 */
int pow_two(const int r)
{
  int i = 0;
  int p = pow(2,i);
  while(p<r)
  {
    p = pow(2,i);
    i++;
  }
  return p;
}

/*
 * This function finds the smallest power of two bigger than the greatest of its arguments.
 */
int pow_two_rectangular(const size_t A_r, const size_t A_c)
{
  if(A_r>=A_c)
  {
    return pow_two(A_r);
  }
  else
  {
    return pow_two(A_c);
  }
}

/*
 * This function performs the sum between two matrices.
 */
void sum_matrix(float **C, float const *const *const A, float const *const *const B, const size_t n)
{
    for (size_t y = 0; y < n; y++) {
      for (size_t x = 0; x < n; x++) {
          C[y][x] = A[y][x] + B[y][x];
      }
    }
}

/*
 * This function is used to eliminate from the matrix obtained as result of the "padded" product
 * the superfluous zeros. In this way the result of such a product has the expected dimensions,
 * and all the auxiliary matrices used are freed, to avoid memory leaks.
 */
void fix(float** CC, float** C, const size_t r_b, const size_t c_b, const size_t r_a)
{
  for(size_t i=0; i<r_b; i++)
  {
    memcpy(CC[i], C[i], sizeof(float)*c_b);
  }
  deallocate_matrix(C, r_a);
}

/*
 * This function is used to immerse a possibly rectangular matrix into a squared matrix,
 * whose size is a power of two. All the entries added to original matrix are filled with 0.
 */
float** padding_submatrix(float const *const *const A, //matrix to be squared
                          const size_t A_f_row, const size_t A_f_col, //starting indexes
                          const size_t A_r, const size_t A_c, //number of rows and columns to be copied
                          const size_t n //dimension of the resulting matrix
                        )
{
  float**Q = allocate_matrix(n,n);
  for(size_t i=0; i<A_r; i++)
  {
    for(size_t j=0; j<A_c; j++)
    {
      Q[i][j]=A[A_f_row+i][A_f_col+j];
    }
  }
  return Q;

}

/*
 * This function is used to decompose a matrix in blocks.
 * Our aim is to perform the matrix product AxB, with A and B possibly rectangular.
 * If we assume that the number of rows of A is less than the number of its columns,
 * then we can think of A as divided in square blocks (with side equal to the number of rows of A)
 * and divide B accordingly, so that the pairwise product between blocks can be performed.
 * Obviously, following the reasoning above, the matrix A has to be cut along its columns
 * (so that each block will have the same number of rows of the original matrix)
 * and the matrix B has to be cut along its rows (so that each block will have the same number of rows
 * of the original matrix.)
 */
float*** divide_matrix(float const *const *const A, //matrix to be divided
                       const size_t A_r, const size_t A_c, //dimensions of A
                       const size_t B_r, const size_t B_c, //dimensions of each block
                       const size_t n, //dimension of the padded block
                       int by_row // 1--> matrix has to be divided by row, 0--> matrix has to be divided by column
                     )
{
  int n_blocks = 0;
  if(by_row==1)
  {
    n_blocks = (A_r%B_r>0)? A_r/B_r+1 : A_r/B_r;
    unsigned int start[n_blocks]; //row at which the i-th block starts
    for(size_t i=0; i<n_blocks; i++)
    {
      start[i]=i*B_r;
    }
    //array of n_blocks matrix
    float*** Pad_A = (float ***)malloc(sizeof(float **) * n_blocks);

    for (size_t i = 0; i < n_blocks; i++) {
        //lowest block
        if(i==(n_blocks-1))
        {
          int end_idx = (start[i]+(B_r-1)>(A_r-1))? A_r-1 : start[i]+(B_r-1);
          Pad_A[i] = padding_submatrix(A, start[i], 0, (end_idx-start[i]+1), B_c, n);

        }
        else
        {
          Pad_A[i] = padding_submatrix(A, start[i], 0, B_r, B_c, n);
        }
    }
    return Pad_A;
  }

  else //by_row == 0
  {
    n_blocks = (A_c%B_c>0)? A_c/B_c+1 : A_c/B_c;
    unsigned int start[n_blocks]; //column at which the i-th block starts
    for(size_t i=0; i<n_blocks; i++)
    {
      start[i]=i*B_c;
    }
    //array of n_blocks matrix
    float*** Pad_A = (float ***)malloc(sizeof(float **) * n_blocks);

    for (size_t i = 0; i < n_blocks; i++)
    {
      //rightmost block
      if(i==(n_blocks-1))
      {
        int end_idx = (start[i]+(B_c-1)>(A_c-1))? A_c-1 : start[i]+(B_c-1);
        Pad_A[i] = padding_submatrix(A, 0, start[i], B_r, (end_idx-start[i]+1), n);
      }
      else
      {
        Pad_A[i] = padding_submatrix(A, 0, start[i], B_r, B_c, n);
      }
    }
    return Pad_A;
  }
}

/*
 * This function is used to perfom the "block-wise" product.
 * If we divide the two matrices as described above, than the product will
 * be given by the sum of all the product of the corresponding blocks.
 * All the auxiliary matrices used inside the function are deallocated before the
 * function returns.
 */
void lego_bricks_multiplication(float** C, float const *const *const A,
                                   const size_t A_r, const size_t A_c, //dimensions of A
                                   float const *const *const B,
                                   const size_t B_r, const size_t B_c //dimensions of B
                                   )
{
  //divide A and B
  //first of all we need the dimension of the blocks of A
  size_t BA_r = A_r; //A will have squared blocks (except the last)
  size_t BB_r = BA_r;
  size_t BB_c = B_c;

  size_t n = (BA_r>BB_c)? pow_two(BA_r) : pow_two(BB_c); //size of the padded blocks

  float*** Pad_A = divide_matrix(A, A_r, A_c, BA_r, BA_r, n, 0);
  float*** Pad_B = divide_matrix(B, B_r, B_c, BB_r, BB_c, n, 1);

  //by contruction length(Pad_A) = length(Pad_B) since they are built in the common dimension
  int n_blocks = A_c/BA_r+1;

  //now we need to multiply correspondent blocks and sum the pairwise results
  float** D = allocate_matrix(n,n); //this will be returned
  float** CC = allocate_matrix(n,n); //this is auxiliary, for partial results

  for(size_t i=0; i<n_blocks; i++)
  {
    naive_matrix_multiplication(CC, (float const *const *const)Pad_A[i], (float const *const *const)Pad_B[i], n); //STRASSEN
    sum_matrix(D,(float const *const *const)D,(float const *const *const)CC,n);
  }

  for (size_t i = 0; i < n_blocks; i++) {
    deallocate_matrix(Pad_A[i], n);
    deallocate_matrix(Pad_B[i], n);
  }
  free(Pad_A);
  free(Pad_B);
  deallocate_matrix(CC, n);
  fix(C, D, A_r, B_c, n);
}

/*
 * This function is used to perform the product between two square matrices,
 * whose side is not a power of two.
 */
void squared_matrix_padding(float** C, float const *const *const A,
                               const size_t A_r, const size_t A_c,  //dimensions of A
                               float const *const *const B,
                               const size_t B_r, const size_t B_c //dimensions of B
                               )
{
    int p = pow_two(A_r);
    if(p==A_r) //dimension is a power of two
    {
      strassen(C, A, B, A_r);
    }
    else //the dimension is not a power of two
    {
      float** Pad_A = padding_submatrix(A, 0, 0, A_r, A_c, p);
      float** Pad_B = padding_submatrix(B, 0, 0, B_r, B_c, p);
      float **CC = allocate_matrix(p,p);
      strassen_matrix_multiplication(CC, (float const *const *const)Pad_A, (float const *const *const)Pad_B, p);
      deallocate_matrix(Pad_A, p);
      deallocate_matrix(Pad_B, p);
      fix(C, CC, A_r, B_c, p);
    }
}

/*
 * This function is used to perform the product of two matrices of any size
 * (provided that their dimensions are compatible for the product).
 */
void strassen_padding(float** C, float const *const *const A,
                         const size_t A_r, const size_t A_c, //dimension of A
                         float const *const *const B,
                         const size_t B_r, const size_t B_c //dimensions of B
                         )
{
  if(A_r==A_c) //A is a squared matrix
  {
    if(B_c==B_r) //B is squared too
    {
      squared_matrix_padding(C, A, A_r, A_c, B, B_r, B_c);
    }
    else //B is not squared
    {
      if(A_r>B_c)
      {
        int p = pow_two_rectangular(B_r, B_c);
        float** B_sq = padding_submatrix(B, 0, 0, B_r, B_c, p);
        squared_matrix_padding(C, A, A_r, A_c, (float const *const *const)B_sq, p, p);
        deallocate_matrix(B_sq, p);
      }
      else //B_c>A_r
      {
        int p = pow_two(B_c);
        float** A_sq = padding_submatrix(A, 0, 0, A_r, A_c, p);
        float** B_sq = padding_submatrix(B, 0, 0, B_r, B_c, p);
        float** CC = allocate_matrix(p,p);
        strassen_matrix_multiplication(CC, (float const *const *const)A_sq, (float const *const *const)B_sq, p);
        fix(C, CC, A_r, B_c, p);
        deallocate_matrix(A_sq,p);
        deallocate_matrix(B_sq,p);
      }
    }
  }
  else //A is not squared
  {
    if(A_r>A_c) //number of rows of A is bigger than the number of its columns
    {
      int p = pow_two_rectangular(A_r, B_c);
      float** A_sq = padding_submatrix(A, 0, 0, A_r, A_c, p);
      float** B_sq = padding_submatrix(B, 0, 0, B_r, B_c, p);
      float** CC = allocate_matrix(p,p);
      strassen_matrix_multiplication(CC, (float const *const *const)A_sq, (float const *const *const)B_sq, p);
      fix(C, CC, A_r, B_c, p);
      deallocate_matrix(A_sq,p);
      deallocate_matrix(B_sq,p);
    }
    else //number of rows of A is smaller than the number of its columns
    {
      lego_bricks_multiplication(C, A, A_r, A_c, B, B_r, B_c);
    }
  }
}
