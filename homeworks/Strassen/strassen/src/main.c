#include <stdio.h>
#include <time.h>

#include "../include/test.h"
#include "../include/matrix.h"
#include "../include/strassen.h"
#include "../include/padding.h"

int main(int argc, char *argv[]) {

  size_t nn = 1 << 11;

  struct timespec start, stop;

  float **A = allocate_random_matrix(nn, nn);
  float **B = allocate_random_matrix(nn, nn);
  float **C0 = allocate_matrix(nn, nn);
  float **C1 = allocate_matrix(nn, nn);
  float **C2 = allocate_matrix(nn, nn);

  printf("\t\tTESTING SQUARE MATRICES\n");

  printf("\tStrassen\tStrassenOpt.\tNaive Alg.\tSame result (naive vs opt.)\n");

  for (size_t j = 1; j <= nn; j *= 2) {

    printf("%ld", j);

    printf("\t%lf", test(strassen_matrix_multiplication, C1, A, B, j));

    printf("\t%lf", test(strassen, C2, A, B, j));
    printf("\t%lf", test(naive_matrix_multiplication, C0, A, B, j));

    printf("\t%d\n", same_matrix((float const *const *const)C1,
                                 (float const *const *const)C2, j, j));
  }


deallocate_matrix(A, nn);
deallocate_matrix(B, nn);
deallocate_matrix(C0, nn);
deallocate_matrix(C1, nn);
deallocate_matrix(C2, nn);


int n;
int l;
int m;

printf("\nPlease insert dimensions of the matrix to be multiplied. \nC(nxm) =  A(nxl) X B(lxm)\n");

printf("\nn = ");
scanf("%u", &n);
printf("\nl = ");
scanf("%u", &l);
printf("\nm = ");
scanf("%u", &m);

float** A1 = allocate_random_matrix(n,l);
float** B1 = allocate_random_matrix(l,m);
float** C = allocate_matrix(n,m);
float**CC = allocate_matrix(n,m);

printf("\t\tTESTING PADDING\n");
printf("\t\t  Same result\n");

strassen_padding(CC,(float const *const *const)A1, n, l, (float const *const *const)B1, l, m);

naive_rectangular(C, (float const *const *const)A1, (float const *const *const)B1, n, l, m);

printf("\t\t       %d\n", same_matrix((float const *const *const)CC,
                             (float const *const *const)C, n, m));



deallocate_matrix(CC, n);
deallocate_matrix(A1, n);
deallocate_matrix(B1, l);
deallocate_matrix(C, n);


  return 0;
}
