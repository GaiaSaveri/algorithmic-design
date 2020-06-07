#include "quick_sort.h"
#include "swap.h"
#include "total_order.h"


//fuction which perfrms the partition
int partition(void *A, const size_t elem_size, size_t left, size_t right, size_t pivot_idx, total_order leq)
{
  swap(A+(pivot_idx)*elem_size, A+(left)*elem_size, elem_size); //swap pivot with the first element of the array
  size_t p_idx = left;
  size_t i = left + 1;
  size_t j = right; //right-1?
  while(i<=j)
  {
    if(leq(A+(p_idx)*elem_size, A+(i)*elem_size)) //A[i]>pivot
    {
      swap(A+i*elem_size, A+j*elem_size, elem_size); //place it in G
      j--; //increase G's size
    }
    else { i++; } //A[i] is already in S
  }

  swap(A+p_idx*elem_size, A+j*elem_size, elem_size); //place the pivot between S and G
  return j;
}

//tail recursion
void quick_sort_aux(void *A, void *A_l, void* A_r, const size_t elem_size, total_order leq)
{
  int l = (A_l - A)/(elem_size); //left index
  int r = (A_r - A)/(elem_size); //right index
  while(l<r)
  {
    int p = partition(A, elem_size, l, r, l, leq);
    quick_sort_aux(A, A+l*elem_size, A+(p-1)*elem_size, elem_size, leq);
    l = p+1;
  }
}

void quick_sort(void *A, const unsigned int n,
                const size_t elem_size,
                total_order leq)
{
  quick_sort_aux(A, A, A+(n-1)*elem_size, elem_size, leq);
}
