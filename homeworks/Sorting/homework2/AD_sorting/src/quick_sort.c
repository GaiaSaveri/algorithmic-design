#include<stdio.h>
#include "quick_sort.h"
#include "swap.h"


pair_type tri_partition(void *A, const size_t elem_size, size_t left, size_t right, size_t pivot_idx, total_order leq)
{
  swap(A+left*elem_size, A+pivot_idx*elem_size, elem_size);
  pivot_idx = left;
  left++;
  unsigned int same = 0;

  while(left<=right)
  {

    if(!leq(A+pivot_idx*elem_size, A+left*elem_size) || !leq(A+left*elem_size, A+pivot_idx*elem_size))
    {
      if(leq(A+left*elem_size, A+pivot_idx*elem_size))
      {
        swap(A+(left)*elem_size, A+(pivot_idx-same)*elem_size, elem_size);
        pivot_idx = left;
        left++;
      }
      else
      {
        swap(A+left*elem_size, A+right*elem_size, elem_size);
        right--;
      }
    }

    else // =
    {
      pivot_idx = left;
      left++;
      same++;
    }
  }
  swap(A+(pivot_idx)*elem_size, A+right*elem_size, elem_size);
  pair_type k;
  k.second= right-same;
  k.first = right; 
  return k;
}

void quick_sort_aux(void *A, void *A_l, void* A_r, const size_t elem_size, total_order leq)
{
  int l = (A_l - A)/(elem_size); //left index
  int r = (A_r - A)/(elem_size); //right index

  while(l<r)
  {
    pair_type p = tri_partition(A, elem_size, l, r, l, leq);
    quick_sort_aux(A, A+l*elem_size, A+(p.first-1)*elem_size, elem_size, leq);
    l = p.second+1;
  }
}

void quick_sort(void *A, const unsigned int n,
                const size_t elem_size,
                total_order leq)
{
    quick_sort_aux(A, A, A+(n-1)*elem_size, elem_size, leq);
}
