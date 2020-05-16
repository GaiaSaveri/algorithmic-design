#include "bubble_sort.h"

void bubble_sort(void *A, const unsigned int n,
                 const size_t elem_size,
                 total_order leq)
{
  for(size_t i = n-1; i>=1; i--)
  {
    for(size_t j = 0; j<i; j++)
    {
      if(!leq(A+j*elem_size,A+(j+1)*elem_size))
      {
        swap(A+j*elem_size, A+(j+1)*elem_size, elem_size);
      }
    }
  }
}
