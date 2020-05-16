#include "insertion_sort.h"
#include "swap.h"

void insertion_sort(void *A, const unsigned int n,
                    const size_t elem_size,
                    total_order leq)
{
  void* j;
  for (void *addr=A+elem_size; addr!=A+(n)*elem_size;
    addr+=elem_size)
  {
    j = addr;
    while((j!=A) && leq(j, j-elem_size)) //need to insert A[i] into the sorted sequence
    {
      swap(j-elem_size, j, elem_size);
      j-=elem_size;
    }
  }
}
