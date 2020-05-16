#include "select.h"

unsigned int select_index(void *A, const unsigned int n, const unsigned int i, const size_t elem_size,  total_order leq);

//median of medians algorithm
unsigned int select_aux(void *A, const unsigned int n, const size_t elem_size, total_order leq)
{
  unsigned int left = 0;
  unsigned int right = n-1;

  if(n<10) //base case
  {
    insertion_sort(A, n, elem_size, leq);
    return n/2;
  }

  int n_chunks = n/5;

  for(size_t i = 0; i<n_chunks; i++)
  {

    int end = (5*i+left+4<right)? 5*i+left+4 : right; //ending index current chunk
    int begin = (5*i+left); //starting index current chunk
    int n = end - begin + 1;
    insertion_sort(A+(begin)*elem_size, n, elem_size, leq);
    int median = (5*i+left+2 < right)? 5*i+left+2 : right; //5 elem chunks mostly
    swap(A+(left+i)*elem_size, A+median*elem_size, elem_size);

  }

  return select_index(A, n_chunks, n_chunks/2, elem_size, leq);
}

unsigned int select_index(void *A, const unsigned int n,
                          const unsigned int i,
                          const size_t elem_size,
                          total_order leq)
{
  if(n<=10)
  {
    insertion_sort(A, n, elem_size, leq);
    return i;
  }

  unsigned int j = select_aux(A, n, elem_size, leq);
  pair_type k = tri_partition(A, elem_size, 0, n-1, j, leq);
  if(i<k.first)
  {
    if(k.first>0)
    {
      return select_index(A, k.first-1, i, elem_size, leq);
    }
    else { return k.first;}
  }

  if(i>k.second)
  {
    if(k.second<n)
    {
      return select_index(A+k.second*elem_size, n-k.second-1, i, elem_size, leq);
    }
    else { return k.second; }
  }

  return i;
}

void quick_sort_select_aux(void *A, size_t l, size_t r, const size_t elem_size, total_order leq)
{

  while(l<r)
  {
    unsigned int j = l + select_aux(A+l*elem_size, r-l, elem_size, leq);
    pair_type k = tri_partition(A, elem_size, l, r-1, j, leq);
    quick_sort_select_aux(A, l, k.first, elem_size, leq);
    l = k.second+1;
  }

}

void quick_sort_select(void *A, const unsigned int n,
                       const size_t elem_size,
                       total_order leq)
{
   quick_sort_select_aux(A, 0, n, elem_size, leq);
}
