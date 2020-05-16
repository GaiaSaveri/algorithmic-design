#include<stdlib.h>

#include "../include/queue.h"

//build array queue
A_queue build_Aq(Node* nodes, size_t n)
{
  Node** ptrs = (Node**)malloc(sizeof(Node*)*n);
  for(size_t i=0; i<n; i++)
  {
    ptrs[i] = &(nodes[i]);
  }
  A_queue q;
  q.n = n;
  q.A = ptrs;
  return q;
}

//swaps two nodes inside the array queue
void Aq_swap(Node **a, Node **b)
{
  Node *temp = *a;
  *a = *b;
  *b = temp;
}

//check if the queue is empty
int is_empty_Aq(A_queue* q) { return q->n == 0; }
//delete the queue
void delete_Aq(A_queue q) { free(q.A); }
//delete a node in the queue given the index
void dequeue_Aq(A_queue* q, int i)
{
  Aq_swap(&q->A[i], &q->A[q->n-1]); //swap with last element
  q->n--;
}
//extract the node having minimum candidate distance
Node* extract_minimum_Aq(A_queue *q)
{
  int min = q->A[0]->d;
  int min_idx = 0;
  for(size_t i=0; i<q->n; i++)
  {
    if(q->A[i]->d < min)
    {
      min = q->A[i]->d;
      min_idx = i;
    }
  }

  Node* node = q->A[min_idx];
  dequeue_Aq(q, min_idx);
  return node;
}
