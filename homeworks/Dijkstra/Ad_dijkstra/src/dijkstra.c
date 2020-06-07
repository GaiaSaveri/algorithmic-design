#include<stdlib.h>
#include<stdio.h>
#include<string.h>

#include "../include/dijkstra.h"

#define INFTY 999999

void init_sssp(Graph* G)
{
  for(size_t v=0; v<(G->n); v++)
  {
    (G->V[v]).idx = v;
    (G->V[v]).d = INFTY;
    (G->V[v]).pred = NULL;
  }

}

void relax(Node* s, Node* d, int weight)
{
  if(s->d+weight < d->d)
  {
    d->d = s->d+weight;
    d->pred = s;
  }
}

//---------> ARRAY BASED PRIORITY QUEUE <--------- //

void Dijkstra_Aq(Graph* g, int source_idx)
{
  init_sssp(g);
  Node* s = node(g, source_idx);
  s->d = 0;

  A_queue Q = build_Aq(g->V, g->n);
  A_queue* q = &Q;

  while(!is_empty_Aq(q))
  {
    Node* u = extract_minimum_Aq(q);
    int nhn = neighbours_number(g, u);
    Node** nh = neighbours(g, u, nhn);

    for(size_t i=0; i<nhn; i++)
    {
      relax(u, nh[i], weight(g, u, nh[i]));
    }

    free(nh);
  }
  delete_Aq(Q);
}

//---------> HEAP BASED PRIORITY QUEUE <--------- //

//define a specific order for nodes
int compare_dist(const void *a, const void *b)
{
  Node* s = (Node*)a;
  Node* d = (Node*)b;
  return leq_int((void*)&(s->d), (void*)&(d->d));
}

void Dijkstra_minheap(Graph* g, int source_idx)
{
  init_sssp(g);
  Node* s = node(g, source_idx);

  s->d = 0;
  binheap_type* H = build_heap(g->V, g->n, g->n, sizeof(Node), compare_dist);
  while(!(is_heap_empty(H)))
  {
    Node* u = (Node*)extract_min(H);
    int nhn = neighbours_number(g, u);
    Node** nh = neighbours(g, u, nhn);
    for(size_t i=0; i<nhn; i++)
    {
      relax(u, nh[i], weight(g, u, nh[i]));
    }
    heapify(H, 0);
    free(nh);
  }
  delete_heap(H);
}
