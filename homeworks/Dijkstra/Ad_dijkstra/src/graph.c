#include<stdlib.h>
#include<stdio.h>

#include "../include/graph.h"
//#include "../include/queue.h"

#define INFTY 999999
//weight associated with the edge connecting node s (source) and node d (destination)
int weight(Graph* g, Node* s, Node* d)
{
  return g->E[s->idx*g->n + d->idx];
}

//number of neighbours of a given node
int neighbours_number(Graph* g, Node* node)
{
  int n = 0;
  for(size_t i=0; i<g->n; i++)
  {
    if((g->V)[i].idx == node->idx) continue;
    n+=(weight(g, node, &(g->V)[i]))<INFTY; //weight is less the INFTY
  }
  return n;
}

//returns the neighbours of a given node
Node** neighbours(Graph* g, Node* node, int n)
{
  Node** neigh = (Node**)malloc(sizeof(Node*)*n);
  int m = 0;
  for(size_t i=0; i<g->n; i++)
  {
    int w = weight(g, node, &(g->V)[i]);
    if(w<INFTY && (g->V)[i].idx != node->idx)
    {
      neigh[m++] = &(g->V)[i];
    }
  }
  return neigh;
}

//returns the node in the graph corresponding to the given idx
Node* node(Graph* g, int idx)
{
  Node* nd = g->V; //points to first node in the graph
  for(size_t i=0; i<g->n; i++)
  {
    if(nd[i].idx==idx) { return &nd[i]; }
  }
  return NULL; //avoid warnings
}

//utility function to show graph's adjacency matrix
//distance of every node from the source and predecessor
void graph_properties(Graph* g)
{
    int size = g->n;
    printf("\nthe number of nodes in the graph is: %d\nthe adjacency matrix of the graph is:\n", g->n);

    for(int i=0; i<size; i++)
    {
        for(int j=0; j<size; j++)
        {
            if(g->E[size*i+j]==INFTY)
                printf("NIL\t");
            else
                printf("%d\t", g->E[size*i+j]);
        }
        printf("\n");
    }

    for(int i=0; i<size; i++)
    {
        Node node = g->V[i];

        printf("\nnode %d: dist=%d", node.idx, node.d);
        if((node.pred)!=NULL)
            printf(", pred=%d", node.pred->idx);
        else
            printf(", ");

    }

    printf("\n");
}
