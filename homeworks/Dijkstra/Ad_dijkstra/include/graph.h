#ifndef __GRAPH__

#include "queue.h"

typedef struct Graph {
  Node* V; //set of the nodes of the graph
  int* E; //adjacency matrix ("weight" version)
  int n; //number of nodes
} Graph;

int weight(Graph* g, Node* s, Node* d);

int neighbours_number(Graph* g, Node* node);

Node** neighbours(Graph* g, Node* node, int n);

Node* node(Graph* g, int idx);

void graph_properties(Graph* g);

#endif //__GRAPH__
