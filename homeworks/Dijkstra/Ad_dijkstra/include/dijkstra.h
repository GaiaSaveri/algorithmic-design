#ifndef __DIJKSTRA__

#include "graph.h"
#include "binheap.h"

void init_sssp(Graph* G);

void relax_Aq(A_queue* q, Node* s, Node* d, int weight);

void relax_minheap(binheap_type* H, Node* s, Node* d, int weight);

void Dijkstra_Aq(Graph* g, int source_idx);

void Dijkstra_minheap(Graph* g, int source_idx);

void show_graph(Graph* g);

#endif //__DIJKSTRA__
