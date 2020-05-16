#ifndef __QUEUE__

//array queue
typedef struct Node {
  int idx; //identifier for the node
  int d; //distance from the source
  struct Node *pred; //predecessor of the node
} Node;

//very simple implementation of array based priority queue
typedef struct A_queue {
  Node **A;
  size_t n; //number of nodes in the queue
} A_queue;

A_queue build_Aq(Node* nodes, size_t n);

int is_empty_Aq(A_queue* q);

Node* extract_minimum_Aq(A_queue *q);

void delete_Aq(A_queue q);

#endif //__QUEUE__
