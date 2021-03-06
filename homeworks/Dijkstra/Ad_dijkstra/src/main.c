#include<stdlib.h>
#include<stdio.h>
#include<time.h>

#include"../include/dijkstra.h"

#define INFTY 999999

double elapsed_time(const struct timespec start,
            const struct timespec stop)
{
  return (stop.tv_sec-start.tv_sec) +
                   (stop.tv_nsec-start.tv_nsec)/1E9;
}

int main() {

//---------> CORRECTNESS TEST <--------- //
  //create the example showed in the slides
  int size = 6;
  int* adj_mat = (int*)malloc(sizeof(int)*size*size);

  for(size_t i=0; i<size*size; i++)
  {
    adj_mat[i] = INFTY;
  }

  adj_mat[size*0+1] = 1;
  adj_mat[size*0+2] = 5;
  adj_mat[size*1+5] = 15;
  adj_mat[size*2+3] = 2;
  adj_mat[size*3+4] = 1;
  adj_mat[size*4+1] = 3;

  Graph g;
  g.n = size;
  g.E = adj_mat;
  g.V = (Node*)malloc(sizeof(Node)*size);

  Graph gg;
  gg.n = size;
  gg.E = adj_mat;
  gg.V = (Node*)malloc(sizeof(Node)*size);

  printf("\nINPUT GRAPH'S PROPERTIES\n");

  graph_properties(&g);

  printf("\nARRAY-BASED DIJKSTRA'S ALGORITHM\n");

  Dijkstra_Aq(&g, 0);
  graph_properties(&g);

  printf("\nHEAP-BASED DIJKSTRA'S ALGORITHM\n");

  Dijkstra_minheap(&gg, 0);
  graph_properties(&gg);

  free(g.V);
  free(gg.V);
  free(adj_mat);


//---------> PERFORMANCE TEST <--------- //

  struct timespec start, stop;

  int size1 = 30000;
  int* adj_mat1 = (int*)malloc(sizeof(int)*size1*size1);

  for(size_t i=0; i<size1; i++)
  {
    for(size_t j=0; j<size1; j++)
    {
      if((i+j)<2500)
        adj_mat1[i*size1+j] = i+j;
      if((i+j)>300)
        adj_mat1[j*size1+i] = i+j+5;
    }
  }

  Graph g1;
  g1.E = adj_mat1;
  g1.V = (Node*)malloc(sizeof(Node)*size1);

  Graph g2;
  g2.E = adj_mat1;
  g2.V = (Node*)malloc(sizeof(Node)*size1);

  printf("\n\n\tTESTING DIJKSTRA\n");

  printf("size\tArray\t\tHeap\n");
  for(size_t i=0; i<10; i++)
  {
    g1.n = size1/(1<<(9-i));
    g2.n = size1/(1<<(9-i));

    printf("\n%d", g1.n);

    clock_gettime(CLOCK_REALTIME, &start);
    Dijkstra_Aq(&g1, 0);
    clock_gettime(CLOCK_REALTIME, &stop);

    printf("\t%lf", elapsed_time(start, stop));

    clock_gettime(CLOCK_REALTIME, &start);
    Dijkstra_minheap(&g2, 0);
    clock_gettime(CLOCK_REALTIME, &stop);

    printf("\t%lf", elapsed_time(start, stop));

  }
  free(adj_mat1);
  free(g1.V);
  free(g2.V);
  return 0;
}
