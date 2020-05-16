#include"binheap.h"
#include <string.h> //memcpy
#include <stdio.h> //printf

//these 3 macros refers to nodes in pos_key now, i.e. node is an index in pos_key
#define PARENT(node) ((node-1)/2)
#define LEFT_CHILD(node) (2*(node)+1)
#define RIGHT_CHILD(node) (2*(node+1))

#define VALID_NODE(H, node) ((H)->num_of_elem>(node))

//these two refers to A
#define ADDR(H, node) ((H)->A+(node)*(H)->key_size)
#define INDEX_OF(H, addr) ((addr-((H)->A))/((H)->key_size))

//this refers to pos_key
//heap_index is the index of the node of interest in the heap structure (not in A)
#define ADDR_KEY(H, heap_index) (&(H->pos_key[heap_index]))
//this refers to rev_pos
#define ADDR_REV(H, pos_key_index) (&(H->rev_key[pos_key_index]))

int is_heap_empty(const binheap_type *H) //not depending on the storing method
{
  return H->num_of_elem==0;
}

const void *min_value(const binheap_type *H)
{
  if(is_heap_empty(H)) {return NULL;}
  //if the heap is not empty, the min is in the root
  //the key of the root can be found in A[pos_key[0]]
  return ADDR(H, *(H->pos_key)); //H->pos_key is the address to the first cell of pos_key, aka pos_key[0]
}

//every time we swap, we need to fix also the elements of rev_key
//in order to mantain the overall coherence of the structure
void fix_rev_key(binheap_type *H)
{
  //recall that pos_key contains indexes of pos_key
  //in particular it holds rev_pos[i] := index of pos_key handling the i-th element of A
  for(int i=0; i<(H)->num_of_elem; i++) //looping pos_key
  {
    unsigned int index = H->pos_key[i]; //content of the i-th cell of pos_key
    //put myself in the index-th cell of rev_key
    //and I put there the index of pos_key accordingly
    H->rev_key[index] = i;
  }
}

//we swap elements of pos_key (aka indeces of A) instead of swapping directly elements of A
void swap_pos_key(binheap_type *H, unsigned int pos_key_a, unsigned int pos_key_b)
{
  //pos_key
  unsigned int *index_a = ADDR_KEY(H, pos_key_a);
  unsigned int *index_b = ADDR_KEY(H, pos_key_b);
  unsigned int *tmp_index = malloc(sizeof(unsigned int));
  //swap in pos_key
  memcpy(tmp_index, index_a, sizeof(unsigned int));
  memcpy(index_a, index_b, sizeof(unsigned int));
  memcpy(index_b, tmp_index, sizeof(unsigned int));
  free(tmp_index);
  //in order to maintain the coherence of the structure, everytime we swap in pos_key
  //we need to adjust rev_key
  fix_rev_key(H);
}

void heapify(binheap_type *H, unsigned int node) //node is an index in pos_key now
{
  //heapify works on 3 nodes: current node and its children
  unsigned int dst_node=node; //node containing the min among node and its children --> in pos_key
  unsigned int child; //in pos_key
  unsigned int A_child; //index in A, namely content of pos_key[child]
  unsigned int A_node; //index in A, namely content of pos_key[node]
  unsigned int A_dst; //auxiliary, index in A corresponding to the current dst node
  do {
    node = dst_node; //update node we are interested in
    //i.e. the minimum identified the previous step
    A_node = *(ADDR_KEY(H,node));
    A_dst = *(ADDR_KEY(H,dst_node));
    //identify minimum among current node and its children
    child = LEFT_CHILD(node); //index of the left child of node in pos_key
    A_child = *(ADDR_KEY(H,child));
    if(VALID_NODE(H, child) &&
       H->leq(ADDR(H, A_child), ADDR(H, A_node))) {
         //this child is a good candidate to be placed as root
         dst_node = child;
         A_dst = *(ADDR_KEY(H,dst_node));
       }

    child = RIGHT_CHILD(node); //index of the right child of node in pos_key
    A_child = *(ADDR_KEY(H,child));
    //need to repeat this operation also for the left child
    if(VALID_NODE(H, child) &&
       H->leq(ADDR(H, A_child), ADDR(H, A_node))) {
         //this child is a good candidate to be placed as root
         if(H->leq(ADDR(H, A_child), ADDR(H, A_dst))) //comparison with the current node set as dst
         {
         dst_node = child; //index in pos_key
         }
       }
    if(dst_node!=node) //heap structure is not preserved
    {
      swap_pos_key(H, dst_node, node); //swapping in pos_key and rev_key
    }
  } while(dst_node != node);
}

const void *extract_min(binheap_type *H)
 {
    if(is_heap_empty(H)) {return NULL;}
    //if the heap is not empty, recall that the min is in the root
    //swapping the key among the root and the right-most leaf of the last level
    //A[0] <-> A[num_of_elem-1]
    //in this structure we need to swap pos_key[0]<->pos_key[num_of_elem-1]
    swap_pos_key(H, 0, ((H)->num_of_elem)-1);
    unsigned int index_min = *(ADDR_KEY(H, ((H)->num_of_elem)-1)); //index of the minimum in A
    void *min = ADDR(H, index_min); //this should be returned
    //delete min
    if(H->num_of_elem>1)
    {
     H->num_of_elem--;
     heapify(H, 0); //preserve the whole heap structure
    }
    else {H->num_of_elem--;} //if the heap is left with only one element--> no heapify
    return min;
    //NOTE: the minimum key is actually still in A, even if no element of pos_key is managing it
    //some attempts have been made to fix this, all unsuccessful :(
 }

const void *find_the_max(void *A,
                         const unsigned int num_of_elem,
                         const size_t key_size,
                         total_order_type leq)
{
  if(num_of_elem==0) {return NULL;}
  //assume max is in the first cell of the array
  const void *max_value = A;
  for(const void *addr=A+key_size; addr!=A+num_of_elem*key_size; addr+=key_size)
  {
    if(!leq(addr, max_value))
    {
      max_value = addr;
    }
  }
  return max_value;
}

binheap_type *build_heap(void *A,
                         const unsigned int num_of_elem,
                         const unsigned int max_size,
                         const size_t key_size,
                         total_order_type leq)
{
  binheap_type *H = (binheap_type *)(malloc)(sizeof(binheap_type));
  H->A = A;
  H->pos_key = (unsigned int *)malloc(max_size*(sizeof(unsigned int)));
  H->rev_key = (unsigned int *)malloc(max_size*(sizeof(unsigned int)));
  //initialize pos_key with the first num_of_elem-1 integers
  for(unsigned int i=0; i<num_of_elem; i++)
  {
    H->pos_key[i] = i;
  }
  H->num_of_elem = num_of_elem;
  H->max_size = max_size;
  H->key_size = key_size;
  H->leq = leq;
  H->max_order_value = malloc(key_size); //void*, no type casting of malloc

  if(num_of_elem==0)
  {
    return H; //nothing to fix
  }

  const void* value = find_the_max(A, num_of_elem, key_size, leq);
  memcpy(H->max_order_value, value, key_size);
  //fix the heap property from the second last level, up to the root
  for(unsigned int  i=num_of_elem/2; i>0; i--)
  {
    heapify(H, i);
  }
  heapify(H, 0); //not inlcuded in the loop
  //with heapify I fixed indexes of pos_key
  //still need to fix the ones of rev_key
  fix_rev_key(H);
  return H;
}

void delete_heap(binheap_type *H)
{
    free(H->max_order_value);
    free(H->pos_key);
    free(H->rev_key);
    free(H);
}

const void *decrease_key(binheap_type *H, void *node, const void *value)
{
    //node is the address of an element in A
    unsigned int node_idx = INDEX_OF(H, node); //index in A of the address in input
    //need to individuate the node in pos_key --> need to pass through rev_key
    //rev_key[node_idx] --> index of rev key I'm interested in
    unsigned int *pos_key_addr = ADDR_REV(H, node_idx);
    unsigned int pos_key_index = *(pos_key_addr); //index of the starting node in pos_key

    //if node does not belong to H, or value>=node we return null
    if(!(VALID_NODE(H, node_idx) || !(H->leq(value, node))))
    {
      return NULL;
    }

    memcpy(node, value, H->key_size);
    //we may end up in a situation in which heap property doesn't hold anymore

    if((H)->num_of_elem>1){ //if we do this machinery when num_of_elem = 1, we get segfault
    unsigned int parent_idx = PARENT(pos_key_index); //index of the parent in pos_key
    unsigned int *parent_addr = ADDR_KEY(H, parent_idx); //address in pos_key of the parent

    void *parent = ADDR(H, *(parent_addr)); //address of the parent in A

    while((pos_key_index!=0)&(!H->leq(parent, node)))
    {
    //if the parent has value which is greater than that of the node
    //swap parent and node keys
    swap_pos_key(H, pos_key_index, parent_idx);
    //moving up the check --> focus on the node's parent
    node = ADDR(H, *(parent_addr));
    pos_key_index = parent_idx;

    parent_idx = PARENT(pos_key_index);  //index in pos_key
    if(pos_key_index!=0)
    {
    parent_addr = ADDR_KEY(H, parent_idx); //address in pos_key
    parent = ADDR(H, *(parent_addr)); //address in A
    }
    }
   fix_rev_key(H); //necessary if no swap has been performed yet
   }
   return node;
}

const void *insert_value(binheap_type *H, const void *value)
{
  if(H->max_size == H->num_of_elem) //heap already full
  {
    return NULL;
  }
  //if the value is greater than max_value
  if(H->num_of_elem == 0 || !H->leq(value, H->max_order_value))
  {
    memcpy(H->max_order_value, value, H->key_size);
  }
  void *new_node_addr = ADDR(H, H->num_of_elem); //get position of the new node in A
  //need to fix also pos_key
  unsigned int *new_node_pos_key = ADDR_KEY(H, H->num_of_elem); //address of the cell where insert in pos_key
  unsigned int index_in_pos_key = INDEX_OF(H, new_node_addr); //element to be inserted in pos_key
  memcpy(new_node_addr, H->max_order_value, H->key_size); //max order value is void*
  memcpy(new_node_pos_key, &index_in_pos_key, sizeof(unsigned int));
  //fix rev_key accordingly
  H->rev_key[H->num_of_elem] = H->num_of_elem;
  H->num_of_elem++;

  //decrease the key of this new node to the desired value
  return decrease_key(H, new_node_addr, value);
}

void print_heap(const binheap_type *H, void (*key_printer)(const void *value))
{
    unsigned int next_level_node = 1;
    for(unsigned int node = 0; node < H->num_of_elem; node++) //now this index runs in pos_key
    {
      if(node == next_level_node)
      {
        printf("\n");
        next_level_node = LEFT_CHILD(node); //index in pos_key
      }
      else
      {
        printf("\t");
      }
      unsigned int *A_pos = ADDR_KEY(H, node); //retrieve the actual key to be printed
      key_printer(ADDR(H, *(A_pos)));
    }

    printf("\n");
}
