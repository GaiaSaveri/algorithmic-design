#include <binheap.h>
#include <string.h> //memcpy
#include <stdio.h> //printf

//these 3 macros refers to nodes in pos_key now
//node is an index in pos_key
#define PARENT(node) ((node-1)/2)
#define LEFT_CHILD(node) (2*(node)+1)
#define RIGHT_CHILD(node) (2*(node+1))

#define VALID_NODE(H, node) ((H)->num_of_elem>(node))

//these two refers to A
#define ADDR(H, node) ((H)->A+(node)*(H)->key_size)
#define INDEX_OF(H, addr) ((addr-((H)->A))/((H)->key_size))

int is_heap_empty(const binheap_type *H) //not depending on the storing method
{
  return H->num_of_elem==0;
}

const void *min_value(const binheap_type *H)
{
  //return the minimum value stored as a key in our binheap
  //if the heap is empty, NULL should be returned
  if(is_heap_empty(H)) {return NULL;}
  //if the heap is not empty, the min is in the root
  return ADDR(H,0);
}


void swap_keys(binheap_type *H, unsigned int n_a, unsigned int n_b)
{
  void *p_a = ADDR(H, n_a);
  void *p_b = ADDR(H, n_b);
  void *tmp = malloc(H->key_size); //temporary space to perform the swap

  memcpy(tmp, p_a, H->key_size); //copy p_a in tmp
  memcpy(p_a, p_b, H->key_size);
  memcpy(p_b, tmp, H->key_size);

  free(tmp);
}

void heapify(binheap_type *H, unsigned int node)
{
  //heapify works on 3 nodes: current node and its children
  unsigned int dst_node=node; //node containing the min among node and its children
  unsigned int child;

  do {
    node = dst_node; //update node we are interested in
    //i.e. the minimum identified the previous step

    //identify minimum among current node and its children
    child = RIGHT_CHILD(node); //index of the right child of node
    //need to test if child is a valid node
    //indeed node can have no child
    if(VALID_NODE(H, child) &&
       H->leq(ADDR(H, child), ADDR(H, node))) {
         //this child is a good candidate to be placed as root
         dst_node = child;
       }

    child = LEFT_CHILD(node); //not valid if it is a leaf, or no child
    //need to repeat this operation also for the left child
    if(VALID_NODE(H, child) &&
       H->leq(ADDR(H, child), ADDR(H, node))) {
         //this child is a good candidate to be placed as root
         dst_node = child;
       }

    if(dst_node!=node)
    {
      swap_keys(H, dst_node, node);
    }
  } while(dst_node != node);
}

const void *extract_min(binheap_type *H)
 {
    if(is_heap_empty(H)) {return NULL;}
    //if the heap is not empty, swapping the key among the root
    //and the right-most leaf of the last level
    //A[0] <-> A[num_of_elem-1]

    //out as index of the root in the pos_key the actual index of pos_key[num_of_elem-1]
    //don't want to modify anything in A
    //*((H)->pos_key) = *(((H)->pos_key)+(((H)->num_of_elem-1)*sizeof(int)));
    swap_keys(H, 0, H->num_of_elem-1);

    //deleting the right-most leaf of the last level
    //what I don't care know  is pos_key[num_of_elem]
    //A[num_elem-1]
    H->num_of_elem--; //should be fine

    heapify(H, 0);

    return ADDR(H, H->num_of_elem+1); //min before these operations
}

const void *find_the_max(void *A,
                         const unsigned int num_of_elem,
                         const size_t key_size,
                         total_order_type leq)
{
  if(num_of_elem==0) {return NULL;}
  //assume max is in the first cell of the array
  const void *max_value = A; //address of the first value of A
  //for all the values in A,
  for(const void *addr=A+key_size; addr!=A+num_of_elem*key_size; addr+=key_size)
  {
    //if addr>max_value
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
  H->num_of_elem = num_of_elem;
  H->max_size = max_size;
  H->key_size = key_size;
  H->leq = leq;
  H->max_order_value = malloc(key_size); //void*, no type casting of malloc

  if(num_of_elem==0)
  {
    return H; //nothing to fix
  }

  //get the maximum among A[0:num_of_elem-1]
  //and store it in max_order_value
  const void* value = find_the_max(A, num_of_elem, key_size, leq);
  memcpy(H->max_order_value, value, key_size);

  //fix the heap property from the second last level, up to the root
  for(unsigned int  i=num_of_elem/2; i>0; i--)
  {
    heapify(H, i);
  }
  heapify(H, 0); //not inlcuded in the loop
  return H;
}

void delete_heap(binheap_type *H)
{
    free(H->max_order_value);
    free(H);
}

const void *decrease_key(binheap_type *H, void *node, const void *value)
{
    unsigned int node_idx = INDEX_OF(H, node);
    //we wanto to decrease the key of node of value

    //if node does not belong to H, or value>=node we return null
    if(!(VALID_NODE(H, node_idx) || !(H->leq(value, node))))
    {
      return NULL;
    }

    memcpy(node, value, H->key_size);
    //we may end up in a situation in which heap property doesn't hold anymore
    unsigned int parent_idx = PARENT(node_idx);
    void *parent = ADDR(H, parent_idx);

    while((node_idx!=0)&&(!H->leq(parent, node)))
    {
    //if the parent has value which is greater than that of the node
    //swap parent and node keys
    swap_keys(H, parent_idx, node_idx);
    //moving up the check --> focus on the node's parent
    node = parent;
    node_idx = parent_idx;

    parent_idx = PARENT(node_idx);
    parent = ADDR(H, parent_idx);
    }

   return node;
}

const void *insert_value(binheap_type *H, const void *value)
{ //if the heap is already full
  if(H->max_size == H->num_of_elem)
  {
    //cannot change the size of the array
    return NULL; //only possible behaviour for our code
  }
  //if the value is greater than max_value
  if(H->num_of_elem == 0 || !H->leq(value, H->max_order_value))
  {
    memcpy(H->max_order_value, value, H->key_size);
    //now for sure the max_order_value is the maximum value considered in our computation
  }
  //address in A of the key of the new node
  //we insert the new node as the leftmost available node in the last level
  //the position in the array of this node is the last one
  void *new_node_addr = ADDR(H, H->num_of_elem); //get position of the new node
  memcpy(new_node_addr, H->max_order_value, H->key_size);

  //increase size of the heap by one
  H->num_of_elem++;
  //decrease the key of this new node to the desired value
  return decrease_key(H, new_node_addr, value);
}

void print_heap(const binheap_type *H,
                void (*key_printer)(const void *value))
{
    //we want to print the heap by level
    //to store the index of the leftmost node of the next level
    unsigned int next_level_node = 1;
    for(unsigned int node = 0; node < H->num_of_elem; node++)
    {
      if(node == next_level_node)
      {
        printf("\n"); //new line
        //compute leftmost node of the next level as the left child
        next_level_node = LEFT_CHILD(node);
      }
      else
      {
        printf("\t");
      }

      key_printer(ADDR(H, node));
    }
    printf("\n");
}
