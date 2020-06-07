#include <binheap.h>
#include <string.h>
#include <stdio.h>

#define PARENT(node) ((node-1)/2)
#define LEFT_CHILD(node) (2*(node)+1)
#define RIGHT_CHILD(node) (2*(node+1))

#define VALID_NODE(H,node) ((H)->num_of_elem>(node))

#define ADDR(H,node) ((H)->A+node*(H)->key_size)

#define INDEX_OF(H,addr) (((addr)-((H)->A))/(H)->key_size)


int is_heap_empty(const binheap_type *H)
{
    return H->num_of_elem==0;
}


const void *min_value(const binheap_type *H)
{
    if(is_heap_empty(H)) return NULL;
    return ADDR(H,0);
}


void swap_keys(binheap_type *H, unsigned int n_a, unsigned int n_b)
{
    void *p_a=ADDR(H,n_a);
    void *p_b=ADDR(H,n_b);
    void *tmp=malloc(H->key_size);

    memcpy(tmp,p_a,H->key_size);
    memcpy(p_a,p_b,H->key_size);
    memcpy(p_b,tmp,H->key_size);

    free(tmp);
}


void heapify(binheap_type *H,unsigned int node)
{
    unsigned int dst_node=node, child;

    do
    {
        node=dst_node;
        child=RIGHT_CHILD(node);
        if(VALID_NODE(H,child) && H->leq(ADDR(H,child), ADDR(H,dst_node))) dst_node=child;
        child=LEFT_CHILD(node);
        if(VALID_NODE(H,child) && H->leq(ADDR(H,child), ADDR(H,dst_node))) dst_node=child;
        if(dst_node!=node) swap_keys(H,dst_node,node);

    } while(dst_node != node);
}


void *extract_min(binheap_type *H)
{
    if(is_heap_empty(H)) return NULL;

    swap_keys(H,0, H->num_of_elem-1);

    H->num_of_elem--;

    heapify(H,0);
    return ADDR(H,H->num_of_elem);
}


const void *find_the_max(void *A,const unsigned int num_of_elem,const size_t key_size,total_order_type leq)
{
    if(num_of_elem==0) return NULL;

    const void *max_value=A;

    for(const void *addr=A+key_size;addr!=A+num_of_elem*key_size;addr+=key_size)
    {
        if(!leq(addr,max_value)) max_value=addr;
    }

    return max_value;
}


binheap_type *build_heap(void *A,const unsigned int num_of_elem,const unsigned int max_size,const size_t key_size,total_order_type leq)
{
    binheap_type *H=(binheap_type*)malloc(sizeof(binheap_type));

    H->A=A;
    H->num_of_elem=num_of_elem;
    H->max_size=max_size;
    H->key_size=key_size;
    H->leq=leq;
    H->max_order_value=malloc(key_size);

    if(num_of_elem==0) return H;

    const void *value=find_the_max(A,num_of_elem,key_size,leq);
    memcpy(H->max_order_value, value, key_size);

    for(unsigned int i=num_of_elem/2;i>0;i--) heapify(H,i);

    heapify(H,0);

    return H;
}


void delete_heap(binheap_type *H)
{
    free(H->max_order_value);
    free(H);
}


void *decrease_key(binheap_type *H, void *node, const void *value)
{
    unsigned int node_idx=INDEX_OF(H,node);

    if(!VALID_NODE(H,node_idx) || !(H->leq(value,node))) return NULL;

    memcpy(node, value, H->key_size);

    unsigned int parent_idx=PARENT(node_idx);
    void *parent=ADDR(H,parent_idx);

    while((node_idx!=0) && (!H->leq(parent,node)))
    {
        swap_keys(H,parent_idx,node_idx);

        node=parent;
        node_idx=parent_idx;
        parent_idx=PARENT(node_idx);
        parent=ADDR(H,parent_idx);
    }

    return node;

}


const void *insert_value(binheap_type *H, const void *value)
{

    if(H->max_size==H->num_of_elem) return NULL;

    if(H->num_of_elem==0 || !H->leq(value,H->max_order_value)) memcpy(H->max_order_value,value,H->key_size);

    void *new_node_addr=ADDR(H,H->num_of_elem);
    memcpy(new_node_addr,H->max_order_value,H->key_size);

    H->num_of_elem++;

    return decrease_key(H,new_node_addr,value);
}


void print_heap(const binheap_type *H,
                void (*key_printer)(const void *value))
{
    unsigned int next_level_node=1;
    for(unsigned int node=0;node<H->num_of_elem;node++){
        if(node==next_level_node){
            printf("\n");
            next_level_node=LEFT_CHILD(node);
        }else{
            printf("\t");
        }
     key_printer(ADDR(H,node));
    }
    printf("\n");
}
