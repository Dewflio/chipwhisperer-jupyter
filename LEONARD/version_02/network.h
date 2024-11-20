#include <stdio.h>
#include <stdlib.h>
#include <time.h>

typedef struct neuron_struct {
    int num_weights;
    float *weights;
    float bias;
    float z;
    float a;

    int *mul_indices;
} neuron;

typedef struct layer_struct {
    int num_neurons;
    neuron *neurons;
} layer;

typedef struct network_struct {
    int num_layers;
    layer *layers;
} network;

typedef struct neuron_indices_struct {
    int prev_layer_neuron_num;
    int *prev_layer_neuron_indices; 
} neuron_indices;

typedef struct weight_indices_struct {
    int prev_layer_neuron_num;


} weight_indices;

/** @brief

    Description of the data structure:

    neuron_indices - array of arrays of int
        for each layer in the neural network except the first one holds one array
        each of them describes the order in which multiplications with the previous layer are executed
        the size of each of these arrays equal to the number of neurons in the previous layer
    
    weight_indices - array of arrays of arrays of int
        for each layer
        for each neuron in that layer
        for each neuron in previous layer

*/
typedef struct random_indices_struct {
    neuron_indices *neuron_indices;
    int ***weight_indices;
} random_indices;


void print_network(network net);
void free_network(network *net);


neuron create_neuron(void* weights, int num_out_weights, int layer_idx, int neuron_idx); //legacy - neuron create_neuron(int num_out_weights);
layer create_layer(int num_neurons);
network create_network(int num_layers);
network init_network(int num_layers, int *num_neurons, void* weights);  //legacy - network construct_network(int num_outputs, int num_layers, int *num_neurons);

network forward(network net);   //legacy - void forward(network net);
network forward_shuffled(network net);

//void forward_shuffled(network net);
// No Overhead
//void forward_shuffled_NO(network net, int**** random_indices);
// No Overhead, Activations At the End
//void forward_shuffled_NO_AAE(network net, int**** random_indices);
// No Overhead, Activations At the End, Random Dummy Operations
//void forward_shuffled_NO_AAE_RDO(network net, int**** random_indices, int ***random_dummy_operations);


//RANDOM SHUFFLING
void swap(int *a, int *b);
void shuffleArray(int arr[], int size);
int* get_random_indices(int size);
//int ****generate_random_indices(network net);




