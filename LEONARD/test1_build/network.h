#include <stdio.h>
#include <stdlib.h>
#include <time.h>

typedef struct neuron_struct {
    int num_weights;
    float *weights;
    float bias;
    float z;
    float a;
} neuron;

typedef struct layer_struct {
    int num_neurons;
    neuron *neurons;
} layer;

typedef struct network_struct {
    int num_layers;
    layer *layers;
} network;

neuron create_neuron(int num_out_weights);

layer create_layer(int num_neurons);

network create_network(int num_layers);

network construct_network(int num_outputs, int num_layers, int *num_neurons);

void forward_layer(network net, int layer_idx);

void forward_shuffled(network net);

void forward_shuffled_without_overhead(network net, int**** random_indices);

void forward_shuffled_without_overhead_activations_at_end(network net, int**** random_indices);

int ****generate_random_indices(network net);

void forward(network net);

//RANDOM SHUFFLING
void swap(int *a, int *b);

void shuffleArray(int arr[], int size);

int* get_random_indices(int size);

