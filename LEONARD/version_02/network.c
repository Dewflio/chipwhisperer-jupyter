#include "network.h"
#include <stdint.h>
#include <math.h>

void swap(int *a, int *b){
    int temp = *a;
    *a = *b;
    *b = temp;
}

void free_network(network *net){
    // free the dynamically allocated fields inside the network struct
    for (int i=0; i < net->num_layers; i++){
        
        for(int j=0; j< net->layers[i].num_neurons; j++){
            free(net->layers[i].neurons[j].weights);
        }
        free(net->layers[i].neurons);
    }
    free(net->layers);
}

/** 
* @brief
* Shuffles the array using the Fisher-Yates shuffle
*/
void shuffleArray(int arr[], int size){
    srand(time(NULL));
    for (int i = size - 1; i > 0; i--) {
        int j = rand() % (i + 1);
        swap(&arr[i], &arr[j]);
    }
}

int* get_random_indices(int size){
    int *arr = malloc(size * sizeof(int));
    for (int i=0; i<size; i++){
        arr[i] = i;
    }
    shuffleArray(arr, size);
    return arr;
}

/*
* Prints out the values contained withing the network - number of layers, and then for each layer prints out each neuron
* for each neuron this prints out the a value, the z value, and the weight values (an array of size equal to the number of neurons in the previous layer - except for the first layer - the input layer ) 
*/
void print_network(network net){
    printf("\n");
    printf("-----------------------------------------------------------------------------------------------------------------------------------------------------------------\n");
    printf("Network - num_layers = %d\n", net.num_layers);
    for (int i = 0; i < net.num_layers; i++){
        printf("Layer %d:\n", i);
        for (int j = 0; j < net.layers[i].num_neurons; j++){
            printf("\tNeuron %d | a=%f z=%f\t| ", j, net.layers[i].neurons[j].a,  net.layers[i].neurons[j].z );
            
            if (i >= 1){
                printf("Weights: ");
                for (int k = 0; k < net.layers[i - 1].num_neurons; k++){
                    printf("w%d=%f\t", k, net.layers[i].neurons[j].weights[k]);
                }
            }
            printf("\n");
        }
    }
    printf("-----------------------------------------------------------------------------------------------------------------------------------------------------------------\n");
}


neuron create_neuron(void* weights, int num_in_weights, int layer_idx, int neuron_idx){
    neuron new_neuron;
    new_neuron.a = 0.5;
    new_neuron.z = 0.5;
    new_neuron.bias = 0.0;
    new_neuron.weights = (float*) malloc(num_in_weights * sizeof(float));
    new_neuron.num_weights = num_in_weights;
    
    if (weights != NULL){
        //TODO: Dont question it... it works.
        float (*layer_weights)[num_in_weights] = ((float (*)[num_in_weights])((float**)weights)[layer_idx]);

        for (int i=0; i<num_in_weights; i++){
            new_neuron.weights[i] = layer_weights[neuron_idx][i];
        }
    }
    else {
        for (int i=0; i<num_in_weights; i++){
            new_neuron.weights[i] = ((float)rand())/((float)RAND_MAX);
        }
    }
    
    return new_neuron;
}

layer create_layer(int num_neurons){
    layer lay;
    lay.num_neurons = num_neurons;
    lay.neurons = (neuron*) malloc(num_neurons * sizeof(neuron));
    return lay;
}

network create_network(int num_layers){
    network net;
    net.num_layers = num_layers;
    net.layers = (layer*) malloc(num_layers * sizeof(layer));
    return net;
}

network init_network(int num_layers, int *num_neurons, void* weights) {

    network net = create_network(num_layers);
    int curr_layer_idx, curr_neuron_idx;
    for (curr_layer_idx = 0; curr_layer_idx < num_layers; curr_layer_idx++){
        net.layers[ curr_layer_idx ] = create_layer(num_neurons[ curr_layer_idx ]);
    }

    // create neurons for the first (input) layer - they dont have weights
    for (curr_neuron_idx = 0; curr_neuron_idx < net.layers[0].num_neurons; curr_neuron_idx++){
        net.layers[0].neurons[ curr_neuron_idx ] = create_neuron(NULL, 1, 0, curr_layer_idx);
    }

    // For each following layer create neurons with number of weights eqaual to the number of neurons in the previous layer
    for (curr_layer_idx = 1; curr_layer_idx < num_layers; curr_layer_idx++){
        int prev_layer_idx = curr_layer_idx - 1;

        // pointer to this layers weights

        for (curr_neuron_idx = 0; curr_neuron_idx <net.layers[ curr_layer_idx ].num_neurons; curr_neuron_idx++){
            net.layers[ curr_layer_idx ].neurons[ curr_neuron_idx ] = 
                //create_neuron(net.layers[ prev_layer_idx ].num_neurons);
                create_neuron(weights, net.layers[ prev_layer_idx ].num_neurons, curr_layer_idx, curr_neuron_idx );
        }
    }
    return net;
}



int ****generate_random_indices(network net) {
    int i, j; //, k;

    int **rand_n_indices = malloc((net.num_layers - 1) * sizeof(int*));
    int ***rand_ws_indices  = malloc((net.num_layers - 1) * sizeof(int**));
    
    for (i=1; i<net.num_layers; i++){
        // for each layer
        int *rand_n_idx = get_random_indices(net.layers[i].num_neurons);
        rand_n_indices[i - 1] = rand_n_idx;


        int **rand_w_indices = malloc((net.layers[i].num_neurons) * sizeof(int*));
        // for each neuron in this layer
        for (j=0; j<net.layers[i].num_neurons; j++){
            int *rand_w_idx = get_random_indices(net.layers[i - 1].num_neurons);
            rand_w_indices[j] = rand_w_idx;
        }
        rand_ws_indices[i - 1] = rand_w_indices;
    }

    int ***rand_n_indices_ptr = &rand_n_indices;

    int ****returned_ptr = malloc(2 * sizeof(int***));
    returned_ptr[0] = rand_n_indices_ptr;
    returned_ptr[1] = rand_ws_indices;

    return returned_ptr;
}


void forward_shuffled_NO(network net, int**** random_indices) {
    //int i, j, k,
    int nidx = 0;
    int *rand_n_idx, *rand_w_idx;

    int **rand_n_indices = *random_indices[0];
    int ***rand_ws_indices = random_indices[1];

    //uint8_t result, scmd = 16;
    // for each layer
    for (volatile int i=1; i<net.num_layers; i++){
        

        rand_n_idx = rand_n_indices[i - 1];
        // for each neuron in this layer
        for (volatile int j=0; j<net.layers[i].num_neurons; j++){
            nidx = rand_n_idx[j];  
            //printf("mmmkay\n");
            net.layers[i].neurons[nidx].z = net.layers[i].neurons[nidx].bias;
            //printf("mmmkay\n");

            rand_w_idx = rand_ws_indices[i - 1][j];
            // for all neurons on the previous layer
            
            for (volatile int k=0; k<net.layers[i - 1].num_neurons; k++){
                
                net.layers[i].neurons[nidx].z = net.layers[i].neurons[nidx].z +
                ((net.layers[i-1].neurons[rand_w_idx[k]].weights[nidx]) * (net.layers[i-1].neurons[rand_w_idx[k]].a));
                // We are looking for THIS MULTIPLICATION
            }
            //get a values
            net.layers[i].neurons[nidx].a = net.layers[i].neurons[nidx].z;
            //apply relu
            if(i < net.num_layers-1){
                if((net.layers[i].neurons[nidx].z) < 0)
                {
                    net.layers[i].neurons[nidx].a = 0;
                }

                else
                {
                    net.layers[i].neurons[nidx].a = net.layers[i].neurons[nidx].z;
                }
            }
            //apply sigmoid to the last layer
            else{
                net.layers[i].neurons[nidx].a = 1/(1+exp(-net.layers[i].neurons[nidx].z));
            }
        }
    }
}

void forward_shuffled_NO_AAE(network net, int**** random_indices) {
    //int i, j, k,
    int nidx = 0;
    int *rand_n_idx, *rand_w_idx;

    int **rand_n_indices = *random_indices[0];
    int ***rand_ws_indices = random_indices[1];

    //uint8_t result, scmd = 16;
    // for each layer
    for (volatile int i=1; i<net.num_layers; i++){
        
        rand_n_idx = rand_n_indices[i - 1];
        // for each neuron in this layer
        for (volatile int j=0; j<net.layers[i].num_neurons; j++){
            nidx = rand_n_idx[j];  
            net.layers[i].neurons[nidx].z = net.layers[i].neurons[nidx].bias;


            rand_w_idx = rand_ws_indices[i - 1][j];
            // for all neurons on the previous layer
            for (volatile int k=0; k<net.layers[i - 1].num_neurons; k++){
                
                net.layers[i].neurons[nidx].z = net.layers[i].neurons[nidx].z +
                ((net.layers[i-1].neurons[rand_w_idx[k]].weights[nidx]) * (net.layers[i-1].neurons[rand_w_idx[k]].a));
                // We are looking for THIS MULTIPLICATION
            }
            //get a values
            net.layers[i].neurons[nidx].a = net.layers[i].neurons[nidx].z;
        }

        for (volatile int j=0; j<net.layers[i].num_neurons; j++) {
            //apply relu
            if(i < net.num_layers-1){
                if((net.layers[i].neurons[nidx].z) < 0)
                {
                    net.layers[i].neurons[nidx].a = 0;
                }

                else
                {
                    net.layers[i].neurons[nidx].a = net.layers[i].neurons[nidx].z;
                }
            }
            //apply sigmoid to the last layer
            else{
                net.layers[i].neurons[nidx].a = 1/(1+exp(-net.layers[i].neurons[nidx].z));
            }
        }
    }
}

void forward_shuffled(network net) {
    int i, j, k, nidx;
    int *rand_n_idx, *rand_w_idx;
    //uint8_t result, scmd = 16;
    // for each layer
    for (i=1; i<net.num_layers; i++){
        
        rand_n_idx = get_random_indices(net.layers[i].num_neurons);
        // for each neuron in this layer
        for (j=0; j<net.layers[i].num_neurons; j++){
            nidx = rand_n_idx[j];  
            net.layers[i].neurons[nidx].z = net.layers[i].neurons[nidx].bias;


            rand_w_idx = get_random_indices(net.layers[i - 1].num_neurons);
            // for all neurons on the previous layer
            for (k=0; k<net.layers[i - 1].num_neurons; k++){
                
                net.layers[i].neurons[nidx].z = net.layers[i].neurons[nidx].z +
                ((net.layers[i-1].neurons[rand_w_idx[k]].weights[nidx]) * (net.layers[i-1].neurons[rand_w_idx[k]].a));
                // We are looking for THIS MULTIPLICATION
            }
            //get a values
            net.layers[i].neurons[nidx].a = net.layers[i].neurons[nidx].z;
            //apply relu
            if(i < net.num_layers-1){
                if((net.layers[i].neurons[nidx].z) < 0)
                {
                    net.layers[i].neurons[nidx].a = 0;
                }

                else
                {
                    net.layers[i].neurons[nidx].a = net.layers[i].neurons[nidx].z;
                }
            }
            //apply sigmoid to the last layer
            else{
                net.layers[i].neurons[nidx].a = 1/(1+exp(-net.layers[i].neurons[nidx].z));
            }
        }
    }
}

network forward2(network net){
    volatile int curr_layer_idx, curr_neuron_idx, prev_layer_neuron_idx;
    //uint8_t result, scmd = 16;
    // for each layer
    for (curr_layer_idx=1; curr_layer_idx < net.num_layers; curr_layer_idx++){
        
        int prev_layer_idx = curr_layer_idx - 1;
        // for each neuron in this layer
        for (curr_neuron_idx=0; curr_neuron_idx < net.layers[ curr_layer_idx ].num_neurons; curr_neuron_idx++){   
            net.layers[ curr_layer_idx ].neurons[ curr_neuron_idx ].z = net.layers[ curr_layer_idx ].neurons[ curr_neuron_idx ].bias;

            // for all neurons on the previous layer
            for (prev_layer_neuron_idx = 0; prev_layer_neuron_idx <net.layers[ prev_layer_idx ].num_neurons; prev_layer_neuron_idx++){
                net.layers[ curr_layer_idx ].neurons[ curr_neuron_idx ].z =
                    net.layers[ curr_layer_idx ].neurons[ curr_neuron_idx ].z
                    +
                    (
                        (net.layers[ curr_layer_idx ].neurons[ curr_neuron_idx ].weights[ prev_layer_neuron_idx ])
                        *
                        (net.layers[ prev_layer_idx ].neurons[ prev_layer_neuron_idx ].a)
                    );
                // We are looking for THIS MULTIPLICATION
            }
            //get a values
            net.layers[curr_layer_idx].neurons[ curr_neuron_idx ].a = net.layers[curr_layer_idx].neurons[ curr_neuron_idx ].z;
            //apply relu
            if(curr_layer_idx < net.num_layers-1){
                if((net.layers[curr_layer_idx].neurons[ curr_neuron_idx ].z) < 0)
                {
                    net.layers[curr_layer_idx].neurons[ curr_neuron_idx ].a = 0;
                }
                else
                {
                    net.layers[curr_layer_idx].neurons[ curr_neuron_idx ].a = net.layers[curr_layer_idx].neurons[ curr_neuron_idx ].z;
                }
            }
            //apply sigmoid to the last layer
            else{
                net.layers[curr_layer_idx].neurons[ curr_neuron_idx ].a = 1/(1+exp(-net.layers[curr_layer_idx].neurons[ curr_neuron_idx ].z));
            }
        }
    }
    return net;
}

void forward_shuffled_NO_AAE_RDO(network net, int**** random_indices, int ***random_dummy_operations) {
    int i, j, k, nidx = 0;
    int *rand_n_idx, *rand_w_idx;

    int **rand_n_indices = *random_indices[0];
    int ***rand_ws_indices = random_indices[1];

    int *rand_dummy_ops_idx;

    uint8_t result, scmd = 16;
    // for each layer
    //printf("Network - num layers = %d\n", net.num_layers);
    for (volatile int i=1; i<net.num_layers; i++){
        //printf("Neuron %d - num of neurons = %d\n", j, net.layers[i].num_neurons);
        rand_n_idx = rand_n_indices[i - 1];
        
        // for each neuron in this layer
        for (volatile int j=0; j<net.layers[i].num_neurons; j++){
            nidx = rand_n_idx[j];
            //printf("OK\n");
            net.layers[i].neurons[nidx].z = net.layers[i].neurons[nidx].bias;
            //printf("OK2\n");
            
            rand_w_idx = rand_ws_indices[i - 1][j];
            rand_dummy_ops_idx = random_dummy_operations[i - 1][j];

            // for all neurons on the previous layer
            for (volatile int k=0; k<net.layers[i - 1].num_neurons; k++){
                //printf("OK3\n");
                net.layers[i].neurons[nidx].z = net.layers[i].neurons[nidx].z +
                ((net.layers[i-1].neurons[rand_w_idx[k]].weights[nidx]) * (net.layers[i-1].neurons[rand_w_idx[k]].a));
                // We are looking for THIS MULTIPLICATION

                // Insert a dummy operations if there is a 1 at this index, don't if there is a 0
                if (rand() %2 == 1){
                //if (rand_dummy_ops_idx[k] == 1){
                    result = scmd * scmd;
                    result = scmd * scmd;
                }
            }
            //get a values
            net.layers[i].neurons[nidx].a = net.layers[i].neurons[nidx].z;
        }

        for (volatile int j=0; j<net.layers[i].num_neurons; j++) {
            //apply relu
            if(i < net.num_layers-1){
                if((net.layers[i].neurons[nidx].z) < 0)
                {
                    net.layers[i].neurons[nidx].a = 0;
                }

                else
                {
                    net.layers[i].neurons[nidx].a = net.layers[i].neurons[nidx].z;
                }
            }
            //apply sigmoid to the last layer
            else{
                net.layers[i].neurons[nidx].a = 1/(1+exp(-net.layers[i].neurons[nidx].z));
            }
        }
    }
}