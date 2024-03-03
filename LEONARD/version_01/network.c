#include <network.h>
#include <stdint.h>
#include <math.h>

void swap(int *a, int *b){
    int temp = *a;
    *a = *b;
    *b = temp;
}

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

neuron create_neuron(int num_out_weights){
    neuron new_neuron;
    new_neuron.bias = 0.0;
    new_neuron.z = 0.0;
    new_neuron.a = 0.0;
    new_neuron.weights = (float*) malloc(num_out_weights * sizeof(float));
    new_neuron.num_weights = num_out_weights;

    for (int i=0; i<num_out_weights; i++){
        new_neuron.weights[i] = ((float)rand())/((float)RAND_MAX);
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

network construct_network(int num_outputs, int num_layers, int *num_neurons) {

    network net = create_network(num_layers);
    int i, j;
    for (i=0; i<num_layers; i++){
        net.layers[i] = create_layer(num_neurons[i]);
    }

    // For each layer create neurons with number of weights eqaual to the number of neurons in the following layer
    for (i=1; i<num_layers; i++){
        for (j=0; j<net.layers[i - 1].num_neurons; j++){
            net.layers[i - 1].neurons[j] = create_neuron(net.layers[i].num_neurons);
        }
    }
    // Create neurons for the last layer - the number of outputs is given "by the user"
    for (j=0; j<net.layers[num_layers - 1].num_neurons; j++){
        net.layers[num_layers - 1].neurons[j] = create_neuron(num_outputs);
    }
    return net;
}

void forward_layer(network net, int layer_idx){
    
}

int ****generate_random_indices(network net) {
    int i, j, k;

    int **rand_n_indices = malloc((net.num_layers - 1) * sizeof(int*));
    int ***rand_ws_indices  = malloc((net.num_layers - 1) * sizeof(int**));
    
    for (i=1; i<net.num_layers; i++){
        // for each neuron in this layer
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

void forward_shuffled_without_overhead(network net, int**** random_indices, int dummy_operations) {
    int i, j, k, nidx;
    int *rand_n_idx, *rand_w_idx;

    int **rand_n_indices = *random_indices[0];
    int ***rand_ws_indices = random_indices[1];

    uint8_t result, scmd = 16;
    // for each layer
    for (volatile i=1; i<net.num_layers; i++){
        
        rand_n_idx = rand_n_indices[i - 1];
        // for each neuron in this layer
        for (volatile j=0; j<net.layers[i].num_neurons; j++){
            nidx = rand_n_idx[j];  
            net.layers[i].neurons[nidx].z = net.layers[i].neurons[nidx].bias;


            rand_w_idx = rand_ws_indices[i - 1][j];
            // for all neurons on the previous layer
            for (volatile k=0; k<net.layers[i - 1].num_neurons; k++){
                
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

        //OPTIONAL DUMMY OPERATIONS
        if (dummy_operations > 0){
            for (volatile int dummy_i = 0; dummy_i<10; dummy_i++) {
                result = scmd *scmd;
                result = scmd *scmd;
                result = scmd *scmd;
            }
        }   
    }
}

void forward_shuffled_without_overhead_activations_at_end(network net, int**** random_indices, int dummy_operations) {
    int i, j, k, nidx;
    int *rand_n_idx, *rand_w_idx;

    int **rand_n_indices = *random_indices[0];
    int ***rand_ws_indices = random_indices[1];

    uint8_t result, scmd = 16;
    // for each layer
    for (volatile i=1; i<net.num_layers; i++){
        
        rand_n_idx = rand_n_indices[i - 1];
        // for each neuron in this layer
        for (volatile j=0; j<net.layers[i].num_neurons; j++){
            nidx = rand_n_idx[j];  
            net.layers[i].neurons[nidx].z = net.layers[i].neurons[nidx].bias;


            rand_w_idx = rand_ws_indices[i - 1][j];
            // for all neurons on the previous layer
            for (volatile k=0; k<net.layers[i - 1].num_neurons; k++){
                
                net.layers[i].neurons[nidx].z = net.layers[i].neurons[nidx].z +
                ((net.layers[i-1].neurons[rand_w_idx[k]].weights[nidx]) * (net.layers[i-1].neurons[rand_w_idx[k]].a));
                // We are looking for THIS MULTIPLICATION
            }
            //get a values
            net.layers[i].neurons[nidx].a = net.layers[i].neurons[nidx].z;
            if (dummy_operations > 1){
            for (volatile int dummy_i = 0; dummy_i<50; dummy_i++) {
                result = scmd *scmd;
                result = scmd *scmd;
                result = scmd *scmd;
            }
        } 
        }

        for (volatile j=0; j<net.layers[i].num_neurons; j++) {
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
        //OPTIONAL DUMMY OPERATIONS
        if (dummy_operations > 0){
            for (volatile int dummy_i = 0; dummy_i<50; dummy_i++) {
                result = scmd *scmd;
                result = scmd *scmd;
                result = scmd *scmd;
            }
        } 
    }
}

void forward_shuffled(network net, int dummy_operations) {
    int i, j, k, nidx;
    int *rand_n_idx, *rand_w_idx;
    uint8_t result, scmd = 16;
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
        //OPTIONAL DUMMY OPERATIONS
        if (dummy_operations > 0){
            for (volatile int dummy_i = 0; dummy_i<10; dummy_i++) {
                result = scmd *scmd;
                result = scmd *scmd;
                result = scmd *scmd;
            }
        }
    }
}

void forward(network net, int dummy_operations){
    int i, j, k;
    uint8_t result, scmd = 16;
    // for each layer
    for (i=1; i<net.num_layers; i++){
        
        // for each neuron in this layer
        for (j=0; j<net.layers[i].num_neurons; j++){   
            net.layers[i].neurons[j].z = net.layers[i].neurons[j].bias;

            // for all neurons on the previous layer
            for (k=0; k<net.layers[i - 1].num_neurons; k++){
                net.layers[i].neurons[j].z = net.layers[i].neurons[j].z +
                ((net.layers[i-1].neurons[k].weights[j]) * (net.layers[i-1].neurons[k].a));
                // We are looking for THIS MULTIPLICATION
            }
            //get a values
            net.layers[i].neurons[j].a = net.layers[i].neurons[j].z;
            //apply relu
            if(i < net.num_layers-1){
                if((net.layers[i].neurons[j].z) < 0)
                {
                    net.layers[i].neurons[j].a = 0;
                }

                else
                {
                    net.layers[i].neurons[j].a = net.layers[i].neurons[j].z;
                }
            }
            //apply sigmoid to the last layer
            else{
                net.layers[i].neurons[j].a = 1/(1+exp(-net.layers[i].neurons[j].z));
            }
        }
        //OPTIONAL DUMMY OPERATIONS
        if (dummy_operations > 0){
            for (volatile int dummy_i = 0; dummy_i<10; dummy_i++) {
                result = scmd *scmd;
                result = scmd *scmd;
                result = scmd *scmd;
            }
        }
    }
}