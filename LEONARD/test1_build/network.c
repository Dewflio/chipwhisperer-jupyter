#include <network.h>
#include <stdint.h>
#include <math.h>

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

void forward(network net){
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
                ((net.layers[i-1].neurons[k].weights[j]) * (net.layers[i-1].neurons[j].a));
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

        for (volatile int test = 0; test<40; test++) {
            result = scmd *scmd;
            result = scmd *scmd;
            result = scmd *scmd;
            result = scmd *scmd;
        }
    }
}