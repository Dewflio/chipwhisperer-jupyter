/*
 * SimpleSerial V2 Template C code
 * Can be freely used to implement ChipWhisperer target binaries.
 *
 * Date: 14th March 2021
 */

/*
 * When debugging locally compile using `gcc -o debug-app.exe main.c network.c debug-source.c -DDEBUGGING=1`
 */

#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include "network.h"
#include "network_config.h"

#include "hal/hal.h"
#include "hal/stm32f3/stm32f3_hal.h"

#define SS_VER SS_VER_2_1

#include "simpleserial/simpleserial.h"


/// This function will handle the 'p' command send from the capture board.
/// It returns the squared version of the scmd given.
/// It does this in approximately equal time, which allows us to see clear
/// differences between different scmd values.
uint8_t handle(uint8_t cmd, uint8_t scmd, uint8_t len, uint8_t *buf)
{

  int num_layers = NET_NUM_LAYERS;
  int *num_neurons_arr = NET_NUM_NEURONS;
  srand(time(NULL));
  
  network net = init_network(num_layers, num_neurons_arr, net_config_layer_weights);
  
  //Change the input of the first neuron in the first layer to the provided number
  //convert to float
  float input_value;
  uint8_t input_buffer[4] = {buf[0], buf[1], buf[2], buf[3]};
  memcpy(&input_value, input_buffer, sizeof(float)); 
  net.layers[0].neurons[0].a = input_value;

  #ifdef DEBUGGING
  #if 0
  print_network(net);
  #endif
  #endif
  //int ****random_indices = generate_random_indices(net);
  //int ***random_dummy_operations_indices = generate_random_dummy_operations(net);
  //int ***random_dummy_operations_indices = NULL;
  // Start measurement.
  trigger_high();

  net = forward_shuffled(net);

  //forward_shuffled(net);
  //forward_shuffled_NO(net, random_indices);
  //forward_shuffled_NO_AAE(net, random_indices, 2);
  //forward_shuffled_NO_AAE_RDO(net, random_indices, random_dummy_operations_indices);
  // Stop measurement.
  trigger_low();

  #ifdef DEBUGGING
  #if 1
  //print network (a, z, and w values)
  print_network(net);
  #endif
  #endif
  
  //free dynamically allocated memory
  free_network(&net);
  
  simpleserial_put('r', len, buf);

  return 0;
}

int main(void) {
  //Initialize network weights
  init_weights();
  // Setup the specific chipset.
  platform_init();
  // Setup serial communication line.
  init_uart();
  // Setup measurement trigger.
  trigger_setup();

  simpleserial_init();

  // Insert your handlers here.
  simpleserial_addcmd('p', 16, handle);

  // What for the capture board to send commands and handle them.
  while (1)
    simpleserial_get();
}