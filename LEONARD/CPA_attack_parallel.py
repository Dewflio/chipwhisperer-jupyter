import numpy as np
import chipwhisperer as cw
import cupy as cp
import struct
import time
import math
from datetime import datetime
from decimal import Decimal
from concurrent.futures import ProcessPoolExecutor, as_completed
#import multiprocessing as mp

def float_to_binary_str(f):
    # Pack the float into 4 bytes (32-bit) using IEEE 754 standard
    [packed] = struct.unpack('!I', struct.pack('!f', f))
    # Convert the packed number to a binary string
    return f"{packed:032b}"

def HW_float32(f):
    # Get the binary representation of the 32-bit float
    binary_str = float_to_binary_str(f)
    # Count and return the number of '1' bits
    return binary_str.count('1')

def get_weights_arr():
    weights_arr = []
    w_val = Decimal('-2.0')
    step = Decimal('0.01')
    for i in range (401):
        weights_arr.append(float(w_val))
        w_val += step
    return weights_arr


class CPA_Attack_Obj:
    def __init__(self, trace_waves, inputs, weights):
        self.trace_waves = trace_waves
        self.inputs = inputs
        self.weights = weights

        self.hypothetical_leakages_allw = []
        self.hypothetical_products_allw = []
        self.average_hypothetical_leakages_allw = []

        self.M = len(self.trace_waves)
        self.N = len(self.weights)

        self.trace_waves_gpu = cp.array(self.trace_waves)
        self.inputs_gpu = cp.array(self.inputs)
        self.weights_gpu = cp.array(self.weights)

        # Calculate hypothetical products and leakages
        self.hypothetical_products_gpu = cp.outer(self.weights_gpu, self.inputs_gpu)  # Shape (N, M)
        self.hypothetical_leakages_gpu = cp.vectorize(HW_float32)(self.hypothetical_products_gpu)
        self.average_hypothetical_leakage_gpu = cp.mean(self.hypothetical_leakages_gpu, axis=1)  # Shape (N,)

    def compute_hypotheticals(self):
        self.hypothetical_leakages_allw = []
        self.hypothetical_products_allw = []
        self.average_hypothetical_leakages_allw = []

        M = len(self.inputs)
        N = len(self.weights)
        
        for i in range(N):
            hypothetical_products = []
            hypothetical_leakages = []
            hypothetical_leakage_sum = 0
            for j in range(M):
                hypothetical_product = self.weights[i] * self.inputs[j]
                hypothetical_products.append(hypothetical_product)
                hypothetical_leakage = HW_float32(hypothetical_product)
                hypothetical_leakages.append(hypothetical_leakage)
                hypothetical_leakage_sum += hypothetical_leakage
            average_hypothetical_leakage = hypothetical_leakage_sum / M

            self.hypothetical_leakages_allw.append(hypothetical_leakages)
            self.hypothetical_products_allw.append(hypothetical_products)
            self.average_hypothetical_leakages_allw.append(average_hypothetical_leakage)
    
    def CPA_attack(self, time_sample_index):
        r_abs = []
        n_trace_waves = len(self.trace_waves)
        leakage_sum = 0
        for trace_wave in self.trace_waves:
            leakage_sum += trace_wave[time_sample_index]    
        average_leakage = leakage_sum / n_trace_waves
        M = len(self.inputs)
        N = len(self.weights)
        
        for i in range(N):
            numerator= 0
            denominator1 = 0
            denominator2 = 0
            for j in range(M):
                numerator += (self.hypothetical_leakages_allw[i][j] - self.average_hypothetical_leakages_allw[i]) * (self.trace_waves[j][time_sample_index] - average_leakage)
                denominator1 += math.pow((self.hypothetical_leakages_allw[i][j] - self.average_hypothetical_leakages_allw[i]), 2)
                denominator2 += math.pow((self.trace_waves[j][time_sample_index] - average_leakage), 2)
            denominator1 = math.sqrt(denominator1)
            denominator2 = math.sqrt(denominator2)

            corr_coef = numerator / (denominator1 * denominator2)
            r_abs.append(abs(corr_coef))
        return (r_abs, time_sample_index)

    def CPA_attack2(self, time_sample_index):
        r_abs = []
        n_trace_waves = len(self.trace_waves)
        leakage_sum = 0
        for trace_wave in self.trace_waves:
            leakage_sum += trace_wave[time_sample_index]    
        average_leakage = leakage_sum / n_trace_waves
        M = len(self.inputs)
        N = len(self.weights)
        
        for i in range(N):
            hypothetical_products = []
            hypothetical_leakages = []
            hypothetical_leakage_sum = 0
            for j in range(M):
                hypothetical_product = self.weights[i] * self.inputs[j]
                hypothetical_products.append(hypothetical_product)
                hypothetical_leakage = HW_float32(hypothetical_product)
                hypothetical_leakages.append(hypothetical_leakage)
                hypothetical_leakage_sum += hypothetical_leakage
            average_hypothetical_leakage = hypothetical_leakage_sum / M

            numerator= 0
            denominator1 = 0
            denominator2 = 0
            for j in range(M):
                numerator += (hypothetical_leakages[j] - average_hypothetical_leakage) * (self.trace_waves[j][time_sample_index] - average_leakage)
                denominator1 += math.pow((hypothetical_leakages[j] - average_hypothetical_leakage), 2)
                denominator2 += math.pow((self.trace_waves[j][time_sample_index] - average_leakage), 2)
            denominator1 = math.sqrt(denominator1)
            denominator2 = math.sqrt(denominator2)

            
            corr_coef = numerator / (denominator1 * denominator2)
            r_abs.append(abs(corr_coef))
        return (r_abs, time_sample_index)
    
    def CPA_attack_gpu(self, time_sample_index):
        r_abs = []
        
        average_leakage = cp.mean(self.trace_waves_gpu[:, time_sample_index])

        r_abs = cp.zeros(len(self.weights_gpu))

        for i in range(self.N):
            # Numerator and denominator calculations, vectorized
            num = cp.sum(
                (self.hypothetical_leakages_gpu[i] - self.average_hypothetical_leakage_gpu[i]) * 
                (self.trace_waves_gpu[:, time_sample_index] - average_leakage)
            )
            den1 = cp.sqrt(cp.sum((self.hypothetical_leakages_gpu[i] - self.average_hypothetical_leakage_gpu[i]) ** 2))
            den2 = cp.sqrt(cp.sum((self.trace_waves_gpu[:, time_sample_index] - average_leakage) ** 2))

        # Compute correlation coefficient
        r_abs[i] = cp.abs(num / (den1 * den2))

        # Convert result back to CPU (optional, if you need it on CPU)
        return r_abs.get(), time_sample_index
    
def run_parallel_CPA_attack(cpa_attack_instance, start_time_sample, n_time_samples):
    start = time.time()
    num_completed = 0

    with ProcessPoolExecutor(max_workers=8) as executor:
        futures = {executor.submit(cpa_attack_instance.CPA_attack_gpu, i): i for i in range(start_time_sample, start_time_sample + n_time_samples)}
    
        results = [None] * n_time_samples
        
        for future in as_completed(futures):
            result, i = future.result()
            results[i - start_time_sample] = result
            num_completed += 1
            if num_completed % 10 == 0:
                print(f'completed {num_completed} time samples in {(time.time() - start)} seconds total')

    
    end = time.time()
    print(f'time elapsed:\t{end - start}')

    return results

if __name__ == "__main__":        
    proj = cw.open_project('project-06')
    trace_waves_arr = []
    inputs_arr = []
    for trace in proj.traces:
        trace_waves_arr.append(trace.wave)
        inputs_arr.append(trace.textin)


    weights_arr = []
    w_val = Decimal('-2.0')
    step = Decimal('0.01')
    for i in range (401):
        weights_arr.append(float(w_val))
        w_val += step

    print(f"number of traces:\t{len(trace_waves_arr)}")
    print(f"number of inputs:\t{len(inputs_arr)}")
    print(f"number of weights:\t{len(weights_arr)}")

    r_abs_all_time_samples = []
    start_time_sample = 400
    n_time_samples = 100

    cpaa_obj = CPA_Attack_Obj(trace_waves=trace_waves_arr, inputs=inputs_arr, weights=weights_arr)
    print(f'computing hypotheticals...')
    cpaa_obj.compute_hypotheticals()
    print(f'running the attack...')

    r_abs_all_time_samples = run_parallel_CPA_attack(cpa_attack_instance=cpaa_obj, start_time_sample=start_time_sample, n_time_samples=n_time_samples)

    print(f'number of r absolute arrays (should match num of traces):\t{len(r_abs_all_time_samples)}')


    # Construct correlation graphs for each weight value
    corr_graphs = []
    print(f'constructing a correlation graph for each weight value...')
    for i in range(len(weights_arr)):
        corr_graph = []
        for u in range(n_time_samples):
            corr_graph.append(r_abs_all_time_samples[u][i])
        corr_graphs.append(corr_graph)

    print(f'getting the maxumal correlation coefficient graphs...')
    max_corr = 0
    max_corr_index = 0
    for i in range(len(corr_graphs)):
        if (max(corr_graphs[i]) > max_corr): #and i!=200):
            max_corr = max(corr_graphs[i])
            max_corr_index = i

    print(f'weight value with the highest correlation coefficient value - w={weights_arr[max_corr_index]} (indexed with {max_corr_index}) - correlation coeff val={max_corr}')


    # Saving the calculated data
    now = datetime.now()
    formatted_now = now.strftime("%d-%m-%Y_%H-%M-%S")
    output_folder = "./np_saved_data/"
    savefile_r_abs = output_folder + formatted_now + "_r_abs.npy"
    savefile_corr_graphs = output_folder + formatted_now + "_corr_graphs.npy"

    try:
        print(f'saving r_abs data into {savefile_r_abs}')
        print(f'saving corr_graphs data into {savefile_corr_graphs}')
        np.save(savefile_r_abs, np.array(r_abs_all_time_samples))
        np.save(savefile_corr_graphs, np.array(corr_graphs))
    except Exception as e:
        print(e)