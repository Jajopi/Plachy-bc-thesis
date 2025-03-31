#!/usr/bin/env python3

import subprocess
import math
# from threading import Thread
from concurrent.futures import ThreadPoolExecutor, wait
from compare import *

RESULTS_FILE_NAME = "results_penalties.txt"

KS = [31]
RUN_PENALTIES = [0, 1, 5, 10, 15, 20, 25, 30]

def run_command(command):
    print(*command)
    subprocess.run(command, text=True, check=True)

def compute_missing():
    results = load_all_results(RESULTS_FILE_NAME)
    
    thread_inputs = []
    for inp in load_all_inputs(INPUT_FILE_NAME)[:-3]:
        print(inp)
        for run_penalty in RUN_PENALTIES:
            for alg in ALGORITHMS:
                k = KS[0]
                for complements in (False, True):
                    if (alg, inp, k, complements, run_penalty) in results.keys():
                        continue

                    if MULITHREADING:
                        thread_inputs.append((inp, alg, k, complements, run_penalty))
                    else: run_with_parameters(inp, alg, k, complements, run_penalty)
    if MULITHREADING:
        print(len(thread_inputs))
        with ThreadPoolExecutor(max_workers=MAX_WORKERS) as exe:
            futures = [exe.submit(run_with_parameters, *i) for i in thread_inputs]
            wait(futures)
    print("Done")

if __name__ == "__main__":
    compute_missing()
