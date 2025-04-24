#!/usr/bin/env python3

import subprocess
import math
# from threading import Thread
from concurrent.futures import ThreadPoolExecutor, wait
from compare import *

KS = [23, 31]
RESULTS_FILE_NAME = "results_penalty.txt"

MAX_WORKERS = 12

def run_command(command):
    print(*command)
    subprocess.run(command, text=True, check=True)

def compute_missing(limits=None):
    results = load_all_results(RESULTS_FILE_NAME)
    
    thread_inputs = []
    for inp in load_all_inputs(INPUT_FILE_NAME):
        if limits is not None: thread_inputs = []
        inp = inp.split()[0]
        if "subsampled" in inp: continue
        print(inp)
        for alg in ALGORITHMS:
            for k in KS:
                RUN_PENALTIES = list(range(KS[0]))
                for run_penalty in RUN_PENALTIES:
                    rp = run_penalty if alg != ALG_OLD else None
                    for complements in (True, ):
                        if (alg, inp, k, complements, rp) in results.keys():
                            continue

                        args = (inp, alg, k, complements, rp, RESULTS_FILE_NAME)
                        if MULITHREADING:
                            thread_inputs.append(args)
                            # print(args)
                        else: run_with_parameters(*args)
                    if alg == ALG_OLD: break # only count the same thing once
        
        if len(thread_inputs) == 0: continue
        if limits is not None:
            if limits[0] == 0: continue
            print(len(thread_inputs))
            with ThreadPoolExecutor(max_workers=limits[0]) as exe:
                futures = [exe.submit(run_with_parameters, *i) for i in thread_inputs]
                wait(futures)
            
            limits.pop(0)
            if len(limits) == 0:
                limits = None
                thread_inputs = []

    if MULITHREADING and limits is None:
        print(len(thread_inputs))
        with ThreadPoolExecutor(max_workers=MAX_WORKERS) as exe:
            futures = [exe.submit(run_with_parameters, *i) for i in reversed(thread_inputs)]
            wait(futures)
    print("Done")

if __name__ == "__main__":
    limits = None
    if len(sys.argv) > 1:
        limits = list(map(int, sys.argv[1:]))
    compute_missing(limits)
