#!/usr/bin/env python3

import subprocess
import math
# from threading import Thread
from concurrent.futures import ThreadPoolExecutor, wait
from compare import *

KS = [31]
RUN_PENALTIES = list(range(KS[0] + 1))
RESULTS_FILE_NAME = "results_penalty.txt"

MAX_WORKERS = 4

def run_command(command):
    print(*command)
    subprocess.run(command, text=True, check=True)

def compute_missing():
    results = load_all_results(RESULTS_FILE_NAME)
    
    thread_inputs = []
    for inp in load_all_inputs(INPUT_FILE_NAME):
        inp = inp.split()[0]
        print(inp)
        for alg in ALGORITHMS:
            for run_penalty in RUN_PENALTIES:
                rp = run_penalty if alg != ALG_OLD else None
                k = KS[0]
                for complements in (False, True):
                    if (alg, inp, k, complements, rp) in results.keys():
                        continue

                    args = (inp, alg, k, complements, rp, RESULTS_FILE_NAME)
                    if MULITHREADING:
                        thread_inputs.append(args)
                        # print(args)
                    else: run_with_parameters(*args)
                if alg == ALG_OLD: break # only count the same thing once

    if MULITHREADING:
        print(len(thread_inputs))
        with ThreadPoolExecutor(max_workers=MAX_WORKERS) as exe:
            futures = [exe.submit(run_with_parameters, *i) for i in thread_inputs]
            wait(futures)
    print("Done")

if __name__ == "__main__":
    compute_missing()
