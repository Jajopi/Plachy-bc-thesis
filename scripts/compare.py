#!/usr/bin/env python3

import subprocess
import math
# from threading import Thread
from concurrent.futures import ThreadPoolExecutor, wait
import sys

INPUT_DIR = "data"
INPUT_FILE_NAME = "compare_inputs.txt"
RESULTS_FILE_NAME = "results.txt"
LABELS = ["length", "runs", "time", "memory", "objective", "relative"] # mandatory order
LABELS_COMPRESSIONS = ["gzip", "bzip2", "xz", "zpaq"]
FIG_DIR = "figures"

KS = [23, 31, 47, 63, 95, 127] # 15 was excluded
ALG_OLD = "global"
ALGORITHMS = [ALG_OLD, "loac"] # "csac"

MULITHREADING = True
MAX_WORKERS = 6

def run_command(command):
    print(*command)
    subprocess.run(command, text=True, check=True)

def run_with_parameters(input_name, algorithm, k, complements=False, run_penalty=None, result_file=RESULTS_FILE_NAME):
    run_command(["./scripts/measure_run.sh",
                 result_file,
                 "-p", INPUT_DIR + "/" + input_name,
                 "-k", str(k),
                 "-a", algorithm,
                 "-c" if complements else "  ",
                 f"--run-penalty {run_penalty}" if algorithm != ALG_OLD and run_penalty is not None else ""
                ])

def load_all_inputs(file_name):
    with open(file_name) as file:
        return list(filter(lambda x: len(x) > 0 and x[0] != "#",
                           map(lambda x: x.strip(),
                               file.readlines())))

def parse_header(header):
        params = header.split()[1:]
        
        c = False
        alg = ALG_OLD
        run_p = None
        for i, p in enumerate(params):
            if p == "-a":
                alg = params[i + 1]
            elif p == "-p":
                inp = params[i + 1].split("/")[-1]
            elif p == "-k":
                k = int(params[i + 1])
            elif p == "-c":
                c = True
            elif p == "--run-penalty":
                run_p = int(params[i + 1])

        return (alg, inp, k, c, run_p)

def compute_objective(length, runs):
    return 2 * int(length + runs * (1 + math.log2(length / runs)))

def parse_data(input):
    data = dict()
    for d, l in zip(map(lambda x: int(x.split('.')[0]), (input + " 1 1 1 1").split()),
                    LABELS[:4] + LABELS_COMPRESSIONS):
        data[l] = d
    data[LABELS[4]] = compute_objective(data[LABELS[0]], data[LABELS[1]])
    return data

def load_all_results(file_name):
    try:
        with open(file_name) as file:
            data = list(filter(lambda x: len(x) > 0 and x[0] != "#",
                           map(lambda x: x.strip(),
                               file.readlines())))
    except:
        with open(file_name, "w") as file:
            return dict()
    
    results = dict()
    for result in data:
        header, values = result.split(":=")
        key = parse_header(header)
        results[key] = parse_data(values)

    for result in data: # compute relative objective function results and relative compressions
        header, values = result.split(":=")
        key = parse_header(header)
        alg, inp, k, c, rp = key

        results[key][LABELS[5]] = 1
        if alg != ALG_OLD:
            if (ALG_OLD, inp, k, c, None) in results.keys():
                results[key][LABELS[5]] = results[key][LABELS[4]] / results[(ALG_OLD, inp, k, c, None)][LABELS[4]]
            else: results[key][LABELS[5]] = 0
        
        for label in LABELS_COMPRESSIONS:
            if alg != ALG_OLD:
                if (ALG_OLD, inp, k, c, None) in results.keys():
                    results[key][label] = results[key][label] / results[(ALG_OLD, inp, k, c, None)][label]
                else: results[key][label] = 1
        
    # for result in data: # set realtive compressions of ALG_OLD to 1
    #     header, values = result.split(":=")
    #     key = parse_header(header)
    #     alg, inp, k, c, rp = key
        
    #     for label in LABELS_COMPRESSIONS:
    #         if alg == ALG_OLD:
    #             results[(ALG_OLD, inp, k, c, None)][label] = 1

    return results

def compute_missing(limits=None):
    results = load_all_results(RESULTS_FILE_NAME)
    
    thread_inputs = []
    for inp in load_all_inputs(INPUT_FILE_NAME):
        print(inp)
        if limits is not None: thread_inputs = []
        limit = int(inp.split()[1]) if len(inp.split()) > 1 else 128
        run_penalty = int(inp.split()[2]) if len(inp.split()) > 2 else None
        inp = inp.split()[0]
        for alg in (reversed(ALGORITHMS) if MULITHREADING else ALGORITHMS):
            rp = run_penalty if alg != ALG_OLD else None
            for k in KS:
                if k >= limit: continue
                for complements in (True, ):
                    if (alg, inp, k, complements, rp) in results.keys():
                        continue

                    if MULITHREADING:
                        thread_inputs.append((inp, alg, k, complements, rp))
                    else: run_with_parameters(inp, alg, k, complements, rp)
        
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
