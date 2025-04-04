#!/usr/bin/env python3

import subprocess
import math
# from threading import Thread
from concurrent.futures import ThreadPoolExecutor, wait

INPUT_DIR = "data"
INPUT_FILE_NAME = "compare_inputs.txt"
RESULTS_FILE_NAME = "results.txt"
LABELS = ["length", "runs", "time", "memory", "objective", "relative"] # mandatory order
FIG_DIR = "figures"

KS = [23, 31, 47, 63, 95, 127] # 15 was excluded
ALG_OLD = "global"
ALGORITHMS = [ALG_OLD, "loac"] # "csac"

MULITHREADING = True
MAX_WORKERS = 2

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
    # return int(length + runs * math.log2(length)) # TODO switch from RLE to EF penalty
    return int(length + runs * (3 + math.log2(length / runs)))

def parse_data(input):
    data = dict()
    for d, l in zip(map(lambda x: int(x.split('.')[0]), input.split()),
                    LABELS[:4]):
        data[l] = d
    data[LABELS[4]] = compute_objective(data[LABELS[0]], data[LABELS[1]])
    return data

def load_all_results(file_name):
    try:
        with open(file_name) as file:
            results = list(filter(lambda x: len(x) > 0 and x[0] != "#",
                           map(lambda x: x.strip(),
                               file.readlines())))
    except:
        with open(file_name, "w") as file:
            return dict()
    
    done = dict()
    for result in results:
        header, data = result.split(":=")
        key = parse_header(header)
        done[key] = parse_data(data)

    for result in results: # compute relative objective function results
        header, data = result.split(":=")
        key = parse_header(header)
        alg, inp, k, c, rp = key
        done[key][LABELS[5]] = (done[key][LABELS[4]] / done[(ALG_OLD, inp, k, c, None)][LABELS[4]]
                                if alg != ALG_OLD
                                else 1)

    return done

def compute_missing():
    results = load_all_results(RESULTS_FILE_NAME)
    
    thread_inputs = []
    for inp in load_all_inputs(INPUT_FILE_NAME):
        limit = int(inp.split()[1]) if len(inp.split()) > 1 else 128
        run_penalty = int(inp.split()[2]) if len(inp.split()) > 2 else None
        inp = inp.split()[0]
        for alg in ALGORITHMS:
            rp = run_penalty if alg != ALG_OLD else None
            for k in KS:
                if k >= limit: continue
                for complements in (False, True):
                    if (alg, inp, k, complements, rp) in results.keys():
                        continue

                    if MULITHREADING:
                        thread_inputs.append((inp, alg, k, complements, rp))
                    else: run_with_parameters(inp, alg, k, complements, rp)
    if MULITHREADING:
        print(len(thread_inputs))
        with ThreadPoolExecutor(max_workers=MAX_WORKERS) as exe:
            futures = [exe.submit(run_with_parameters, *i) for i in thread_inputs]
            wait(futures)
    print("Done")

if __name__ == "__main__":
    compute_missing()
