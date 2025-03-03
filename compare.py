#!/usr/bin/env python3

import subprocess
import matplotlib.pyplot as plt
import math

INPUT_DIR = "data"
INPUT_FILE_NAME = "compare_inputs.txt"
LABELS = ["length", "runs", "time", "memory", "objective", "relative obj."]
KS = [15, 23, 31, 47, 63, 95, 127]

def get_results_file_name(file_name):
    return file_name.split(".txt")[0] + "_results.txt"

def compute_objective(result):
    return result[0] + result[1] * math.log2(result[0])

def run_command(command, output_file_name):
    with open(output_file_name, "a") as output_file:
        subprocess.run(command, text=True, check=True, stdout=output_file)

def run_with_parameters(output_name, input_name, k, complements=False):
    print(input_name, k, complements)

    run_command([f"./test_parameters.sh",
                 INPUT_DIR + "/" + input_name,
                 str(k),
                 "F" + ("C" if complements else "")],
                output_name)

def compute(file_name, complements=False):
    with open(file_name) as file:
        lines = list(map(lambda x: x.strip(), file.readlines()))
    for input_name in lines:
        for k in KS:
            run_with_parameters(get_results_file_name(file_name),
                                input_name, k, complements)

if __name__ == "__main__":
    compute(INPUT_FILE_NAME, True)
