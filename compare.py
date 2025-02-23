#!/usr/bin/env python3

import subprocess
import matplotlib.pyplot as plt
import math

INPUT_DIR = "data"
LABELS = ["length", "runs", "time", "memory", "objective", ""]
KS = [15, 23, 31, 47, 63, 95, 127]

def get_results_file_name(file_name):
    return file_name.split(".txt")[0] + "_results.txt"

def compute_objective(result, k):
    return result[0] + result[1] * math.log2(result[0] * k)

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

def plot(file_name):
    with open(get_results_file_name(file_name)) as file:
        lines = list(map(lambda x: x.strip(), file.readlines()))
    data = dict()
    for line in lines:
        key, value = line.split(":")
        input_name, k, complements = key.split()
        k = int(k)
        complements = (complements == "C")

        if not input_name in data.keys(): data[input_name] = dict()
        data[input_name][k] = value

    for input_name in sorted(data.keys()):
        d = data[input_name]
        results = [list() for i in range(4)]
        for k in d.keys():
            value = d[k]
            for i, part in enumerate(value.split(";")):
                if len(part.strip()) == 0: continue
        
                results[i][-1].append(list(map(lambda x: int(float(x)), part[0].split())))
                results[i][-1].append(compute_objective(results[i][-1], k))
                results[i][-1].append(0) # TODO use

        fig, ax = plt.subplots(2, 3)
        fig.set_figwidth(12)
        fig.set_figheight(8)
        for i, label in enumerate(LABELS):
            ax[i // 3, i % 3].set_title(label)
            ax[i // 3, i % 3].set_ylim(ymin=0)
            for r in results[i % 2::2]:
                ax[i // 3, i % 3].plot(KS, list(map(lambda x: x[i], r)))
        plt.show()

if __name__ == "__main__":
    compute("compare_inputs.txt", False)
