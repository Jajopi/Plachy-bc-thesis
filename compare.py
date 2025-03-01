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

def plot(file_name):
    with open(get_results_file_name(file_name)) as file:
        lines = list(map(lambda x: x.strip(), file.readlines()))
    data = dict()
    for line in lines:
        if len(line.strip()) == 0: continue
        if not ":" in line or not ";" in line: continue
        key, value = line.split(":")
        input_name, k, complements = key.split()
        k = int(k)
        complements = (complements == "C")

        if not input_name in data.keys(): data[input_name] = dict()
        data[input_name][k] = value

    for input_name in sorted(data.keys()):
        d = data[input_name]
        results = [list() for i in range(4)]
        for k in sorted(d.keys()):
            value = d[k]
            for i, part in enumerate(value.split(";")):
                if len(part.strip()) == 0: continue
        
                results[i].append(list(map(lambda x: int(float(x)), part.split())))
                results[i][-1].append(compute_objective(results[i][-1]))
                results[i][-1].append(0)
                results[i][-1][-1] = results[i][-1][-2] / results[0][-1][-2]

        fig, ax = plt.subplots(2, 3)
        fig.set_figwidth(12)
        fig.set_figheight(8)
        fig.suptitle(input_name)
        for i, label in enumerate(LABELS):
            # y_max = 0
            for r in results:
                if len(r) == 0: continue
                y_data = list(map(lambda x: x[i], r))
                # y_max = max(y_max, max(y_data))
                ax[i // 3, i % 3].plot(KS[:len(r)], y_data, )
                for x, y in zip(KS, y_data):
                    ax[i // 3, i % 3].annotate("%s" % y, xy=(x, y), textcoords="data")
            ax[i // 3, i % 3].set_title(label)
            # ax[i // 3, i % 3].set_ylim(ymin=0, ymax=y_max * 1.25)
        plt.show()

if __name__ == "__main__":
    # compute(INPUT_FILE_NAME, False)
    plot(INPUT_FILE_NAME)
