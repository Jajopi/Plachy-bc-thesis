#!/usr/bin/env python3

import subprocess
import matplotlib.pyplot as plt
import seaborn as sns
import math

INPUT_DIR = "data"
INPUT_FILE_NAME = "compare_inputs.txt"
RESULTS_FILE_NAME = "results.txt"
LABELS = ["length", "runs", "time", "memory", "objective", "relative"] # mandatory order
FIG_DIR = "figures"

RECOMPUTE_ALL = False

KS = [15, 23, 31, 47, 63, 95, 127]
ALG_OLD, ALG_NEW = "global", "csac"
ALGORITHMS = [ALG_OLD, ALG_NEW]

def compute_objective(result):
    return result[0] + result[1] * math.log2(result[0])

def run_command(command):
    print(*command)
    subprocess.run(command, text=True, check=True)

def run_with_parameters(input_name, algorithm, k, complements=False):
    run_command(["./scripts/measure_run.sh",
                 "-p", INPUT_DIR + "/" + input_name,
                 "-k", str(k),
                 "-a", algorithm,
                 "--precision 100" if algorithm == ALG_NEW else "",
                 "-c" if complements else ""
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
        for i, p in enumerate(params):
            if p == "-a":
                alg = params[i + 1]
            elif p == "-p":
                inp = params[i + 1].split("/")[-1]
            elif p == "-k":
                k = int(params[i + 1])
            elif p == "-c":
                c = True

        return (alg, inp, k, c)

def compute_objective(length, runs):
        return int(length + runs * math.log2(length))

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
        alg, inp, k, c = key
        done[key][LABELS[5]] = (done[key][LABELS[4]] / done[(ALG_OLD, inp, k, c)][LABELS[4]]
                                if alg == ALG_NEW
                                else 1)

    return done

def compute_missing():
    results = load_all_results(RESULTS_FILE_NAME)

    for inp in load_all_inputs(INPUT_FILE_NAME):
        big = len(inp.split()) > 1 and inp.split()[1] == "-"
        inp = inp.split()[0]
        for alg in ALGORITHMS:
            for k in (KS[:5] if big else KS):
                for complements in (False, True):
                    if (alg, inp, k, complements) in results.keys():
                        if not RECOMPUTE_ALL:
                            continue
                        if alg != ALG_NEW: # no need to recompute gg
                            continue

                    run_with_parameters(inp, alg, k, complements)

def plot_all():
    plt.rcParams["font.family"] = "serif"

    results = load_all_results(RESULTS_FILE_NAME)

    for inp in load_all_inputs(INPUT_FILE_NAME):
        inp = inp.split()[0]
        print(inp)

        fig, axs = plt.subplots(2, 3)
        fig.set_figheight(12)
        fig.set_figwidth(18)
        fig.suptitle(inp)
        for alg in ALGORITHMS:
            for complements in (False, True):
                ys, xs = [list() for label in LABELS], []
                for k in KS[1:]: # exclude k = 15
                    if not (alg, inp, k, complements) in results.keys(): continue

                    d = results[(alg, inp, k, complements)]
                    # print((alg, inp, k, complements))
                    xs.append(k)
                    for i in range(len(LABELS)):
                        ys[i].append(d[LABELS[i]])
        
                for i, label in enumerate(LABELS):
                    sns.lineplot(y=ys[i], x=xs, ax=axs[i // 3, i % 3],
                                 marker="o" if complements else "s", markeredgewidth=0,
                                 label=alg + ("(c)" if complements else ""))
                    # for y, x in zip(ys[i], xs):
                    #     ax.text(x, y, "{:.2E}".format(y).split("E")[0],
                    #             alpha=0.5,
                    #             horizontalalignment="left" if complements else "right",
                    #             verticalalignment="bottom" if alg == ALG_NEW else "top",
                    #             stretch="ultra-condensed",
                    #             rotation=45 * (1 if complements == (alg == ALG_NEW) else -1))
        
        for i, label in enumerate(LABELS):
            ax = axs[i // 3, i % 3]
            ax.set_title(label)
            ax.set(xticks=xs)
        plt.legend()
        # plt.show()
        fig.savefig(f"./{FIG_DIR}/{inp}.svg")

if __name__ == "__main__":
    # compute_missing()
    plot_all()
