#!/usr/bin/env python3

from compare import *
import matplotlib.pyplot as plt
import seaborn as sns
import os

PENALTY_SCALE = 2.5
LAST_PENALTY_TYPE = f"manual"

def plot():
    plt.rcParams["font.family"] = "serif"

    results = load_all_results(RESULTS_FILE_NAME)

    fig, axs = plt.subplots(2, 2)
    fig.set_figheight(12)
    fig.set_figwidth(18)
    fig.suptitle(f"Optimal run penalty from penalty: {LAST_PENALTY_TYPE}")

    alg = ALGORITHMS[1]
    for complements in (False, True):
        for inp in load_all_inputs(INPUT_FILE_NAME):
            manual_penalty = int(inp.split()[2]) if len(inp.split()) > 2 else 0
            inp = inp.split()[0]

            ys, xs, yys = [], [], []
            for k in KS:
                if not (alg, inp, k, complements, manual_penalty) in results.keys(): continue
                xs.append(k)

                d = results[(alg, inp, k, complements, manual_penalty)]
                ys.append(
                    3 + math.log2(d[LABELS[0]] / d[LABELS[1]])
                )
                yys.append(
                    (manual_penalty if manual_penalty != 0 else math.log2(d[LABELS[0]])) / (3 + math.log2(d[LABELS[0]] / d[LABELS[1]]))
                )

            sns.lineplot(y=ys, x=xs, ax=axs[0][1 if complements else 0],
                        marker="o" if complements else "s", markeredgewidth=0,
                        label=inp[:25] + ": " + str(manual_penalty))
            sns.lineplot(y=yys, x=xs, ax=axs[1][1 if complements else 0],
                        marker="o" if complements else "s", markeredgewidth=0,
                        label=inp[:25])
    for i in range(4): axs[i // 2][i % 2].set(xticks=KS)
        
    plt.legend()
    plt.show()

    file_name = f"./{FIG_DIR}/run_penalty_{LAST_PENALTY_TYPE}.svg"
    if os.path.isfile(file_name):
        a = input("\n===================\nOverwrite old file? (y/...) ")
        if a.strip().lower() == "y": fig.savefig(file_name)
    else:
        with open(file_name, "w") as file: pass
        fig.savefig(file_name)

if __name__ == "__main__":
    plot()
