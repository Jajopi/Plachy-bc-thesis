#!/usr/bin/env python3

from compare import *
import matplotlib.pyplot as plt
import seaborn as sns

def plot_all():
    plt.rcParams["font.family"] = "serif"

    results = load_all_results(RESULTS_FILE_NAME)

    for inp in load_all_inputs(INPUT_FILE_NAME):
        inp = inp.split()[0]
        print(inp)
        empty = True
        max_xs = []

        fig, axs = plt.subplots(2, 3)
        fig.set_figheight(12)
        fig.set_figwidth(18)
        fig.suptitle(inp)
        for alg in ALGORITHMS:
            for complements in (False, True):
                ys, xs = [list() for label in LABELS], []
                for k in KS[1:]: # exclude k = 15
                    if not (alg, inp, k, complements) in results.keys(): continue
                    empty = False

                    d = results[(alg, inp, k, complements)]
                    # print((alg, inp, k, complements))
                    xs.append(k)
                    for i in range(len(LABELS)):
                        ys[i].append(d[LABELS[i]])

                if len(xs) > len(max_xs): max_xs = xs
        
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
        
        if not empty:
            for i, label in enumerate(LABELS):
                ax = axs[i // 3, i % 3]
                ax.set_title(label)
                ax.set(xticks=max_xs)
            plt.legend()
            # plt.show()
            fig.savefig(f"./{FIG_DIR}/{inp}.svg")
        else: print("(skipped)")

if __name__ == "__main__":
    plot_all()
