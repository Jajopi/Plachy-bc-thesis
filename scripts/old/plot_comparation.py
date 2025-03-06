#!/usr/bin/env python3

from compare import *

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
                y_data = list(map(lambda x: x[i], r))[1:]
                # y_max = max(y_max, max(y_data))
                ax[i // 3, i % 3].plot(KS[1:len(r)], y_data, )
                for x, y in zip(KS[1:], y_data):
                    ax[i // 3, i % 3].annotate("%s" % y, xy=(x, y), textcoords="data")
            ax[i // 3, i % 3].set_title(label)
            # ax[i // 3, i % 3].set_ylim(ymin=0, ymax=y_max * 1.25)
        plt.show()

if __name__ == "__main__":
    plot(INPUT_FILE_NAME)
