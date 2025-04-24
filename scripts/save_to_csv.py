#!/usr/bin/env python3

from compare import *

CSV_SEPARATOR = ','

def print_to_csv(results_or_resources=True, k=31, complements=True, subsampled=False, pangenomes=False, path_prefix=""):
    results = load_all_results(RESULTS_FILE_NAME)

    if subsampled and k != 31: return

    out_file = path_prefix + "outputs/"
    out_file += "results-length-runs" if results_or_resources else "results-time-memory"
    if not complements: out_file += "-noncomplement"
    if pangenomes: out_file += "-pan"
    if subsampled: out_file += "-sub"
    out_file += "-" + str(k)
    out_file += ".csv"

    with open(out_file, "w") as file:
        if results_or_resources:
            file.write(CSV_SEPARATOR.join(("dataset", "method", "length", "runs", "objective", "relative")) + '\n')
        else:
            file.write(CSV_SEPARATOR.join(("dataset", "method", "time", "relativetime", "memory", "relativememory")) + '\n')
    
        for inp in load_all_inputs(INPUT_FILE_NAME):
            run_penalty = int(inp.split()[2]) if len(inp.split()) > 2 else None
            
            name = inp.split()[3].replace('_', ' ')
            if name == "-": continue
            if not subsampled:
                if ("(p)" in name or "RNA" in name) and not pangenomes: continue
                if not ("(p)" in name or "RNA" in name) and pangenomes: continue
            name = name.replace("(p)", "")

            inp = inp.split()[0]
            if ("subsampled" in inp) != subsampled: continue
            # print(name)
            
            for alg in (reversed(ALGORITHMS)):
                rp = run_penalty if alg != ALG_OLD else None
                if not (alg, inp, k, complements, rp) in results.keys(): continue
                result = results[(alg, inp, k, complements, rp)]
                
                if results_or_resources:
                    file.writelines(map(lambda x: str(x).strip(), (
                        name if alg != ALG_OLD else "", CSV_SEPARATOR,
                        "GGMO" if alg == ALG_OLD else f"LOAC ({rp})", CSV_SEPARATOR,
                        result[LABELS[0]], CSV_SEPARATOR,
                        result[LABELS[1]], CSV_SEPARATOR,
                        result[LABELS[4]], CSV_SEPARATOR,
                        f"{result[LABELS[5]]:3f}" if alg != ALG_OLD else ""
                    )))
                else:
                    relative_time = 1 if alg == ALG_OLD else result[LABELS[2]] / results[(ALG_OLD, inp, k, complements, None)][LABELS[2]]
                    relative_memory = 1 if alg == ALG_OLD else result[LABELS[3]] / results[(ALG_OLD, inp, k, complements, None)][LABELS[3]]
                    file.writelines(map(lambda x: str(x).strip(), (
                        name if alg != ALG_OLD else "", CSV_SEPARATOR,
                        "GGMO" if alg == ALG_OLD else "LOAC", CSV_SEPARATOR,
                        result[LABELS[2]], CSV_SEPARATOR,
                        "" if alg == ALG_OLD else f"{relative_time:3f}", CSV_SEPARATOR,
                        result[LABELS[3]], CSV_SEPARATOR,
                        "" if alg == ALG_OLD else f"{relative_memory:3f}"
                    )))
                file.write('\n')
    print(out_file)

if __name__ == "__main__":
    for k in (23, 31, 63):
        for res_or_res in (True, False):
            print_to_csv(results_or_resources=res_or_res, k=k, pangenomes=True)
            print_to_csv(results_or_resources=res_or_res, k=k, pangenomes=False)
            print_to_csv(results_or_resources=res_or_res, k=k, subsampled=True)
            print_to_csv(results_or_resources=res_or_res, k=k)
