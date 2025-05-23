# This repo README

## Running and testing scripts

There are quite many testing scripts now in the `scripts` directory.

If you want to compare both `gg` and `loac` with one specific input,
use `test_parameters.sh` with three to four arguments in defined order:
- the method to compare ('loac' or 'csac' -- deprecated)
- the path to input file (usually `data/...`)
- k
- optionally string of one or more uppercase letters:
  - `C` to run also for complements
  - `N` to also count the number of kmers in resulting masked superstrings
  (is written in python and takes awfully lot of time for larger sets, uses `count[_noncomplement]_kmers.py`)
  - `F` to output all four paramters of each of the program runs on a single line
  (in order length of MS, number of runs, time spent, max memory used) --
  this option is priamry intended for legacy results plotting
  - `L` to use local directory instead a `/tmp` one
  - `S` with `L` to read the results previously computed in the local directory

If you want to compare huge amounts inputs, define them in file `compare_inputs.txt`
and use `compare.py`.
The structure of the text file is explained inside, with comments.
Some parameters of the computation, such as number of threads to use or value of `k` to test
can be changed in the python script.
These are then used for plotting in the jupyter notebook (`visualization.ipynb`), too.
Internal shellscript for usage of this script is `measure_run.sh`.

Searching for the optimal run penalties in the file `compare_penalties.py`
works similarly as `compare.py`. Both can be visualized using the jupyter notebook.

`save_to_csv.py` is used to create csv tables to be displayed in the thesis.

`copy_to_kam.sh` is used to copy files to *kamenac* and invoke recompilation
of current version of 🐫.
You probably won't use it.
By default, the set of files to be copied is specified in the script
(just as the target directory).
You can provide any number of arguments, then those files are copied instead.
If you provide at least one argument, the recompilation is not invoked.
If the first argument is `-`, default files are copied without recompilation.

# README from the forked version of [KmerCamel🐫](https://github.com/OndrejSladky/kmercamel)
[![KmerCamel test](https://github.com/OndrejSladky/kmercamel/actions/workflows/ci.yml/badge.svg)](https://github.com/OndrejSladky/kmercamel/actions/)

<!-- vim-markdown-toc GFM -->

* [Introduction](#introduction)
* [Prerequisites](#prerequisites)
* [Getting started](#getting-started)
  * [Installation](#installation)
  * [Compression for k-mer set storage](#compression-for-k-mer-set-storage)
  * [k-mer set indexing](#k-mer-set-indexing)
* [Detailed instructions](#detailed-instructions)
  * [Arguments](#arguments)
  * [Converting k-mer set superstring representation to the (r)SPSS representations](#converting-k-mer-set-superstring-representation-to-the-rspss-representations)
* [How it works](#how-it-works)
* [How to test](#how-to-test)
* [Issues](#issues)
* [Changelog](#changelog)
* [Licence](#licence)
* [Contact](#contact)

<!-- vim-markdown-toc -->

## Introduction

KmerCamel🐫 is a tool for efficiently representing a set of k-mers a [masked superstring](https://doi.org/10.1101/2023.02.01.526717).

It is based on the following paper:

> Ondřej Sladký, Pavel Veselý, and Karel Břinda: Masked superstrings as a unified framework for textual *k*-mer set representations. *bioRxiv* 2023.02.01.526717, 2023.
[https://doi.org/10.1101/2023.02.01.526717](https://doi.org/10.1101/2023.02.01.526717)

See [supplementary materials](https://github.com/karel-brinda/masked-superstrings-supplement) of the aforementioned paper for experimental results with KmerCamel🐫.

The computation of masked superstring using KmerCamel🐫 is done in two steps -
first a superstring is computed with its default mask and then its mask can be optimized.

The computation of the masked superstring works as follows. KmerCamel🐫 reads an input FASTA file (optionally `gzip`ed), retrieves the associated k-mers (with supported $k$ up to 127), and outputs
a fasta file with a single record - a masked-cased superstring, which is in the nucleotide alphabet with case of the letters determining the mask symbols.
KmerCamel🐫 implements two different algorithms for computing the superstring:
global greedy and local greedy. Global greedy produces more compact superstrings and therefore is the default option,
but local greedy requires less memory and hence can be more suitable in use cases where memory is the main limitation.

To compute masked superstrings takes about 4-10s / 1M k-mers, which means about 5h to compute masked superstrings for the human genome.

All algorithms can be used to either work in the unidirectional model or in the bidirectional model
(i.e. treat $k$-mer and its reverse complement as the same; in this case either of them appears in the result).

Additionally, KmerCamel🐫 can optimize the mask of the superstring via the `optimize`subcommand. The implemented mask optimization algorithms are the following:
- Minimize the number of 1s in the mask.
- Maximize the number of 1s in the mask.
- Minimize the number of runs of 1s in the mask.

## Prerequisites

* GCC
* Zlib
* GLPK (can be installed via `apt-get install libglpk-dev` on Ubuntu or `brew install glpk` on macOS)

## Getting started

### Installation

Download and compile KmerCamel🐫 by running the following commands:

```
git clone --recursive https://github.com/OndrejSladky/kmercamel
cd kmercamel && make
```

Alternatively, you can install KmerCamel from Bioconda:
```
   conda install kmercamel
```

### Compression for k-mer set storage

```
./kmercamel -p yourfile.fa -k 31 -c > ms.fa        # Compute MS with the default mask
cat ms.fa | tr acgt 0000 | tr ACGT 1111 > mask.txt # Extract mask
cat ms.fa | tr acgt ACGT > superstring.txt         # Extract superstring
bzip2 --best mask.txt
xz -T1 -9 superstring.txt
```

For a super efficient compression of the superstring (often <2 bits / bp), you use some of the specialized tools based on statistical compression such as [GeCo3](https://github.com/cobilab/geco3) or [Jarvis3](https://github.com/cobilab/jarvis3).


### k-mer set indexing

Example with [FMSI](https://github.com/OndrejSladky/fmsi/activity?ref=main):
```
kmercamel -p yourfile.fa -k 31 -c > ms.fa          # Compute MS with the default mask
kmercamel optimize -p ms.fa -k 31 -c -o ms-opt.fa  # Maximize the number of 1s in the mask
fmsi index -p ms-opt.fa                            # Create a k-mer index
```

## Detailed instructions

Computing masked superstrings:
```
./kmercamel -p ./spneumoniae.fa -k 31 -c                # From a fasta file
./kmercamel -p - -k 31 -c                               # Read from stdin
./kmercamel -p ./spneumoniae.fa.gz -k 31 -c             # From a gzipped fasta file
./kmercamel -p ./spneumoniae.fa -k 127 -c               # Largest supported k
./kmercamel -p ./spneumoniae.fa -k 31 -a local -d 5 -c  # Use local greedy
./kmercamel -p ./spneumoniae.fa -k 31 -c -o out.fa      # Redirect output to a file
./🐫 -p ./spneumoniae.fa -k 31 -c                        # An alternative if your OS supports it
```

Optimizing masks:
```
./kmercamel optimize -p ./masked-superstring.fa -k 31 -a runs -c        # Minimize the number of runs of 1s
./kmercamel optimize -p ./masked-superstring.fa -k 31 -a ones -c        # Maximize the number of 1s
./kmercamel optimize -p ./masked-superstring.fa -k 31 -a zeros -c       # Maximize the number of 0s
./kmercamel optimize -p ./masked-superstring.fa -k 31 -a runapprox -c   # Approximately minimize the number of runs of 1s
```

Compute lower bound on the minimum possible superstring length of a k-mer set:
```
./kmercamel -l -p ./spneumoniae.fa -k 31
```

Additionally, KmerCamel🐫 experimentally implements both algorithms in their Aho-Corasick automaton versions. To use them, add `AC` to the algorithm name.
Note that they are slower than the original versions, but they can handle arbitrarily large *k*s.

### Arguments

The program has the following arguments:

- `-p path_to_fasta` - the path to fasta file (can be `gzip`ed). This is a required argument.
- `-k value_of_k` - the size of one k-mer (up to 127). This is a required argument.
- `-a algorithm` - the algorithm which should be run. Either `global` or `globalAC` for Global Greedy, `local` or `localAC` for Local Greedy.
The versions with AC use Aho-Corasick automaton. Default `global`.
- `-o output_path` - the path to output file. If not specified, output is printed to stdout.
- `-d value_of_d` - d_max used in Local Greedy. Default 5. Increasing `d` beyond `k` has no effect.
- `-c` - treat k-mer and its reverse complement as equal.
- `-l` - compute lower bound on the superstring length instead of the superstring.
- `-m` - turn off memory optimizations for `global`.
- `-h` - print help.
- `-v` - print version.

For mask optimization, run the subcommand `optimize` with the following arguments:

- `p path_to_fasta` - the path to fasta file (can be `gzip`ed). This is a required argument.
- `k k_value` - the size of one k-mer. This is a required argument.
- `a algorithm` - the algorithm for mask optimization. Either `ones` for maximizing the number of 1s, `runs` for minimizing the number of runs of 1s, `runsapprox` for approximately minimizing the number of runs of 1s, or `zeros` for maximizing the number of 0s. Default `ones`.
- `o output_path` - the path to output file. If not specified, output is printed to stdout.
- `c` - treat k-mer and its reverse complement as equal.
- `h` - print help.
- `v` - print version.


### Converting k-mer set superstring representation to the (r)SPSS representations

We provide a Python script for converting any masked superstring to a (r)SPSS representation.
Run `./convert_superstring.py < input.fa`. This runs a Python script which inputs a fasta file with masked-cased superstring and outputs the (r)SPSS representation.

## How it works

For details about the algorithms and their implementation, see the [Code README](./src/README.md).

## How to test

To ensure correctness of the results, KmerCamel🐫 has two levels of tests - unit tests and file-specific integration tests.

For integration tests  install [jellyfish (v2)](https://github.com/gmarcais/Jellyfish)
and add it to PATH.

You can verify all the algorithms for `1 < k < 128` on a *S. pneumoniae* by running `make verify`.
To run it on another dataset, see the [verification script](./verify.py).

You can run the C++ unittests by `make cpptest`.

To run all the test, simply run `make test`.

## Issues

Please use [Github issues](https://github.com/OndrejSladky/kmercamel/issues).


## Changelog

See [Releases](https://github.com/OndrejSladky/kmercamel/releases).


## Licence

[MIT](https://github.com/OndrejSladky/kmercamel/blob/master/LICENSE.txt)


## Contact

[Ondrej Sladky](https://iuuk.mff.cuni.cz/~sladky/) \<ondra.sladky@gmail.com\>\

