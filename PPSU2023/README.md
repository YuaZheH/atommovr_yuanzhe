# Engineering fast algorithms for the bottleneck matching problem

## Organization

The relavent source files are located in the top-level directory of the repository.
The source files for MatchMaker are located in the [*extern*](extern) folder.

The directory structure of the repository is as below.
```
.
├── bottleneckBipartiteMatching.c
├── bottleneckBipartiteMatching_runner.c
├── bvnGreedy.c
├── compareBisectionBasedOnMC64J3_Bottled_runner.c
├── doMatchingExps.c
├── extern
│   ├── cheap.c
│   ├── Makefile
│   ├── matching.c
│   └── matchmaker.h
├── Makefile
├── Makefile.inc
├── Makefile.inc.tmp
├── matrixUtils.c
├── matrixUtils.h
├── mc64a.F
├── mc64main.c
├── measureMMakerTime.c
├── mmio.c
├── mmio.h
├── myExpsForItersBttldBsctn.m
└── portdefs.h
```

## Building and Usage

### Prerequisites
* gcc version 10.2.1 or higher

### Building

- The codes are written in C

- Within the top-level directory of the repository, build the code by running
  ```
  make
  ```

### Usage

For running experiments, the executables of interest are:

- `doMatchingExps`: compares **mc64**, **bottled** and **thresh** methods for the bottleneck matching problem
- `bvnGreedy`: performs Birkhoff-von Neumann (BvN) decomposition using the proposed heuristic

<pre><code>
./doMatchingExps
matrix      : the input matrix
scaling     : 1 => no scaling; 2 => pattern scaling; 3 => value scaling, of the matrix
perm        : 1 => no permutation; 2 => column permutation; 3 => row permutation; 4 => column and row permutation, of the matrix
mc64-job-id : 2 => SAP; 3 => threshold-based 
</code></pre>

<pre><code>
./bvnGreedy
matrix      : the input matrix
scaling     : 1 => no scaling; 2 => pattern scaling; 3 => value scaling, of the matrix
perm        : 1 => no permutation; 2 => column permutation; 3 => row permutation; 4 => column and row permutation, of the matrix
[customZero]: a value close to 0, which is to be considered as 0. Presently ignored and computed in the code
[customOne] : a value close to 1, which is to be considered as 1. Presently ignored and computed in the code
[numPerm]   : number of requested permutation matrices 
</code></pre>

#### Example

```
./doMatchingExps ecology1.mtx 2 1 2
```
Produces the following ouput:

<pre><code>
Performing experiments on matrix ecology1.mtx with scaling 2 and permutation 1
...
read-in time 0.92
Should scale the matrix pattern
No permutation
...
total mc64 time 31.06
total bottled time 1.69, threshold 0.1999 iters 5
total thresh time 3.25, threshold 0.1999 iters 19
</code></pre>

```
./bvnGreedy atmosmodm.mtx 2 1 0.1 0.9 50
```
Produces the following ouput:

<pre><code>
...
read-in time 2.95
Should scale the matrix pattern
No permutation
the params zero, one, numf 0.000001 1.00 50
this is the current threshold 0.14286 (0.00000 1.00000)
total bvnGreedy time 202.03, numfactors computed 50, sum alphas 0.9997
</code></pre>
