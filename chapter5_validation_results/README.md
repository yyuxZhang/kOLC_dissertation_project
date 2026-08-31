# Fixed-Radius k-Orthogonal Line Center

This repository contains the implementation and experimental materials for the MSc dissertation:

**The Fixed-Radius k-Orthogonal Line Center Problem: Implementation and Extension**

**Author:** Yuxin Zhang  
**Programme:** MSc Advanced Computer Science  
**Institution:** University of Leeds  
**Academic year:** 2025/2026

## Project Overview

The project implements Algorithm 1 for the fixed-radius **k-Orthogonal Line Center (k-OLC)** problem.

For point inputs, horizontal and vertical candidate lines are constructed with spacing `2r`. Each point is assigned to one vertical and one horizontal candidate, producing an edge in a bipartite graph. A minimum vertex cover is then recovered from a maximum matching to select the output lines.

If the original point instance can be covered by `k` horizontal or vertical lines at radius `r`, the implemented routine returns at most `2k` lines at the same radius.

The implementation is extended to equal-radius disks using the effective centre threshold:

```text
R = r + rho
```

where `rho` is the common disk radius. The disk version therefore uses candidate-line spacing `2R` and verifies the returned solution using the original disk-to-line distance.

## Repository Structure

```text
kOLC_dissertation_project/
├── orthogonal_line_center_complete_experiments.ipynb
├── chapter5_independent_validation.py
├── kolc_experiment_results/
└── chapter5_validation_results/
```

### `orthogonal_line_center_complete_experiments.ipynb`

The main notebook contains:

- the point-input implementation of Algorithm 1;
- the equal-radius disk extension;
- sparse candidate-grid and bipartite-graph construction;
- built-in Hopcroft-Karp matching and minimum vertex-cover recovery;
- an exact endpoint-enumeration solver for small instances;
- point and disk visualisation;
- interactive comparison interfaces;
- deterministic tests;
- randomised cross-validation;
- approximation-quality experiments; and
- runtime and sparse-graph experiments.

### `chapter5_independent_validation.py`

This script provides an independent validation path using:

- a separate depth-first augmenting-path matching routine;
- independent minimum vertex-cover recovery;
- exhaustive minimum-cover checks for small bipartite graphs; and
- scalar point-to-line and disk-to-line distance calculations.

### `kolc_experiment_results/`

This directory contains the Chapter 4 figures, raw CSV results, grouped summaries and experiment metadata.

### `chapter5_validation_results/`

This directory contains the independent graph-validation and geometry-validation records, together with a JSON summary.

## Requirements

The project requires Python 3 and the packages listed in `requirements.txt`.

Install them with:

```bash
python -m pip install -r requirements.txt
```

The core solver uses NumPy and Matplotlib. `ipywidgets` is used for the interactive interfaces, while NetworkX is supported as an optional comparison backend.

## Running the Notebook

Open the notebook in Jupyter Notebook, JupyterLab or Google Colab:

```bash
jupyter notebook orthogonal_line_center_complete_experiments.ipynb
```

Run the definition cells before running the tests or experiments.

The notebook contains a smoke-test section for quick verification. Full experiments are disabled by default:

```python
RUN_EXPERIMENT_SUITE = False
FULL_EXPERIMENTS = False
```

To reproduce the complete dissertation experiment suite, change both values to:

```python
RUN_EXPERIMENT_SUITE = True
FULL_EXPERIMENTS = True
```

The generated files are written to `kolc_experiment_results/`.

## Running the Independent Validation

From the repository directory, run:

```bash
python chapter5_independent_validation.py \
  --notebook orthogonal_line_center_complete_experiments.ipynb
```

The outputs are written to `chapter5_validation_results/`.

## Experimental Configuration

The global random seed is:

```text
20260822
```

The reported evaluation includes:

- 200 randomised point/disk cross-validation records;
- 180 exact approximation-quality comparisons;
- 168 runtime measurements for input sizes up to 5,000 objects;
- 1,000 independent bipartite-graph checks, including 200 exhaustive small-graph checks; and
- 500 independent scalar geometry checks.

## Reported Results

Across the 180 exact comparisons, the mean line-count ratio was `1.362` and the median was `1.5`. Algorithm 1 matched the exact minimum in 65 cases, and no observed ratio exceeded the theoretical fixed-radius bound of `2`.

For the scalability experiment, the largest recorded median runtime at `n = 5,000` was below `8.3 ms` in the reported environment.

## Scope and Limitations

The implementation is limited to:

- two-dimensional Euclidean geometry;
- horizontal and vertical lines;
- a supplied fixed coverage radius; and
- disks with one common radius.

It does not search for the minimum possible radius for a fixed line budget, handle unequal-radius disks, or allow unrestricted line orientations.

## Academic Attribution

The point-input algorithm implemented in this repository is based on:

Das, A.K., Das, S. and Mukherjee, J. (2023) 'Approximation algorithms for orthogonal line centers', *Discrete Applied Mathematics*, 338, pp. 69-76. doi: 10.1016/j.dam.2023.05.014.

The Python implementation, equal-radius disk adaptation, exact small-instance baseline, visualisation, experiments and validation framework were developed for this MSc project.
