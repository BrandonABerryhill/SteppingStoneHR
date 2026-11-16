# Unstable heteroresistance as a stepping stone to stable resistance

## About this repository

This repository holds the code, the simulation results, and the source data behind the preprint

Berryhill, B. A. and Andersson, D. I. (2026). Unstable Heteroresistance as a Stepping Stone to Stable Resistance via Compensatory Mutations. bioRxiv. https://doi.org/10.64898/2026.10.08.756632

The model follows twenty bacterial states (R0 to R9 carry 0 to 9 copies of a β-lactamase gene; R0C to R9C are the same states with a compensatory mutation), a limiting resource, an antibiotic and β-lactamase under serial transfer. It is integrated by fixed-step Euler steps with a Monte Carlo layer for the transitions between states. There are two implementations of the same model: an R reference engine and a compiled C++ engine, the R package `hrsim`. The simulations in the paper were run with the compiled engine, and a subset with the R reference.

One R Markdown notebook, `SteppingStone.Rmd`, is the entry point. By default it rebuilds the seven data figures and the source data workbook from the shipped results. With one parameter changed it runs the simulations again first.

## Contents

| Path | What it is |
|---|---|
| `README.md` | This file |
| `LICENSE` | The MIT licence |
| `CITATION.cff` | Citation metadata for the preprint |
| `.gitignore` | Files a knit or an R session writes that are not part of the repository |
| `SteppingStone.Rmd` | The notebook: settings, engine comparison, batches, figures, workbook |
| `R/simulation/` | The model and the simulation pipeline, sourced by the notebook in file order: `00_config.R` (classification rule, paths), `01_parameters.R` (states, founder population, `make_condition()`), `02_transitions.R` (the 150 transitions), `03_engine_R.R` (R reference engine), `04_classify.R` (outcome classifier, frequencies, Clopper-Pearson intervals), `05_run.R` (batch runner with checkpoint and resume), `06_sensitivity.R` (Latin hypercube sweep and PAWN indices), `07_engines.R` (selects the R or the compiled engine), `08_notebook.R` (the notebook's steps: engine comparison, load or run a batch, frequency summary, sweep, run dates) |
| `R/figures/` | `00_style.R` (fonts, palette, theme, sizes, save helper); `01_panel_data.R` (computes the 17 panel data frames from `results/`); `fig1.R` to `fig4.R`, `figS3.R`, `figS4.R`, `figS5.R` (one script per figure); `source_data.R` (writes the workbook) |
| `hrsim/` | Source of the compiled engine package: `DESCRIPTION`, `NAMESPACE`, `LICENSE`, `LICENSE.md`, `R/`, `src/` (`simulate_hr.cpp`, `Makevars`, `RcppExports.cpp`), `man/`, `tests/` (verification tests; in `tests/ref/`, the R reference implementation the compiled engine was tested against; fixtures in `tests/fixtures/`), `inst/CITATION` |
| `hrsim_0.1.0.tar.gz` | The same package as a source tarball, built from `hrsim/` with `R CMD build` |
| `results/` | The shipped simulation results the notebook reads by default (layout below) |
| `figures/` | `Fig1.png` to `Fig4.png` and `FigS3.png` to `FigS5.png`, written by the notebook |
| `source_data/` | `Source_Data.xlsx`, the source data workbook, and `panel_data.rds`, the 17 panel data frames; both written by the notebook |
| `tests/resimulate.R` | A short test of the notebook's re-simulation mode |
| `session_info.txt` | R and package versions of the knit that wrote the shipped figures and workbook |

Layout of `results/`:

| Path | What it is |
|---|---|
| `results/run_dates.csv` | The date on which each batch was run, and the source of that date |
| `results/frequency_summary.rds` | Outcome counts, percentages and 95% Clopper-Pearson intervals of the frequency batches |
| `results/S1/`, `results/S2/`, `results/omega0/` | One folder per batch, each with `sims.rds` (one data frame per replicate, named by seed) and `condition.rds` (the parameters). Scenario 1: `S1_freq_cpp`, `S1_freq_R`, `S1_traj300`, `S1_long`. Scenario 2: `S2_freq_cpp`, `S2_freq_R`, `S2_long`. Control without compensation (ω = 0): `omega0_freq_cpp`, `omega0_freq_R`, `omega0_traj300`, `omega0_long` |
| `results/sensitivity/sweep_S1/` | The sensitivity sweep: `sweep.rds` and `sweep.csv` (300 parameter sets with their outcome frequencies) and `pawn.rds` (PAWN indices) |
| `results/long_run/` | `long_run_raw.csv` and `long_run_percycle.csv`, the series behind Supplementary Fig. 4 (see below) |

Batch names: `freq_cpp` is 4,000 replicates on the compiled engine and `freq_R` 100 replicates on the R reference, both run to 300 h and stored on the 12 transfer rows from 24 h to 288 h; `traj300` is 100 replicates to 300 h recorded from 0 h; `long` is 100 replicates to 3,000 h. Outcomes are classified from the stored replicates each time the notebook runs; classifications are not shipped.

## System requirements

- R. The version tested is 4.6.0.
- A C++17 compiler for the `hrsim` package. On Windows this is Rtools for the installed R version; on macOS, the Xcode command line tools.
- Pandoc, which `rmarkdown` needs to knit the notebook (RStudio includes it).
- The Arial font for the figures. Without it the figure scripts fall back on the default sans-serif font and say so.
- The R packages below.

One configuration was tested: an Apple silicon Mac (14 cores) running macOS 27.0.1, R 4.6.0 (aarch64-apple-darwin23), Apple clang 21.0.0 and pandoc 3.8.3. These are the package versions tested:

| Package | Version tested | Needed for |
|---|---|---|
| ggplot2 | 4.0.3 | Figures |
| patchwork | 1.3.2 | Figures |
| dplyr | 1.2.1 | Panel data, figures |
| tidyr | 1.3.2 | Panel data, figures |
| scales | 1.4.0 | Figures |
| ragg | 1.5.2 | PNG output of the figures |
| systemfonts | 1.3.2 | Font lookup of the figures |
| openxlsx | 4.2.8.1 | Workbook |
| zip | 2.3.3 | Unpacking and repacking the saved workbook file (installed with openxlsx, which imports it) |
| rmarkdown | 2.31 | Knitting the notebook |
| knitr | 1.51 | Knitting the notebook; used by the re-simulation test |
| Rcpp | 1.1.1-1.1 | Building and running `hrsim` |
| hrsim | 0.1.0 | Compiled engine (this repository) |
| future | 1.70.0 | Re-simulation only |
| furrr | 0.4.0 | Re-simulation only |
| lhs | 1.3.0 | Re-simulation only |

The notebook checks for the packages it needs when it starts and stops with a message naming any that are missing. The check does not name zip: openxlsx cannot be loaded without it, so if zip is missing the message names openxlsx. It installs nothing.

## Installation

Install the R packages from CRAN, in R:

```r
install.packages(c("ggplot2", "patchwork", "dplyr", "tidyr", "scales", "ragg",
                   "systemfonts", "openxlsx", "rmarkdown", "knitr", "Rcpp",
                   "future", "furrr", "lhs"))
```

Then install the compiled engine from the tarball, in a terminal, from the repository root:

```sh
R CMD INSTALL hrsim_0.1.0.tar.gz
```

On the tested configuration, compiling and installing `hrsim` took 2 s (two measurements, 2.1 and 2.2 s, each into an empty library). The time to install the CRAN packages was not measured; it depends on the network and on whether binary packages are available.

`hrsim` is optional for the quick start below: without it the notebook skips the engine comparison and still builds the figures and the workbook.

## Quick start: figures and workbook from the shipped results

From the repository root:

```sh
Rscript -e 'rmarkdown::render("SteppingStone.Rmd")'
```

This reads `results/`, classifies the stored replicates, runs the 24 h engine comparison if `hrsim` is installed, and writes the seven figures to `figures/`, `panel_data.rds` and `Source_Data.xlsx` to `source_data/`, `session_info.txt`, and the knitted document `SteppingStone.html`. The figure scripts also try to save a PDF copy of each figure beside its PNG file (`figures/Fig1.pdf` and so on) through R's cairo PDF device; if that device cannot be opened or the export stops with an error, the copy is skipped with a message. On the tested configuration the device could not be opened and no PDF copy was written. Existing files with these names are overwritten. It runs no stochastic simulation.

Render the notebook as shown, into R's global environment, which is the default of `rmarkdown::render()` called from a fresh R session and of RStudio's Knit button: the notebook sources its scripts into the global environment, and they read the notebook's settings there.

Run time on the tested configuration, with `hrsim` installed: 29.2 s and 29.7 s (two measurements, each a knit of a fresh copy of the repository).

## Re-running the simulations

From the repository root:

```sh
Rscript -e 'rmarkdown::render("SteppingStone.Rmd", params = list(resimulate = TRUE))'
```

This runs the batches behind the figures into `results_resimulated/`, then builds the figures, the workbook and `session_info.txt` from that folder, overwriting the ones built from the shipped results. The workbook's run dates are then those of the new run. The folder `results/` is not changed. The batch runner saves a checkpoint after each wave of replicates, so a knit that is interrupted continues from the saved replicates when it is started again with the same output folder. The notebook's parameters `resim_dir` (output folder, default `results_resimulated`) and `n_workers` (parallel workers, default the number of cores minus 2) can be set in the same `params` list. The number of workers does not change any result: each replicate's random stream depends only on its own seed.

A re-simulation is long, and the full run has not been timed. The table gives the number of runs in each batch and the time of one run measured on the tested configuration, one run at a time at low load: a compiled run to 300 h took 11.3 s (two later timings: 11.6 and 11.8 s), an R-reference run to 300 h 126 s, and a compiled run to 3,000 h 107 s. Under load from other jobs the same runs took 12.0, 159 and 152 s. The last column is arithmetic, not a measurement: runs × time per run ÷ 12 workers. It leaves out the time between runs and the uneven last wave when the runs do not divide evenly among the workers. The sweep column uses the time of a Scenario 1 run; runs with the other parameter sets of the sweep were not timed.

| Batch | Engine and horizon | Runs | Time per run | Runs × time ÷ 12 workers |
|---|---|---|---|---|
| Engine comparison | R reference and compiled, 300 h, random layer off | 2 | see below | not divided |
| `S1_freq_cpp`, `S2_freq_cpp`, `omega0_freq_cpp` | Compiled, 300 h | 4,000 each | 11.3 s | 3,767 s (63 min) each |
| `S1_freq_R`, `S2_freq_R`, `omega0_freq_R` | R reference, 300 h | 100 each | 126 s | 1,050 s (18 min) each |
| `S1_traj300`, `omega0_traj300` | Compiled, 300 h | 100 each | 11.3 s | 94 s each |
| `S1_long`, `S2_long`, `omega0_long` | Compiled, 3,000 h | 100 each | 107 s | 892 s (15 min) each |
| Sensitivity sweep | Compiled, 300 h | 300 sets × 100 = 30,000 | 11.3 s | 28,250 s (7.8 h) |

The same arithmetic over all batches gives about 45,600 s, or 12.7 h, on 12 workers. The engine comparison at 300 h took 225 s on the tested configuration (one measurement, while other jobs were running).

### Test of the re-simulation path

`tests/resimulate.R` runs a cut-down set of batches through the functions the notebook uses in re-simulation mode, writing to a temporary folder, and compares the result with the shipped files. It runs the engine comparison at 300 h and halts if the largest relative error exceeds 1e-9. It compares, column by column:

- Scenario 1, seeds 1 to 20, with the transfer rows of the same seeds in `results/S1/S1_traj300/sims.rds`;
- Scenario 2, seeds 1 to 20, with the transfer rows from 24 h to 288 h of the same seeds in `results/S2/S2_long/sims.rds`;
- the control, seeds 1 to 20, with the same seeds in `results/omega0/omega0_freq_cpp/sims.rds`;
- one Scenario 1 run to 3,000 h, seed 1, with seed 1 in `results/S1/S1_long/sims.rds`;
- sweep sets 1 and 2, 100 seeds each, with rows 1 and 2 of `results/sensitivity/sweep_S1/sweep.rds`.

It also runs Scenario 1 seeds 1 and 2 on the R reference engine and checks that both are returned. It prints one PASS or FAIL line per comparison and exits with status 1 if any fails. It needs `hrsim`, future, furrr, lhs and knitr. From the repository root:

```sh
Rscript tests/resimulate.R
```

Run time: 12 min on the tested configuration with 12 workers (two measurements, 12.3 and 12.4 min, each the time the test prints on its last line).

## How the C++ engine is verified against the R reference

The compiled engine is `hrsim/src/simulate_hr.cpp`. The R reference engine that the notebook sources is `R/simulation/03_engine_R.R`; the package tests use the R reference implementation in `hrsim/tests/ref/`. The compiled engine is compared with the R reference in three places.

**The package tests.** `hrsim/tests/run_all_gates.R` runs four check scripts. The three that call the compiled code (`gate1b_check.R`, `gate1_check.R` and `gate3_check.R`) build the Scenario 1 condition, its initial state and the transition table with the files in `hrsim/tests/ref/`, and the stored values named below are in `hrsim/tests/fixtures/`. Each file's header gives what it verifies, its tolerance and the seeds it uses.

- `gate1b_check.R` compares the rates and transition probabilities of the first step with the values of the R reference expressions stored in `tests/fixtures/step1_internals.rds` (largest relative error at most 1e-9), and checks over seeds 1 to 200,000 that each transition fires with a frequency within 5 binomial standard errors of its probability.
- `gate1_check.R` runs the compiled engine with the random layer off for Scenario 1 over 24 h and 300 h and compares it, cell by cell, with R-reference trajectories stored in `tests/fixtures/`. A cell passes when |C++ − R| ≤ 1e-12 + rtol · |R|, with rtol = 1e-9 over 24 h and 1e-6 over 300 h.
- `gate3_check.R` checks invariants of recorded trajectories: densities and antibiotic are never negative, values are finite, and fired transitions add as many cells as they remove.
- `gate2_check.R` reads the outcomes of a stored sample of nine Scenario 1 runs and compares their frequencies with the reference percentages written in the script. It runs no simulation, and with nine seeds it is a weak test.

With `hrsim` installed, from the repository root:

```sh
cd hrsim
Rscript -e 'source("tests/run_all_gates.R")'
```

On the tested configuration this took 17 s and printed PASS for gates 1b, 1 and 3.

**The notebook's deterministic comparison.** The section "Engine verification" of `SteppingStone.Rmd` runs both engines for Scenario 1 with the random layer switched off and reports the largest relative error over all 25 output columns (time, the 20 states, resource, antibiotic, β-lactamase and the transfer flag) and all recorded rows. Default mode runs it to 24 h and only reports it. Re-simulation mode runs it to 300 h and halts if the error exceeds 1e-9. The Methods of the paper report the value obtained on the machine that ran the published simulations. On the configuration tested here the notebook's comparison gives a maximum relative error of 0, at 24 h and at 300 h.

**The statistical comparison.** For each condition, the outcome frequencies of the 4,000 compiled replicates are compared with those of the 100 R-reference replicates through their 95% Clopper-Pearson intervals. The notebook's section "Results in brief" prints both tables and, for each condition and outcome, whether the two intervals overlap. For the shipped batches they overlap for all three conditions and all three outcomes.

**Compiler flag.** `hrsim/src/Makevars` sets `-ffp-contract=off`, which stops the compiler from fusing a multiplication and an addition into a single fused multiply-add. A fused multiply-add rounds once, where the separate multiplication and addition round twice. The agreement reported in the Methods was obtained with this flag.

## Data behind Supplementary Fig. 4

The data behind Supplementary Fig. 4 come from an old iteration of the model which is not regenerated by this code. They are in `results/long_run/`: `long_run_raw.csv` holds the hourly series and `long_run_percycle.csv` its reduction to one row per 24 h cycle, from which the figure is drawn. The notebook reads these files in both modes, and in both modes it derives the per-cycle reduction again from the raw file with `check_long_run_reduction()` and stops if it differs from the shipped one.

## Provenance of the C++ engine

The C++ engine, `src/simulate_hr.cpp` in `hrsim/`, was translated from the R reference implementation (`R/simulation/03_engine_R.R`) to C++ with the assistance of a generative AI model (Claude Opus 4.8, Anthropic), as described in the Methods. The translation was verified against the R reference before any results were generated from it: deterministically, with the random layer disabled, and statistically, on outcome frequencies. The checks are in `hrsim/tests/`. The authors reviewed the code and take responsibility for it.

## Figure panels, scripts and workbook tabs

Each panel is drawn by its figure script from one data frame computed by `R/figures/01_panel_data.R`, and its tab in `source_data/Source_Data.xlsx` is written by `R/figures/source_data.R` from the same data frame. The workbook has an `Index` tab and 24 data tabs. Supplementary Figs 1 and 2 are schematics: they have no script and no tab.

| Panel | Script | Workbook tab | Data frame | Shipped results read |
|---|---|---|---|---|
| Fig. 1A | `R/figures/fig1.R` | `Fig 1A` | `Fig1A` | `S1/S1_traj300`, seed 2 |
| Fig. 1B | `R/figures/fig1.R` | `Fig 1B` | `Fig1B` | `S1/S1_traj300`, seed 7 |
| Fig. 1C | `R/figures/fig1.R` | `Fig 1C` | `Fig1C` | `frequency_summary.rds`; per-replicate values from `S1/S1_freq_cpp` |
| Fig. 1D | `R/figures/fig1.R` | `Fig 1D` | `Fig1D` | `omega0/omega0_traj300` |
| Fig. 2A | `R/figures/fig2.R` | `Fig 2A` | `Fig2A` | `S1/S1_long` |
| Fig. 2B | `R/figures/fig2.R` | `Fig 2B` | `Fig2B` | `S1/S1_long` |
| Fig. 2C | `R/figures/fig2.R` | `Fig 2C` | `Fig2C` | `S1/S1_long` |
| Fig. 2D | `R/figures/fig2.R` | `Fig 2D` | `Fig2D` | `S1/S1_long` |
| Fig. 3A | `R/figures/fig3.R` | `Fig 3A` | `Fig3A` | `S1/S1_long`, `S2/S2_long` |
| Fig. 3B | `R/figures/fig3.R` | `Fig 3B` | `Fig3B` | `S2/S2_long` |
| Fig. 3C | `R/figures/fig3.R` | `Fig 3C` | `Fig3C` | None: computed from the Scenario 2 parameters |
| Fig. 4A | `R/figures/fig4.R` | `Fig 4A` | `Fig4A` | `sensitivity/sweep_S1/pawn.rds` |
| Fig. 4B | `R/figures/fig4.R` | `Fig 4B` | `Fig4B` | `sensitivity/sweep_S1/pawn.rds` |
| Fig. 4C | `R/figures/fig4.R` | `Fig 4C` | `Fig4C` | `sensitivity/sweep_S1/pawn.rds` |
| Supplementary Fig. 3 | `R/figures/figS3.R` | `Supp Fig 3` | `FigS3` | `omega0/omega0_long` |
| Supplementary Fig. 4 | `R/figures/figS4.R` | `Supp Fig 4` | `FigS4` | `long_run/long_run_percycle.csv` |
| Supplementary Fig. 5A | `R/figures/figS5.R` | `Supp Fig 5A` | `FigS5` | `sensitivity/sweep_S1/sweep.csv` |
| Supplementary Fig. 5B | `R/figures/figS5.R` | `Supp Fig 5B` | `FigS5` | `sensitivity/sweep_S1/sweep.csv` |
| Supplementary Fig. 5C | `R/figures/figS5.R` | `Supp Fig 5C` | `FigS5` | `sensitivity/sweep_S1/sweep.csv` |
| Supplementary Fig. 5D | `R/figures/figS5.R` | `Supp Fig 5D` | `FigS5` | `sensitivity/sweep_S1/sweep.csv` |
| Supplementary Fig. 5E | `R/figures/figS5.R` | `Supp Fig 5E` | `FigS5` | `sensitivity/sweep_S1/sweep.csv` |
| Supplementary Fig. 5F | `R/figures/figS5.R` | `Supp Fig 5F` | `FigS5` | `sensitivity/sweep_S1/sweep.csv` |
| Supplementary Table 1 | None | `Supp Table 1` | None | The `states` vector of `R/simulation/01_parameters.R` |
| Supplementary Table 2 | None | `Supp Table 2` | None | `make_condition()` for Scenarios 1 and 2 |

Paths in the last column are relative to `results/`, except for the two tables.

## Licence

The code and the data in this repository are released under the MIT licence; see `LICENSE`. The `hrsim` package carries the same licence in `hrsim/LICENSE` and `hrsim/LICENSE.md`.

## Contact

Brandon A. Berryhill, Brandon.Berryhill@imbim.uu.se

## Disclosure

Claude Opus 4.8 (Anthropic) was used to convert the model from R to C++, which reduced the runtime of the model by roughly an order of magnitude. The converted code was verified against the R reference before any results were generated from it, as described above. The authors reviewed all model outputs using non-artificial intelligence. The files in this GitHub were marked-up/commented by the same AI to be human readable.

## How to cite

If you use this code or these data, please cite the preprint:

Berryhill, B. A. and Andersson, D. I. (2026). Unstable Heteroresistance as a Stepping Stone to Stable Resistance via Compensatory Mutations. bioRxiv. https://doi.org/10.64898/2026.10.08.756632

The same reference is in `CITATION.cff`, and `citation("hrsim")` returns it in R.
