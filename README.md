# Black rot field trials, 2021–2025

Data and annotated analysis code accompanying *Five-year field evaluation of biological fungicide programs for black rot management in ‘Niagara’ juice grape* (Kaitlin M. Gold and David B. Combs).

The cluster dataset contains **264 plots from 66 program-year identities** at one vineyard over five seasons. This repository provides the observations, application codebook, saved numerical results, and code for the manuscript analyses and figures.

| Folder | Contents |
|---|---|
| `data/` | Plot and cluster observations, application codebook, field dictionary, and data corrections |
| `scripts/` | Primary and secondary models, calibrated tests, descriptive summaries, and figure rendering |
| `results/` | Recorded numerical results and R session information |
| `analyses/` | Completed 2022 assignment sensitivity and program-screen diagnostics |
| `figures/` | Vector figure sources and Figure 4–5 input tables/code |
| `descriptive_newa_timeline/` | NEWA source tables, dates, counts, and Figure 6 code |

**Start:** `python3 scripts/check_data.py` checks the distributed observations and file checksums without fitting models. See [REPRODUCE.md](REPRODUCE.md) for installation and execution commands, [METHODS.md](METHODS.md) for analysis interpretation, and `data/data_dictionary.csv` for field definitions.

Saved estimates and simulation results are provided so figures and summaries can be inspected without repeating lengthy bootstrap runs. Model objects are generated locally by the analysis scripts. Run analyses in the separate `rerun/` directory to preserve the recorded results.

OxiDate was applied at 1 gal/acre in programs Y2022_P16 and Y2023_P07, as confirmed by the authors and recorded in the annual trial reports. Ensendo contains azoxystrobin and Pseudomonas chlororaphis AFS009. Curezin was preregistration material in the 2023–2024 trials; the current registration is [EPA 102559-1](https://www3.epa.gov/pesticides/chem_search/ppls/102559-00001-20250806.pdf), with copper sulfate pentahydrate and phosphorous acid as active ingredients and nutritional zinc. The 2025 trial report identifies MagnaBon as CS 2005. Supplementary Table S2 provides product identities and sources.

NEWA data are attributed to Cornell’s Network for Environment and Weather Applications in the timeline folder. Cite the accompanying manuscript when using the study data; add its final citation or DOI when available.
