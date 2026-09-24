# Program-screen diagnostics

`run_secondary_diagnostics.R` reproduces the 118 within-year program-versus-control contrasts, checks residuals and influence, and reports HC3 covariance and 512 single-plot omission refits. The original endpoint/year Holm families are retained.

These are sensitivity checks for transformed-percentage program + block screens, not replacements for the primary mixed models. Original results remain in `results/descriptive/`. `screen_sensitivity_summary.csv` summarizes changes; detailed comparisons and residuals are supplied alongside it. See the root REPRODUCE.md for execution.
