# 2022 assignment sensitivity

Eight records with unresolved program/block assignment (Y2022_P02 and Y2022_P03) are removed, leaving 256 plots and 64 program-years. Primary results remain those from all 264 plots.

`run_additive_sensitivity.R` refits the original additive structures. `prepare_interaction_models.R` prepares reduced-data ML models. The two bootstrap scripts calibrate interaction tests. `summarize_public_sensitivity.py` creates the two summary tables and verifies simulation counts and plus-one probabilities from recorded outputs.

Pooled means changed by at most 2.160 incidence and 2.033 severity percentage points; no Holm class-contrast decision changed. Interaction bootstrap P values are 0.358871 (incidence; 1,983 valid of 1,999) and 0.079314 (severity; 9,972 valid of 9,999). The nonregular full incidence fit is not used for ordinary likelihood-ratio inference. This sensitivity does not recover ambiguous identifiers or establish equivalence. Commands are in the root REPRODUCE.md.
