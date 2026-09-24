# Black rot field trials, 2021–2025

Data and annotated analysis code accompanying *Five-year field evaluation of biological fungicide programs for black rot management in ‘Niagara’ juice grape* (Kaitlin M. Gold and David B. Combs).

## Download and run

[Download the complete code and data package](Supplementary_File_S1_v5.2.zip). Extract it and open the `BlackRot_GitHub_2026-09-24` folder. Its `README.md`, `REPRODUCE.md`, and `METHODS.md` describe the data, software requirements, commands, and analysis interpretation.

From the extracted package directory, run:

```sh
python3 scripts/check_data.py
```

The package contains 173 files: cluster observations from 264 plots and 66 program-year identities, application codebook, data dictionary, annotated R/Python scripts, recorded model and bootstrap results, figure sources, NEWA weather tables, and SHA-256 checksums. It preserves the recorded statistical results; model objects can be generated locally using the supplied commands.

The complete folder structure is inside the ZIP. Download and extract that archive before running scripts. Existing application-rate and formulation uncertainties are documented in the package. Cite the accompanying manuscript when using the study data; its final citation or DOI will be added when available.
