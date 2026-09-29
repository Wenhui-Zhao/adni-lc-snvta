# LC and SN-VTA integrity quantified using dual-echo PD-T2-weighted MRI is associated with cognition 10 years later

Analysis code accompanying the paper by Zhao et al.

## Scripts

Start with `run.sh`, which runs the each analyses. The individual analysis and helper scripts do not need to be launched separately.

Install the software below, make `Rscript`, `python3`, and `matlab` available on your `PATH`, and place `rri_boot_check.m` in the `matlab/` subfolder. From the repository directory, replace the paths below and run:

```bash
export PLS_TOOLBOX_ROOT="/path/to/pls_toolbox"

# Check the prepared input tables and software.
bash run.sh check --data "/path/to/prepared_data" --out "/path/to/output_directory"

# Run all analyses locally.
bash run.sh all --data "/path/to/prepared_data" --out "/path/to/output_directory" --mode local
```

To run one figure's analyses and required dependencies, replace `all` with `figure1`, `figure2`, `figure3`, `figure4`, or `figure5`. Summary CSV tables are saved in `<output_directory>/results/`. Keep the output directory outside the repository and separate from the input-data directory.

## Software

The analysis environment is:

- **Linux with Bash**
- **R 4.4.1** with `data.table`, `lme4`, `lmerTest`, `mgcv`, `MASS`, `Matrix`, `digest`, `jsonlite`, and `reformulas`. Mediation additionally requires the study-validated **`lavaan` 0.7.2** library.
- **Python 3.9.21** with **NumPy 2.0.2**, **pandas 2.2.3**, and **SciPy 1.13.1**.
- **MATLAB R2019b** with the **McIntosh PLS toolbox** for PLSC analyses.
