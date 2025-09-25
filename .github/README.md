# GitHub Actions CI/CD

## Automated Testing

This repository includes automated testing via GitHub Actions that validates simulation results against golden files.

### Workflow: `validate-simulation.yml`

**Triggers:**
- Push to `main`, `master`, or `develop` branches
- Pull requests to `main`, `master`, or `develop` branches  
- Manual dispatch via GitHub UI

**What it does:**
1. Sets up MATLAB environment (R2023b)
2. Checks for required simulation files
3. Verifies golden file exists
4. Runs simulation validation with tolerance `1e-10`
5. Reports success/failure with appropriate exit codes

**Requirements:**
- Golden file must be committed: `goldenfiles/goldenfile_Np6_tmax020.mat`
- Required files: `main_2Dcrossing_3DHyQMOM35.m`, `simulation_plots.m`, `validate_against_goldenfile.m`

### Creating/Updating Golden Files

Golden files should be created locally and committed to the repository:

```matlab
% Create golden file
>> run_goldenfile_creation

% Commit the golden file
$ git add goldenfiles/
$ git commit -m "Update golden file"
$ git push
```

### Tolerance Settings

The CI uses a tolerance of `1e-10` which is slightly more lenient than the default `1e-12` to account for:
- Different MATLAB versions
- Different system architectures  
- Numerical precision variations in CI environment

### Troubleshooting Failed Tests

If validation fails:

1. **Check the GitHub Actions logs** for specific error messages
2. **Run validation locally** with the same tolerance:
   ```matlab
   >> validate_against_goldenfile(1e-10)
   ```
3. **If changes are intentional**, update the golden file and commit it
4. **If unexpected**, investigate the differences reported in the logs

### Artifacts

On failure, the workflow uploads simulation artifacts for debugging:
- Log files
- Generated `.mat` files  
- Golden files directory

These are retained for 7 days and can be downloaded from the GitHub Actions page.
