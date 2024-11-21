

# README for Indicator-DeepKriging

## Overview

This repository contains initial code blocks for implementing DeepKriging for indicator variables. Follow the steps below to run the scripts in the correct order:

## Steps to Run

1. **Generate data**  
   ```bash
   Rscript data_generation.R
   ```

2. **Perform 2D-to-1D Projection**  
   ```bash
   Rscript 2d-to-1d_Projection.R
   ```

3. **Run the DeepKriging model**  
   ```bash
   python indicator_deepkriging.py
   ```

4. **Compute and compare results**  
   ```bash
   Rscript Gaussian_Kriging.R
   Rscript Comute-Comparison_results.R
   ```
Note that, running the python codes are comptationally efficient, however, if you run the `Gaussian_Kriging.R` file, it wil require huge computation time and a system with enough number of cores. 

## Additional Resources

There are additional codes available in the `raf_codes/` folder. You can explore this folder to review the following files:

- **`bivariate_cdf.R`**: Focuses on calculating the bivariate cumulative distribution function (CDF).  
- **`projected_script1.R`** and **`projected_script2.R`**: Both scripts explore meaningful 2D-to-1D projections of bivariate data.

Feel free to check out these scripts for further insights and experiments.
