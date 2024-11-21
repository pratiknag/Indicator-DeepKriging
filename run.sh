#!/usr/bin/bash 

Rscript data_generation.R
Rscript 2d-to-1d_Projection.R
python indicator_deepkriging.py
