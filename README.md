# Incomplete Markets

This repository contains MATLAB code for solving a continuous time incomplete markets model under different equilibrium concepts using the Achdou et al. (2022) method. The model can be solved in:

- Partial equilibrium
- General equilibrium Aiyagari-style with capital
- General equilibrium Huggett-style with bonds

The code builds heavily on the Achdou et al. (2022) codes available on [Benjamin Moll's website](https://benjaminmoll.com/codes/). This repository serves as a clean set of codes for my personal use for building on in future projects, and is provided as-is in case it is useful. Of minor interest is the procedure for building the sparse A matrix, which is scales well to larger problems. 

## Structure

- `Matlab/`: Model code and calibration routines.
- `Paper/`: LaTeX source and compiled PDF describing the model.

## Getting Started

1. Open `Matlab/main1_cali.m`.
2. Choose the equilibrium mode using the provided switch.
3. Run the script to solve the model and generate output.

## Documentation

A detailed description of the model and solution methods is available in the [model PDF](https://github.com/alexclymo/incomplete-markets/blob/main/Paper/model.pdf).
