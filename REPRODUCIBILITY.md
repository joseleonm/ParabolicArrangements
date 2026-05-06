# Reproducibility Guide

This document explains how to reproduce the core computational results using a local SageMath setup or Binder.

## 1. Quick Start (Binder)

1. Open the Binder badge from the README.
2. Run notebooks in order:
   - notebooks/01_Cohomology_and_Asphericity.ipynb
   - notebooks/02_Fundamental_Groups_and_Orbifolds.ipynb
   - notebooks/03_Euler_Characteristic_and_Combinatorics.ipynb

## 2. Local Environment (Conda)

1. Create and activate the environment:

   conda env create -f environment.yml
   conda activate parabolic-arrangements

2. Start SageMath:

   sage

3. Load the module:

   load('parabolic_arrangements.sage')

## 3. Built-in Validation in Sage

Validation can be performed directly from the module in a Sage session, without running repository test scripts.

Minimal validation example:

load('parabolic_arrangements.sage')
W, Plist, _ = build_W_P('A', 4)
Delta = ideal_k_parabolic(W, Plist, k=3)
arr = ParabolicArrangement(W, Plist, Delta)
arr.verify_leibniz_rule(ring=GF(2), trials=10)

If this returns True and notebook computations reproduce expected dimensions/results, the setup is considered validated.

## 4. Reproducibility Notes

- GAP is required for the fundamental-group routines.
- Non-crystallographic finite Coxeter types are not currently supported by the implementation.
- For publication snapshots, cite a tagged release and its Zenodo DOI.
