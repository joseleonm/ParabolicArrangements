# Parabolic DGA for Real Parabolic Arrangements

Self-contained Sage module implementing the relative parabolic chain complex,
integer/GF(2) boundary operators, Betti number computation, and a minimal
intersection product that matches the cohomological cup product.

**Reference:**  
J. Cantarero – J. L. León-Medina, *The Cohomology Ring of Real Parabolic Arrangements*.

## Quick Start

### Run doctests
```bash
sage -t parabolic_dga.sage
````

### Minimal usage (Sage console / Jupyter)

```python
load('parabolic_dga.sage')
W,P,Plist = build_W_P("A",5)
Delta = ideal_non_k_equal_A(W, Plist, k=3)  # non-3-equal (type A5)
D = ParabolicDGA(W,P,Plist,Delta)
D.self_test(trials=30, check_betti=True)     # d^2=0, Leibniz, and A5 sanity
D.betti_numbers(GF(2))                       # expected (1, 111, 20) on keys 0,1,2
```

## Notebook

See **`notebooks/demo_parabolic_dga.ipynb`** for an end-to-end example:

* build the arrangement,
* run internal tests,
* compute Betti numbers (GF(2), ZZ),
* try the minimal intersection product.

## Files

* `parabolic_dga.sage` — core module (autocontained, with doctests)
* `notebooks/demo_parabolic_dga.ipynb` — reproducible demonstration

## License

MIT (see `LICENSE`).

## Citation

If you use this software, please cite:

> J. Cantarero – J. L. León-Medina, *The Cohomology Ring of Real Parabolic Arrangements*.