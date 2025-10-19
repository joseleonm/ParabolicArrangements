# Parabolic DGA for Real Parabolic Arrangements

[![Launch in Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/joseleonm/ParabolicArrangements/HEAD)
[![Open demo notebook](https://img.shields.io/badge/launch-demo%20notebook-blue?logo=jupyter)](https://mybinder.org/v2/gh/joseleonm/ParabolicArrangements/HEAD?labpath=demo_parabolic_dga.ipynb)

Self-contained **SageMath 10.7** module implementing the relative parabolic chain complex,
integer/GF(2) boundary operators, Betti number computation, and a minimal
intersection product compatible with the cohomological cup product.

**Reference:**  
J. Cantarero – J. L. León-Medina, *The Cohomology Ring of Real Parabolic Arrangements*.

---

## 🔧 Requirements

This module was designed and tested in **SageMath 10.7**.

To run locally, install SageMath ≥ 10.7 from  
👉 [https://www.sagemath.org/download.html](https://www.sagemath.org/download.html)

Alternatively, you can try it **directly online** using Binder (no installation required):  
→ [Launch interactive Binder environment](https://mybinder.org/v2/gh/joseleonm/ParabolicArrangements/HEAD)

---

## 🚀 Quick Start (Local)

### Run doctests
```bash
sage -t parabolic_dga.sage
````

### Minimal usage (Sage console / Jupyter)

```python
load('parabolic_dga.sage')
W, P, Plist = build_W_P("A", 5)
Delta = ideal_non_k_equal_A(W, Plist, k=3)   # non-3-equal (type A₅)
D = ParabolicDGA(W, P, Plist, Delta)
D.self_test(trials=30, check_betti=True)      # d²=0, Leibniz, and A₅ sanity
D.betti_numbers(GF(2))                       # expected (1, 111, 20) on keys 0,1,2
```

---

## 📓 Notebook Demo

A complete example is available in
**`notebooks/demo_parabolic_dga.ipynb`**,
which demonstrates:

* building the arrangement,
* running internal consistency tests,
* computing Betti numbers over GF(2) and ℤ,
* testing the minimal intersection product.

You can open it instantly online:
👉 [**Launch demo notebook in Binder**](https://mybinder.org/v2/gh/joseleonm/ParabolicArrangements/HEAD?labpath=notebooks%2Fdemo_parabolic_dga.ipynb)

---

## 📁 Files

* `parabolic_dga.sage` — core module (self-contained, with doctests)
* `notebooks/demo_parabolic_dga.ipynb` — reproducible demonstration

---

## 🧾 License

MIT License (see `LICENSE`).

---

## 📚 Citation

If you use this software, please cite:

> J. Cantarero – J. L. León-Medina,
> *The Cohomology Ring of Real Parabolic Arrangements*.
