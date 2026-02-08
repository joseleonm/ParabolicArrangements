# Parabolic Arrangements - Cohomology Ring Computation

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/joseleonm/ParabolicArrangements/HEAD?labpath=notebooks%2Fdemo_parabolic_dga.ipynb)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![SageMath](https://img.shields.io/badge/SageMath-10.7+-blue.svg)](https://www.sagemath.org/)

SageMath implementation for computing the cohomology ring of real parabolic arrangement complements.

## Table of Contents
- [Features](#features)
- [Quick Start](#quick-start)
- [Supported Arrangements](#supported-arrangements)
- [Key Properties](#key-properties)
- [Demo Notebook](#demo-notebook)
- [Citation](#citation)
- [Mathematical Background](#mathematical-background)
- [License](#license)
- [Acknowledgments](#acknowledgments)

## Features

- Compute Betti numbers and cohomology bases for parabolic arrangements
- Cup product computation using the zonotopal diagonal
- DGA structure validation (Leibniz rule)
- Support for all finite Weyl group types

## Quick Start

### Try Online (No Installation Required)

Click the Binder badge above to run the demo notebook in your browser!

### Local Installation

Requires [SageMath](https://www.sagemath.org/download.html) 10.7 or newer.

```bash
git clone https://github.com/joseleonm/ParabolicArrangements.git
cd ParabolicArrangements
sage
```

### Minimal Example

```python
sage: load('parabolic_arrangement_cohomology.sage')
sage: W, Plist, _ = build_W_P('A', 3)
sage: Delta = ideal_k_parabolic(W, Plist, k=3)
sage: cohomology = ParabolicArrangementCohomology(W, Plist, Delta)
sage: [len(cohomology.cohomology_basis(k, ring=GF(2)))
....:  for k in range(cohomology.max_grade + 1)]
[1, 14, 0, 0, 0]  # Betti numbers
```

### Basic Usage

```python
sage: load('parabolic_arrangement_cohomology.sage')
sage:
sage: # Build the 3-equal arrangement in type A_4
sage: # (excludes parabolic cosets where 3 consecutive coordinates are equal)
sage: W, Plist, _ = build_W_P('A', 4)
sage: Delta = ideal_k_parabolic(W, Plist, k=3)
sage: arrangement = ParabolicArrangementCohomology(W, Plist, Delta)
sage:
sage: # Compute cohomology ring over F_2
sage: for k in range(arrangement.max_grade + 1):
....:     basis = arrangement.cohomology_basis(k, ring=GF(2))
....:     print(f"dim H^{k} = {len(basis)}")
dim H^0 = 1
dim H^1 = 35
dim H^2 = 0
dim H^3 = 0
dim H^4 = 0
sage:
sage: # Verify DGA structure (Leibniz rule)
sage: arrangement.verify_leibniz_rule(ring=GF(2), trials=10)
--- Verifying DGA Leibniz Rule over Finite Field of size 2 ---
Testing degrees (1, 1)...
Leibniz Rule Verified: The DGA structure is consistent.
True
```

## Supported Arrangements

The mathematical theory applies to **any parabolic arrangement** (filter in the poset of parabolic cosets) for any finite Weyl group. The implementation provides a constructor for the **k-parabolic arrangements** introduced by Severs and White [[SW2012]](https://www.sciencedirect.com/science/article/pii/S0097316512000490):

| Type | Weyl Group | Constructor Example |
|------|------------|---------------------|
| A_n | Symmetric group S_{n+1} | `build_W_P('A', 5)` |
| B_n | Hyperoctahedral group | `build_W_P('B', 4)` |
| D_n | Even signed permutations | `build_W_P('D', 4)` |
| E_6, E_7, E_8 | Exceptional types | `build_W_P('E', 6)` |
| F_4, G_2 | Other exceptional types | `build_W_P('F', 4)` |

The **k-parabolic ideal** constructor `ideal_k_parabolic(W, Plist, k)` excludes parabolic cosets `wW_I` where the subset `I` contains a connected component of size ≥ k-1 in the Dynkin diagram. For type A, this corresponds to the classical k-equal arrangements.

**Note**: While `ideal_k_parabolic` is currently the primary constructor, users can define custom arrangements by specifying any filter `Delta` in the parabolic coset poset. The cohomology computation works for arbitrary parabolic arrangements.

## Key Properties

### Computational Efficiency
- Precomputed geometric structures for cup products are cached and reused
- Minimal coset representatives cached with LRU eviction policy
- Longest elements in parabolic subgroups computed once and cached
- Direct cellular computation avoids costly barycentric subdivision

### Mathematical Structure
- Exact cellular DGA model based on the restricted Coxeter permutahedron
- Cup product via zonotopal (Serre) diagonal transported from cubical structure
- Satisfies graded Leibniz rule: δ(u ∪ v) = (δu) ∪ v + (-1)^|u| u ∪ (δv)

## Demo Notebook

See `notebooks/demo_parabolic_dga.ipynb` for a complete walkthrough including:

- Building parabolic arrangements from Weyl groups
- Computing cohomology rings and Betti numbers
- Computing cup products using the zonotopal diagonal
- Verifying DGA structure (Leibniz rule)
- Examples in types A, B, D

## Citation

If you use this software in your research, please cite:

```bibtex
@misc{leon2025parabolic,
  author = {León-Medina, J. L. and Cantarero, J.},
  title = {The Topology of Real Parabolic Arrangements},
  year = {2025},
  note = {Manuscript in preparation}
}
```

## Mathematical Background

This implementation is derived from the mathematical theory developed in:

> **J. L. León-Medina and J. Cantarero**,
> *"The Topology of Real Parabolic Arrangements"*,
> manuscript in preparation.

The theory establishes that the complement of a real parabolic arrangement (a filter in the poset of parabolic cosets of a finite Weyl group) is homotopy equivalent to a **restricted Coxeter permutahedron** Perm_A(W). This provides an explicit cellular model via Björner duality between the Coxeter complex and the permutahedron.

The cohomology ring structure arises from a **zonotopal diagonal map** that transports the classical Serre diagonal from the cubical structure of the permutahedron (realized as a root zonotope) to the parabolic coset complex. This yields an explicit cup product formula:

    (u ∪ v)(wW_I) = Σ ε(J,K) · u(wW_J) · v(w·w_{0,J} W_K)

where the sum runs over ordered partitions I = J ⊔ K, w_{0,J} is the longest element in W_J, and ε(J,K) is the shuffle sign.

**Additional references**:
- C. Severs and J. White, "Constructions of k-parabolic subspace arrangements", European J. Combin. 33 (2012), 1002–1022. (k-parabolic arrangements)
- A. Björner and G. Ziegler, "Combinatorial stratification of complex arrangements", J. Amer. Math. Soc. 5 (1992), 105–149. (Björner duality)
- M. Davis, *The Geometry and Topology of Coxeter Groups*, Princeton Univ. Press, 2008.
- J.-P. Serre, "Homologie singulière des espaces fibrés. Applications", Ann. of Math. (1951). (cubical diagonal)

## License

MIT License - see [LICENSE](LICENSE) file.

## Acknowledgments

This code was developed as part of the project **SECIHTI CBF2023-2024-4059: "Interacciones topológico-computacionales"**.

*Note: Portions of the implementation were developed with AI assistance and have been rigorously verified against the mathematical theory.*

---

**Version**: 1.2 | **Updated**: February 2025
