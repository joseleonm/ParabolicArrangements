# Parabolic Arrangements

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/joseleonm/ParabolicArrangements/HEAD?labpath=notebooks%2F01_Cohomology_and_Asphericity.ipynb)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![SageMath](https://img.shields.io/badge/SageMath-10.7+-blue.svg)](https://www.sagemath.org/)

SageMath implementation for computing topological invariants of real parabolic arrangement complements, accompanying the paper:

> **J. Cantarero and J. L. León-Medina**,
> *"The Topology of Real Parabolic Arrangements"*.

## Table of Contents
- [Overview](#overview)
- [Quick Start](#quick-start)
- [What Is Implemented](#what-is-implemented)
- [Scope and Limitations](#scope-and-limitations)
- [Supported Arrangements](#supported-arrangements)
- [Repository Structure](#repository-structure)
- [Tutorial Notebooks](#tutorial-notebooks)
- [Citation](#citation)
- [Mathematical Background](#mathematical-background)
- [License](#license)
- [Acknowledgments](#acknowledgments)

## Overview

The module `parabolic_arrangements.sage` provides algorithms for computing topological invariants of the complement of a parabolic arrangement associated to a finite Coxeter group, using the restricted permutahedral complex as cellular model.

**Implemented invariants:**

- Betti numbers and integral cohomology bases
- Cup product (cohomology ring structure) for arrangements whose cells are *strictly commutative*
- Fundamental groups via Bridson–Haefliger complexes of groups
- Euler characteristics via orbit-weighted formulas
- Geometric asphericity criteria (triangle-free and flag criteria)

## Quick Start

### Try Online (No Installation Required)

Click the Binder badge above to launch the first tutorial notebook in your browser.

### Local Installation

Requires [SageMath](https://www.sagemath.org/download.html) 10.7 or newer.

```bash
git clone https://github.com/joseleonm/ParabolicArrangements.git
cd ParabolicArrangements
sage
```

### Minimal Example

```python
sage: load('parabolic_arrangements.sage')
sage: W, Plist, _ = build_W_P('A', 3)
sage: Delta = ideal_k_parabolic(W, Plist, k=3)
sage: arr = ParabolicArrangement(W, Plist, Delta)
sage: [len(arr.cohomology_basis(k, ring=GF(2)))
....:  for k in range(arr.max_grade + 1)]
[1, 14, 0, 0, 0]
```

### Basic Usage

```python
sage: load('parabolic_arrangements.sage')
sage:
sage: # Build the 3-equal arrangement in type A_4
sage: W, Plist, _ = build_W_P('A', 4)
sage: Delta = ideal_k_parabolic(W, Plist, k=3)
sage: arrangement = ParabolicArrangement(W, Plist, Delta)
sage:
sage: # Cohomology over GF(32003) — Betti numbers
sage: for k in range(arrangement.max_grade + 1):
....:     basis = arrangement.cohomology_basis(k, ring=GF(32003))
....:     print(f"dim H^{k} = {len(basis)}")
dim H^0 = 1
dim H^1 = 35
dim H^2 = 0
dim H^3 = 0
dim H^4 = 0
sage:
sage: # Verify the Leibniz rule (DGA structure)
sage: arrangement.verify_leibniz_rule(ring=GF(2), trials=10)
True
```

## What Is Implemented

### Cellular (co)homology

The coboundary operator of the restricted permutahedral complex is fully implemented for arbitrary parabolic arrangements of any finite Weyl group type.  Cohomology bases are computed by Smith normal form or rank reduction over any coefficient ring.

### Cup product and DGA structure

The paper establishes a general zonotopal diagonal approximation obtained as the cellular pushforward of the classical Serre diagonal on cubes:

    Δ_Z(wW_I) := (p_{I#} ⊗ p_{I#})(Δ_Serre(E_I))

where `p_I : □^{Φ_I⁺} → wW_I` is the canonical affine projection from the ambient hypercube to the zonotope.  This diagonal is theoretically well-defined for any finite Coxeter type and induces the correct cup product on cohomology.

However, evaluating this pushforward for a general non-commutative cell (one where `m_{ij} ≥ 3` for some pair in `I`) requires computing the cellular chain map of a degenerate projection — a task closely related to the Saneblidze–Umble diagonal on permutohedra — and is not implemented in this module.

**What is implemented** is the *cubical reduction* (Theorem: Cubical Diagonal Reduction in the paper): when a cell `wW_I` is **strictly commutative** (`m_{ij} = 2` for all `sᵢ, sⱼ ∈ I`), the projection `p_I` is a cellular isomorphism and the diagonal reduces exactly to the explicit combinatorial formula

    (u ∪ v)(wW_I) = Σ_{J⊔K=I} ε(J,K) · u(wW_J) · v(w·w_{0,J} W_K).

This covers the following important cases in full:

| Scope | Detail |
|-------|--------|
| **All 3-parabolic arrangements** | Every surviving cell is strictly commutative. This is the primary application of the paper. |
| **Right-angled arrangements** | Any arrangement retaining only rank-2 strata with `m_{st} = 2`. The complement is modeled entirely by commutative cells. |
| **Degree-1 products** | Cup products of degree-1 cochains always land on rank-2 cells, which are always cubes (dimension ≤ 2). |

For arrangements containing non-commutative cells of rank ≥ 3 (e.g., the full 4-parabolic arrangement in type D₅ with strata of type A₃), the cup product formula is applied at the cochain level but **does not yield a strict DGA**: associativity may fail at the cochain level.  The induced product on cohomology is still well-defined (it agrees with the singular cup product), but the Leibniz rule may not hold at the cochain level.

### Fundamental groups

The fundamental group π₁(K_A) is computed via three complementary routes:

- **W-route**: Modified Coxeter presentation via the exact sequence `1 → π₁ → Γ_A → W → 1`. Valid for W-invariant arrangements retaining all rank-1 strata.
- **B&H route**: Bridson–Haefliger complex-of-groups for H-invariant arrangements (H ≤ W). Computes group extensions with explicit transport elements.
- **CW 2-skeleton route**: Direct computation from the 2-skeleton. Always valid; robust fallback for non-invariant arrangements.

All three routes invoke GAP for coset enumeration and abelianization, so a working GAP installation is required.

### Geometric asphericity

Two criteria for the K(π,1) property are implemented:

- `check_flag_criterion()`: Checks whether every vertex link is a metric flag complex.
- `is_triangle_free()`: For W-invariant arrangements removing all rank ≥ 3 strata, checks the triangle-free criterion on the retained graph Γ_A.

### Euler characteristic

The orbit-weighted formula `χ = |W| Σ_{J ∈ T} (-1)^|J| / |W_J|` for W-invariant arrangements is implemented in `euler_characteristic_W_invariant()`, alongside a direct cellular computation `euler_characteristic()` valid for any arrangement.

## Scope and Limitations

| Feature | Status |
|---------|--------|
| Coboundary matrix / Betti numbers | Fully implemented for all types |
| Cup product (commutative cells) | Fully implemented; strict DGA |
| Cup product (general cells) | Not implemented (requires Saneblidze–Umble diagonal) |
| Fundamental groups | Fully implemented (requires GAP) |
| Euler characteristic | Fully implemented |
| Asphericity criteria | Fully implemented |
| Non-crystallographic types (H₃, H₄, I₂(m)) | Not supported by current code (theory applies) |

## Supported Arrangements

The mathematical theory applies to **any parabolic arrangement** (filter in the poset of parabolic cosets) for any finite Coxeter group.  The implementation covers all finite **Weyl group** types:

| Type | Weyl Group | Constructor |
|------|------------|-------------|
| A_n | Symmetric group S_{n+1} | `build_W_P('A', n)` |
| B_n | Hyperoctahedral group | `build_W_P('B', n)` |
| D_n | Even signed permutations | `build_W_P('D', n)` |
| E_6, E_7, E_8 | Exceptional types | `build_W_P('E', n)` |
| F_4, G_2 | Other exceptional types | `build_W_P('F', 4)`, `build_W_P('G', 2)` |

The **k-parabolic ideal** constructor `ideal_k_parabolic(W, Plist, k)` generates the standard k-parabolic arrangements introduced by Severs and White.  Custom arrangements may be defined by specifying any filter `Delta` in the parabolic coset poset.

## Repository Structure

```
parabolic_arrangements.sage   # Main module (v3.0)
paper.tex                     # Accompanying paper
notebooks/
  01_Cohomology_and_Asphericity.ipynb    # Cohomology ring and DGA tutorial
  02_Fundamental_Groups_and_Orbifolds.ipynb  # π₁ via complexes of groups
  03_Euler_Characteristic_and_Combinatorics.ipynb  # Euler characteristic formulas
tests/                        # Automated verification scripts
  test_stress_dga.sage        # Stress tests: Leibniz rule across types/ranks
  test_universal_coboundary.sage  # δ² = 0 verification
  test_euler.sage             # Euler characteristic cross-checks
  test_purity.sage            # Homological purity tests
  test_exotic_valid_filters.sage  # Non-standard filter validation
  test_structural_nonflag_D5.sage # D₅ non-flag example (non-commutative cells)
  test_massey_D4.sage         # Massey product explorations (D₄)
  test_C_RACG.sage            # Right-angled Coxeter groups
  test_D_B3_RACG.sage         # B₃ right-angled quotient
  test_E_B3_BSW.sage          # B₃ Barcelo–Severs–White comparison
  test_F_trianglecriterion.sage  # Triangle-free criterion verification
examples/                     # Standalone worked examples
  example_B3_selective.sage
  example_B4_4parabolic.sage
  example_borromean_triples.sage
  example_minimal_nonflag_D5.sage
  example_minimal_nonflag_D6.sage
legacy/                       # Archived earlier code
binder/                       # Binder configuration (Dockerfile)
```

## Tutorial Notebooks

Three self-contained Jupyter notebooks (SageMath kernel) serve as the computational supplement to the paper:

| Notebook | Content |
|----------|---------|
| **01 — Cohomology and Asphericity** | Cellular model, Betti numbers, cup products, Leibniz verification, asphericity criteria.  Includes a detailed case study of D₅ showing the limitations of the cubical diagonal for non-flag complexes. |
| **02 — Fundamental Groups and Orbifolds** | Permutahedral complex visualization, orbit structure, scwol construction, spanning tree, transport elements, Bridson–Haefliger group extensions, π₁ computation via three routes. |
| **03 — Euler Characteristic** | Orbit-weighted Euler characteristic formula, cross-check against alternating Betti sums. |

## Citation

If you use this software in your research, please cite:

```bibtex
@misc{cantarero2025parabolic,
  author  = {Cantarero, Jos\'e and Le\'on-Medina, Jos\'e Luis},
  title   = {The Topology of Real Parabolic Arrangements},
  year    = {2025},
  note    = {Preprint}
}
```

## Mathematical Background

A **parabolic arrangement** is a collection of parabolic cosets forming a filter in the inclusion poset of a finite Coxeter group (W, S).  The complement of such an arrangement in the Coxeter complex is homotopy equivalent to a subcomplex of the Coxeter permutahedron (Björner–Ziegler duality):

    M(A) ≃ K_A ⊂ Perm(W).

Cells of K_A are indexed by parabolic cosets wW_I, graded by |I|.  The boundary operator is:

    ∂(wW_I) = Σ_{s∈I} (-1)^{o_I(s)} Σ_{g∈W_I^{I\{s}}} (-1)^{ℓ(g)} wg W_{I\{s}}.

The cup product is defined via a diagonal approximation transported from the classical Serre diagonal on cubes to the permutahedral cells.  For cells that are *strictly commutative* (m_{ij} = 2 for all pairs in I), this reduces to an exact combinatorial formula involving longest elements of parabolic subgroups and shuffle signs.  This covers all 3-parabolic arrangements in full.

When the arrangement is W-invariant, the fundamental group is obtained from a modified Coxeter presentation where excised rank-2 strata correspond to deleted relations, and asphericity reduces to a graph-theoretic condition (triangle-free criterion).

**Key references:**
- A. Björner and G. Ziegler, *Combinatorial stratification of complex arrangements*, J. Amer. Math. Soc. **5** (1992), 105–149.
- M. Davis, *The Geometry and Topology of Coxeter Groups*, Princeton Univ. Press, 2008.
- M. R. Bridson and A. Haefliger, *Metric Spaces of Non-Positive Curvature*, Springer, 1999.
- C. Severs and J. White, *k-parabolic subspace arrangements*, European J. Combin. **33** (2012), 1002–1022.
- J.-P. Serre, *Homologie singulière des espaces fibrés. Applications*, Ann. of Math. (1951).

## License

MIT License — see [LICENSE](LICENSE).

## Acknowledgments

Developed as part of **SECIHTI CBF2023-2024-4059: "Interacciones topológico-computacionales"**.

*Portions of the implementation were developed with AI assistance and have been rigorously verified against the mathematical theory.*

---

**Version**: 3.0 | **Updated**: March 2026
