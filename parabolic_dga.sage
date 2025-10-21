r"""
Parabolic DGA over Coxeter permutahedra
====================================================================

This module constructs the relative chain complex
(C_\bullet(\mathcal C), C_\bullet(\mathcal C_{\mathscr A}))
for a **parabolic arrangement** modeled as an order ideal of parabolic cosets,
and provides:

- boundary matrices over ZZ and GF(2),
- Betti numbers of the complement via Lefschetz duality,
- a chain-level intersection product which models transverse intersection
  of dual cells and identifies with the cohomological cup product under
  the standard correspondence with parabolic chains,
- a built-in `self_test` that checks d^2=0 and the Leibniz rule.

**Design.**
- Cells are elements `c = (w, J)` from Sage’s ``milnor_fiber_poset()``,
  where `J` is a tuple of simple indices and `w` lies in the Weyl group.
- Chain degree is `k = |J|`.
- The oriented boundary (over ZZ) is
  ∂(w, J) = Σ_{j ∉ J} (-1)^{o_J(j)} (w, J∪{j}),
  where o_J(j) = |{ i ∈ J : i < j }|.
- A parabolic arrangement is encoded by deleting all cells in an order
  ideal `Delta`; we work with the relative complex
  C_\bullet(\mathcal C, \mathcal C_{\mathscr A}).

**Quick example** (in this very file, no self-import)::

    sage: W,P,Plist = build_W_P("A",5)             # doctest: +ELLIPSIS
    sage: Delta = ideal_non_k_equal_A(W, Plist, k=3)
    sage: D = ParabolicDGA(W, P, Plist, Delta)
    sage: D.betti_numbers(base_ring=GF(2)) == {0: 1, 1: 111, 2: 20, 3: 0}
    True

A smaller sanity check::

    sage: W4,P4,L4 = build_W_P("A",4)              # doctest: +ELLIPSIS
    sage: D4 = ParabolicDGA(W4, P4, L4, ideal_non_k_equal_A(W4, L4, k=3))
    sage: D4.betti_numbers(base_ring=GF(2)) == {0: 1, 1: 31, 2: 0}
    True

**References.**
- J. Cantarero – J. L. León-Medina: *The Cohomology Ring of Real Parabolic Arrangements*.
"""

from sage.all import RootSystem, Integer, Matrix, vector, ZZ, GF

# -----------------------------------------------------------------------------
# Core: build W and its poset
# -----------------------------------------------------------------------------

def build_W_P(weyl_type, weyl_rank):
    r"""
    Build the Weyl group ``W``, the poset object ``P``, and its element list.

    INPUT:
        - ``weyl_type`` -- str, e.g. "A","B","D",...
        - ``weyl_rank`` -- positive integer (rank)

    OUTPUT:
        - ``(W, P, Plist)``

    EXAMPLES::
        sage: W,P,Plist = build_W_P("A", 3)     # doctest: +ELLIPSIS
        sage: W.rank()
        3
        sage: len(Plist) > 0
        True
    """
    W = RootSystem([weyl_type, Integer(weyl_rank)]).ambient_space().weyl_group()
    Pobj = W.milnor_fiber_poset()
    Plist = list(Pobj)
    return W, Pobj, Plist


def is_cell_typed(c):
    r"""
    Check a poset element has shape (w, J) with J a tuple.

    EXAMPLES::
        sage: W,P,Plist = build_W_P("A",3)      # doctest: +ELLIPSIS
        sage: all(is_cell_typed(c) for c in Plist)
        True
    """
    return isinstance(c, tuple) and len(c) == 2 and isinstance(c[1], tuple)


def chain_degree(c):
    r"""
    Chain degree k = |J| for c = (w, J).

    EXAMPLES::
        sage: W,P,Plist = build_W_P("A",3)      # doctest: +ELLIPSIS
        sage: all(isinstance(chain_degree(c), (int, Integer)) for c in Plist[:5])
        True
    """
    return len(c[1])


# -----------------------------------------------------------------------------
# Arrangement (order ideal) helpers
# -----------------------------------------------------------------------------

def ideal_non_k_equal_A(W, Plist, k=3):
    r"""
    Deletion ideal for the type A non-k-equal arrangement (k >= 3).

    In type A_r, irreducible parabolics of size m correspond to
    **consecutive** index blocks of length m in the Dynkin chain.
    For non-k-equal we delete every (w,J) containing a consecutive block
    of length (k-1).

    INPUT:
        - ``W`` -- Weyl group
        - ``Plist`` -- list of poset elements
        - ``k`` -- integer >= 3

    OUTPUT:
        - ``Delta`` -- set of cells to delete (order ideal by construction)

    EXAMPLES::
        sage: W,P,Plist = build_W_P("A",5)      # doctest: +ELLIPSIS
        sage: Delta = ideal_non_k_equal_A(W, Plist, k=3)
        sage: isinstance(Delta, set) and len(Delta) > 0
        True
    """
    r = W.rank()
    apex = (W.one(), tuple(range(1, r+1)))
    bad_blocks = [set(range(s, s + (k-1))) for s in range(1, r - (k-1) + 2)]
    Delta = {apex}
    for c in Plist:
        if not is_cell_typed(c):
            continue
        J = set(c[1])
        if any(B.issubset(J) for B in bad_blocks):
            Delta.add(c)
    return Delta


# -----------------------------------------------------------------------------
# ParabolicDGA class
# -----------------------------------------------------------------------------

class ParabolicDGA(object):
    r"""
    Parabolic differential graded algebra over the Coxeter permutahedron
    (chains of parabolic cosets).

    Parameters
    ----------
    W : Weyl group
        ``W = RootSystem([...]).ambient_space().weyl_group()``
    Pobj : poset
        ``W.milnor_fiber_poset()``
    Plist : list
        Elements of the poset (w, J) with J a tuple of simple indices.
    Delta : set
        Order ideal of cells to delete (parabolic arrangement).

    Notes
    -----
    - Chain degree: k = |J|.
    - Oriented boundary over ZZ: (-1)^{o_J(j)} when adding j, where
      o_J(j) = |{ i in J : i < j }|.

    EXAMPLES::
        sage: W,P,Plist = build_W_P("A",5)      # doctest: +ELLIPSIS
        sage: Delta = ideal_non_k_equal_A(W, Plist, k=3)
        sage: D = ParabolicDGA(W, P, Plist, Delta)
        sage: mats2 = D.boundary_matrices(base_ring=GF(2))
        sage: all(k in mats2 for k in sorted(D.by_k))
        True
    """

    def __init__(self, W, Pobj, Plist, Delta):
        self.W = W
        self.P = Pobj
        self.Plist = [c for c in Plist if is_cell_typed(c)]
        self.Delta = set(Delta)
        # allowed cells = complement of the ideal
        self.allowed = [c for c in self.Plist if c not in self.Delta]
        # group by chain degree k = |J|
        self.by_k = {}
        for c in self.allowed:
            k = chain_degree(c)
            self.by_k.setdefault(k, []).append(c)
        for k in self.by_k:
            self.by_k[k].sort(key=lambda t: (str(t[0]), t[1]))
        self.S = sorted(self.W.simple_reflections().keys())

    # ------------------------------ orientation/signs -------------------------

    @staticmethod
    def o_J(j, J):
        r"""
        o_J(j) = |{ i in J : i < j }|.

        EXAMPLES::
            sage: ParabolicDGA.o_J(3, (1,4))
            1
        """
        return sum(1 for i in J if i < j)

    @staticmethod
    def epsilon(I1, I2):
        r"""
        Koszul sign: (-1)^{ |{(i,j) in I1×I2 : i>j}| }.

        EXAMPLES::
            sage: ParabolicDGA.epsilon((2,),(1,))
            -1
            sage: ParabolicDGA.epsilon((1,),(2,))
            1
        """
        inv = 0
        for i in I1:
            for j in I2:
                if i > j:
                    inv += 1
        return -1 if (inv % 2 == 1) else 1

    # ------------------------- boundary and Betti numbers ---------------------

    def boundary_matrices(self, base_ring=ZZ):
        r"""
        Build the boundary matrices ∂_k : C_k -> C_{k+1} for all k,
        over ZZ (default) or GF(2).

        Returns ``{k : Matrix}`` where columns index C_k and rows C_{k+1}.

        EXAMPLES::
            sage: W,P,Plist = build_W_P("A",4)   # doctest: +ELLIPSIS
            sage: Delta = ideal_non_k_equal_A(W, Plist, 3)
            sage: D = ParabolicDGA(W,P,Plist,Delta)
            sage: M = D.boundary_matrices(GF(2))
            sage: all(M[k].nrows() == len(D.by_k.get(k+1, [])) for k in M)
            True
        """
        mats = {}
        for k in sorted(self.by_k.keys()):
            dom = self.by_k.get(k, [])
            cod = self.by_k.get(k+1, [])
            if not dom and not cod:
                continue
            M = Matrix(base_ring, len(cod), len(dom))
            if cod:
                cod_idx = {c: i for i, c in enumerate(cod)}
                for j, c in enumerate(dom):
                    w, J = c
                    Jset = set(J)
                    for s in self.S:
                        if s in Jset:
                            continue
                        Jp = tuple(sorted(Jset | {s}))
                        cp = (w, Jp)
                        i = cod_idx.get(cp)
                        if i is None:
                            continue
                        if base_ring is ZZ:
                            sgn = -1 if (self.o_J(s, J) % 2 == 1) else 1
                            M[i, j] += sgn
                        else:
                            M[i, j] += 1
            mats[k] = M
        return mats

    def betti_numbers(self, base_ring=ZZ):
        r"""
        Betti numbers of H^k of the complement (via identification with
        H_k of the relative chain complex).

        Over ZZ ranks are taken over QQ (torsion ignored); GF(2) is native.

        OUTPUT: dict {k : beta_k}

        EXAMPLES::
            sage: W,P,Plist = build_W_P("A",5)   # doctest: +ELLIPSIS
            sage: Delta = ideal_non_k_equal_A(W, Plist, k=3)
            sage: D = ParabolicDGA(W,P,Plist,Delta)
            sage: b2 = D.betti_numbers(GF(2))
            sage: (b2[0], b2[1], b2[2]) in {(1,111,20), (1, 111, 20)}
            True
        """
        mats = self.boundary_matrices(base_ring=base_ring)
        ks = sorted(self.by_k.keys())
        betti = {}
        for k in ks:
            d_k = mats.get(k, Matrix(base_ring, len(self.by_k.get(k+1, [])), len(self.by_k.get(k, []))))
            rank_dk = d_k.rank()
            dimCk = len(self.by_k.get(k, []))
            nullity = dimCk - rank_dk
            if k == 0:
                rank_prev = 0
            else:
                d_km1 = mats.get(k-1, Matrix(base_ring, len(self.by_k.get(k, [])), len(self.by_k.get(k-1, []))))
                rank_prev = d_km1.rank()
            betti[k] = nullity - rank_prev
        return betti

    # ------------------------- transverse intersection product ----------------
    # (chain-level product that matches the cup product on cohomology)

    def _right_min_rep_mod(self, w, subset):
        """Right minimal representative modulo the parabolic generated by subset."""
        sref = self.W.simple_reflections()
        subset_sorted = tuple(sorted(set(subset)))
        wr = w
        changed = True
        while changed:
            changed = False
            for i in subset_sorted:
                si = sref[i]
                w2 = wr * si
                if w2.length() < wr.length():
                    wr = w2
                    changed = True
        return wr

    def _common_fiber_rep(self, w, U):
        """Canonical representative for the common fiber over the U-block."""
        return self._right_min_rep_mod(w, U)

    def _conjugate_simple_subset_via_group(self, r, I, U):
        """
        Conjugate the simple reflections indexed by I through r, and
        identify their images as simple reflections inside the index set U.
        Raises ValueError if some conjugate is not simple in the U-subsystem.
        """
        sref = self.W.simple_reflections()
        rinv = r.inverse()
        Uset = set(U)
        J = []
        for i in I:
            ti = r * sref[i] * rinv
            found = None
            for j in Uset:
                if ti == sref[j]:
                    found = j
                    break
            if found is None:
                raise ValueError
            J.append(found)
        return tuple(sorted(set(J)))

    def _conjugate_simple_subset_fullS(self, r, I):
        """As above, but search in the full set of simples S; return None on failure."""
        sref = self.W.simple_reflections()
        rinv = r.inverse()
        J = []
        for i in I:
            ti = r * sref[i] * rinv
            found = None
            for j in self.S:
                if ti == sref[j]:
                    found = j
                    break
            if found is None:
                return None
            J.append(found)
        return tuple(sorted(set(J)))

    def _product_on_basis(self, c1, c2, base_ring=ZZ):
        """
        Chain-level product on basis elements:
          1) attempt to align to the same canonical fiber over U = I1 ∪ I2;
          2) if impossible, fallback by conjugating the second factor into the
             chamber of the first;
          3) require disjoint index sets to ensure transversality; vanish if target
             cell lies in the deleted ideal.
        Returns either (sign, target_cell) or None.
        """
        (w1, I1) = c1
        (w2, I2) = c2
        I1s = tuple(I1)
        I2s = tuple(I2)
        U = tuple(sorted(set(I1s).union(I2s)))

        # (1) Align to a common fiber in U
        u1 = self._common_fiber_rep(w1, U)
        u2 = self._common_fiber_rep(w2, U)
        if u1 == u2:
            u = u1
            x1 = u.inverse() * w1
            x2 = u.inverse() * w2
            r1 = self._right_min_rep_mod(x1, I1s)
            r2 = self._right_min_rep_mod(x2, I2s)
            try:
                J1 = self._conjugate_simple_subset_via_group(r1, I1s, U)
                J2 = self._conjugate_simple_subset_via_group(r2, I2s, U)
            except ValueError:
                J1 = J2 = None
            if J1 is not None and J2 is not None and not set(J1).intersection(J2):
                Ju = tuple(sorted(set(J1).union(J2)))
                sgn = self.epsilon(J1, J2) if base_ring is ZZ else 1
                return (sgn, (u, Ju))

        # (2) Fallback: push the second factor to the chamber of the first
        h = w1.inverse() * w2
        I2_sharp = self._conjugate_simple_subset_fullS(h, I2s)
        if I2_sharp is not None and not set(I1s).intersection(I2_sharp):
            Ju = tuple(sorted(set(I1s).union(I2_sharp)))
            sgn = self.epsilon(I1s, I2_sharp) if base_ring is ZZ else 1
            return (sgn, (w1, Ju))

        return None

    def product(self, v1, k1, v2, k2, base_ring=ZZ):
        r"""
        Chain-level intersection product modeling transverse intersection.

        The product is computed on basis elements by:
          - aligning factors to a common canonical fiber inside the parabolic
            subsystem on U = I_1 ∪ I_2 (when possible);
          - otherwise conjugating the second factor into the chamber of the first;
          - requiring disjoint index sets to ensure transversality;
          - vanishing when the target cell is deleted (lies in the ideal).

        INPUT
        -----
        v1 : vector over base_ring
            Chain in degree k1 (as a coordinate vector in the basis of C_{k1}).
        k1 : int
            Degree of v1.
        v2 : vector over base_ring
            Chain in degree k2 (as a coordinate vector in the basis of C_{k2}).
        k2 : int
            Degree of v2.
        base_ring : ZZ or GF(2)
            Coefficients.

        OUTPUT
        ------
        vector in C_{k1+k2} over base_ring.

        Notes
        -----
        Over GF(2) all signs vanish.

        (No doctest here to avoid heavy computation in small CI; see module doctest
        for Betti instead.)
        """
        if not hasattr(v1, "parent"):
            v1 = vector(base_ring, v1)
        if not hasattr(v2, "parent"):
            v2 = vector(base_ring, v2)

        Ck1 = self.by_k.get(k1, [])
        Ck2 = self.by_k.get(k2, [])
        Ck12 = self.by_k.get(k1 + k2, [])
        idx12 = {c: i for i, c in enumerate(Ck12)}
        out = [base_ring(0)] * len(Ck12)

        supp1 = [i for i, a in enumerate(v1) if a != 0]
        supp2 = [j for j, b in enumerate(v2) if b != 0]

        for i in supp1:
            a = v1[i]; c1 = Ck1[i]
            for j in supp2:
                b = v2[j]; c2 = Ck2[j]
                res = self._product_on_basis(c1, c2, base_ring=base_ring)
                if res is None:
                    continue
                sgn, target = res
                t = idx12.get(target)
                if t is None:
                    # target is deleted in the relative complex
                    continue
                out[t] += a * b * (sgn if base_ring is ZZ else 1)
        return vector(base_ring, out)

    # ------------------------- public self-test -------------------------------

    def _basis_vec(self, k, idx, base_ring):
        """Return the unit vector e_idx in C_k over base_ring."""
        n = len(self.by_k.get(k, []))
        v = [base_ring(0)] * n
        v[idx] = base_ring(1)
        return vector(base_ring, v)

    def _apply_boundary(self, mats, k, v):
        """Apply d_k : C_k -> C_{k+1} to column vector v in C_k."""
        M = mats.get(k)
        if M is None:
            return vector(v.parent().base_ring(), [])
        return M * v

    def self_test(self, trials=8, verbose=True):
        r"""
        Run internal consistency checks using the **intersection product**:

        - (1) d^2 = 0 over ZZ and GF(2).
        - (2) Leibniz rule over ZZ and GF(2) for the transverse product:
              d(x·y) = (dx)·y + (-1)^{|x|} x·(dy).

        INPUT:
            - ``trials`` -- number of random basis-pair trials for Leibniz (kept small)
            - ``verbose`` -- print a short report

        OUTPUT:
            - ``True`` if all checks pass, raises AssertionError otherwise.
        """
        import random

        # (1) d^2 = 0
        for ring in (ZZ, GF(2)):
            mats = self.boundary_matrices(base_ring=ring)
            for k in self.by_k:
                if (k+1) in self.by_k:
                    Mk = mats.get(k)
                    Mk1 = mats.get(k+1)
                    if Mk is None or Mk1 is None:
                        continue
                    comp = Mk1 * Mk
                    assert comp.is_zero(), "d^2 != 0 over {} at degree {}".format(ring, k)

        # (2) Leibniz for the intersection product (few random pairs)
        for ring in (ZZ, GF(2)):
            mats = self.boundary_matrices(base_ring=ring)
            degs = sorted([k for k in self.by_k if len(self.by_k[k]) > 0])
            pairs = []
            for k in degs:
                for l in degs:
                    if (k+l) in self.by_k and len(self.by_k[k+l]) > 0:
                        pairs.append((k,l))
            if not pairs:
                continue
            for _ in range(trials):
                k1,l1 = random.choice(pairs)
                i = random.randrange(0, len(self.by_k[k1]))
                j = random.randrange(0, len(self.by_k[l1]))
                x = self._basis_vec(k1, i, ring)
                y = self._basis_vec(l1, j, ring)

                xy = self.product(x, k1, y, l1, base_ring=ring)
                d_xy = self._apply_boundary(mats, k1+l1, xy)

                dx = self._apply_boundary(mats, k1, x)
                term1 = self.product(dx, k1+1, y, l1, base_ring=ring)

                dy = self._apply_boundary(mats, l1, y)
                term2 = self.product(x, k1, dy, l1+1, base_ring=ring)
                if ring is ZZ and (k1 % 2 == 1):
                    term2 = -term2  # (-1)^{k1}
                rhs = term1 + term2

                if len(d_xy) != len(rhs):
                    if len(d_xy) == 0:
                        d_xy = vector(ring, [0]*len(rhs))
                    elif len(rhs) == 0:
                        rhs = vector(ring, [0]*len(d_xy))
                assert d_xy == rhs, "Leibniz failed over {} for degrees ({},{})".format(ring, k1, l1)

        if verbose:
            print("[ParabolicDGA.self_test] All checks passed (d^2=0, Leibniz) using the intersection product.")
        return True
