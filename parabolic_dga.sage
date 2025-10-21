r"""
Parabolic DGA over Coxeter permutahedra (chains of parabolic cosets)
====================================================================

This module constructs the (absolute or relative) chain complex of
**parabolic cosets** over a Coxeter permutahedron and equips it with:

- boundary matrices over ``ZZ`` and ``GF(2)``,
- Betti numbers (via ranks of boundary maps),
- a **transported intersection product** on chains that models
  transverse intersection of dual cells and, on homology/cohomology,
  coincides with the cup product,
- a lightweight ``self_test`` that checks ``d^2 = 0`` and the Leibniz rule.

**Design.**
- Cells are pairs ``c = (w, J)`` coming from Sage’s ``milnor_fiber_poset()``,
  where ``J`` is a tuple of simple indices and ``w`` lies in the Weyl group.
- Chain degree is ``|J|``.
- The oriented boundary over ``ZZ`` is
  :math:`\partial(w,J)=\sum_{j\notin J} (-1)^{o_J(j)} (w, J\cup\{j\})`,
  with :math:`o_J(j)=|\{i\in J: i<j\}|`.
- A parabolic arrangement is encoded by deleting all cells in a given order
  ideal ``Delta``; working relative simply means we restrict to the
  complement of ``Delta`` (the allowed cells).

**Intersection product (chain level).**
Given basis cells ``(w1, I1)`` and ``(w2, I2)``, we
- *transport* the simple directions of the second factor into the
  reference frame of the first via ``h = w1^{-1} w2``; concretely we
  require each conjugate ``h s_i h^{-1}`` (``i in I2``) to be a **global**
  simple ``s_j`` and collect the transported set ``J2``,
- enforce **transversality after transport** by requiring
  ``I1 ∩ J2 = ∅``,
- check the **same-fiber** condition relative to
  ``U = I1 ∪ J2`` by testing that the right-minimal representative of
  ``h`` modulo ``U`` is the identity (equivalently, ``h ∈ W_U``),
- output the oriented target cell ``(w1, U)`` with Koszul sign
  :math:`\epsilon(I_1,J_2)`,
- and return 0 if any step fails or if the target is deleted by ``Delta``.

Over ``GF(2)`` all signs vanish.

**Quick example** (run inside this file)::

    sage: W,P,Plist = build_W_P("A", 5)             # doctest: +ELLIPSIS
    sage: Delta = ideal_non_k_equal_A(W, Plist, k=3)
    sage: D = ParabolicDGA(W, P, Plist, Delta)
    sage: isinstance(D.betti_numbers(base_ring=GF(2)), dict)
    True

References
----------
- D. B. Fuks, *Cohomology of the braid group modulo 2*.
- J. Cantarero – J. L. León-Medina, *The Cohomology Ring of Real Parabolic Arrangements*.
"""

from sage.all import RootSystem, Integer, Matrix, vector, ZZ, GF, QQ

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
    Check a poset element has shape ``(w, J)`` with ``J`` a tuple.

    EXAMPLES::
        sage: W,P,Plist = build_W_P("A",3)      # doctest: +ELLIPSIS
        sage: all(is_cell_typed(c) for c in Plist)
        True
    """
    return isinstance(c, tuple) and len(c) == 2 and isinstance(c[1], tuple)


def chain_degree(c):
    r"""
    Chain degree ``k = |J|`` for a cell ``c = (w, J)``.

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
    Deletion ideal for the type ``A`` non-``k``-equal arrangement (``k >= 3``).

    In type :math:`A_r`, irreducible parabolics of size ``m`` correspond to
    **consecutive** index blocks of length ``m`` in the Dynkin chain.
    For non-``k``-equal we delete every ``(w,J)`` containing a consecutive block
    of length ``k-1``.

    INPUT:
        - ``W`` -- Weyl group
        - ``Plist`` -- list of poset elements
        - ``k`` -- integer >= 3

    OUTPUT:
        - ``Delta`` -- set of cells to delete (an order ideal by construction)

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
    Differential graded algebra on parabolic chains of a Coxeter permutahedron.

    Parameters
    ----------
    W : Weyl group
        ``W = RootSystem([...]).ambient_space().weyl_group()``.
    Pobj : poset
        ``W.milnor_fiber_poset()``.
    Plist : list
        Poset elements of the form ``(w, J)`` with ``J`` a tuple of simple indices.
    Delta : set
        Order ideal of cells to delete (parabolic arrangement). The working
        complex is generated by the complement of ``Delta``.

    Notes
    -----
    - Chain degree: ``k = |J|``.
    - Oriented boundary over ``ZZ`` uses the sign ``(-1)^{o_J(j)}``, where
      ``o_J(j) = |{ i in J : i < j }|``.
    - Intersection product is the **transported** intersection described
      in the module docstring.

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
        self.allowed = [c for c in self.Plist if c not in self.Delta]

        # index allowed cells by chain degree k = |J|
        self.by_k = {}
        for c in self.allowed:
            k = chain_degree(c)
            self.by_k.setdefault(k, []).append(c)
        for k in self.by_k:
            self.by_k[k].sort(key=lambda t: (str(t[0]), t[1]))

        # simple reflections and index set
        self.sref = self.W.simple_reflections()
        self.S = sorted(self.sref.keys())

        # small memoization caches
        self._cache_rmin = {}          # (w_key, subset_tuple) -> wr
        self._cache_common_fiber = {}  # (w_key, U_tuple) -> u
        self._cache_conj_U = {}        # (r_key, I_tuple, U_tuple) -> J_tuple or None
        self._cache_conj_fullS = {}    # (r_key, I_tuple) -> J_tuple or None

    # ------------------------------ orientation/signs -------------------------

    @staticmethod
    def o_J(j, J):
        r"""
        ``o_J(j) = |{ i in J : i < j }|``.

        EXAMPLES::
            sage: ParabolicDGA.o_J(3, (1,4))
            1
        """
        return sum(1 for i in J if i < j)

    @staticmethod
    def epsilon(I1, I2):
        r"""
        Koszul shuffle sign ``(-1)^{ |{(i,j) in I1×I2 : i>j}| }``.

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

    @staticmethod
    def _rank_over_field(M, base_ring):
        """Rank over QQ if ``base_ring`` is ZZ; native otherwise."""
        return Matrix(QQ, M).rank() if base_ring is ZZ else M.rank()

    # ------------------------- boundary and Betti numbers ---------------------

    def boundary_matrices(self, base_ring=ZZ):
        r"""
        Build the boundary matrices ``∂_k : C_k -> C_{k+1}`` for all k,
        over ``ZZ`` (default) or ``GF(2)``.

        Returns ``{k : Matrix}`` where *columns* index ``C_k`` and *rows*
        index ``C_{k+1}``.

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
        Betti numbers ``β_k`` computed from ranks of boundary maps
        over the chosen coefficient ring.

        Over ``ZZ`` ranks are taken over ``QQ`` (torsion ignored);
        over a field (e.g. ``GF(2)``) they are native.

        OUTPUT: ``dict {k : beta_k}``.

        EXAMPLES::
            sage: W,P,Plist = build_W_P("A",5)   # doctest: +ELLIPSIS
            sage: Delta = ideal_non_k_equal_A(W, Plist, k=3)
            sage: D = ParabolicDGA(W,P,Plist,Delta)
            sage: b = D.betti_numbers(GF(2))
            sage: 0 in b and 1 in b
            True
        """
        mats = self.boundary_matrices(base_ring=base_ring)
        ks = sorted(self.by_k.keys())
        betti = {}
        for k in ks:
            d_k = mats.get(k, Matrix(base_ring, len(self.by_k.get(k+1, [])), len(self.by_k.get(k, []))))
            rank_dk = self._rank_over_field(d_k, base_ring)
            dimCk = len(self.by_k.get(k, []))
            nullity = dimCk - rank_dk
            if k == 0:
                rank_prev = 0
            else:
                d_km1 = mats.get(k-1, Matrix(base_ring, len(self.by_k.get(k, [])), len(self.by_k.get(k-1, []))))
                rank_prev = self._rank_over_field(d_km1, base_ring)
            betti[k] = nullity - rank_prev
        return betti

    # ------------------------- memoized helpers (product) ---------------------

    @staticmethod
    def _elt_key(w):
        """Hashable key for Weyl group elements (reduced word when available)."""
        try:
            return tuple(w.reduced_word())
        except Exception:
            return ("str", str(w))

    def _right_min_rep_mod(self, w, subset):
        """
        Right-minimal representative of ``w`` modulo the parabolic generated by ``subset``.
        """
        key = (self._elt_key(w), tuple(sorted(set(subset))))
        wr = self._cache_rmin.get(key)
        if wr is not None:
            return wr
        subset_sorted = tuple(sorted(set(subset)))
        wr = w
        changed = True
        while changed:
            changed = False
            for i in subset_sorted:
                si = self.sref[i]
                w2 = wr * si
                if w2.length() < wr.length():
                    wr = w2
                    changed = True
        self._cache_rmin[key] = wr
        return wr

    def _common_fiber_rep(self, w, U):
        """Canonical representative of the right coset in ``W/W_U``."""
        key = (self._elt_key(w), tuple(sorted(set(U))))
        u = self._cache_common_fiber.get(key)
        if u is not None:
            return u
        u = self._right_min_rep_mod(w, U)
        self._cache_common_fiber[key] = u
        return u

    def _conjugate_simple_subset_via_group(self, r, I, U):
        """
        Conjugate ``{s_i : i in I}`` by ``r`` and identify the images as
        simple reflections **inside** the parabolic on ``U``.
        Return the transported index set ``J ⊂ U``, or raise ``ValueError`` if
        some conjugate is not simple in the ``U``-subsystem.
        """
        key = (self._elt_key(r), tuple(sorted(set(I))), tuple(sorted(set(U))))
        if key in self._cache_conj_U:
            val = self._cache_conj_U[key]
            if val is None:
                raise ValueError
            return val
        rinv = r.inverse()
        Uset = set(U)
        J = []
        for i in I:
            ti = r * self.sref[i] * rinv
            found = None
            for j in Uset:
                if ti == self.sref[j]:
                    found = j
                    break
            if found is None:
                self._cache_conj_U[key] = None
                raise ValueError
            J.append(found)
        val = tuple(sorted(set(J)))
        self._cache_conj_U[key] = val
        return val

    def _conjugate_simple_subset_fullS(self, r, I):
        """
        Conjugate ``{s_i : i in I}`` by ``r`` and recognize them as **global**
        simple reflections in ``S``. Return the transported index set ``J`` or
        ``None`` if some conjugate is not a global simple.
        """
        key = (self._elt_key(r), tuple(sorted(set(I))))
        if key in self._cache_conj_fullS:
            return self._cache_conj_fullS[key]
        rinv = r.inverse()
        J = []
        for i in I:
            ti = r * self.sref[i] * rinv
            found = None
            for j in self.S:
                if ti == self.sref[j]:
                    found = j
                    break
            if found is None:
                self._cache_conj_fullS[key] = None
                return None
            J.append(found)
        val = tuple(sorted(set(J)))
        self._cache_conj_fullS[key] = val
        return val

    # ------------------------- transported intersection product ----------------

    def _product_on_basis(self, c1, c2, base_ring=ZZ):
        """
        Product on basis cells implementing **transported transverse intersection**.

        Steps for cells ``c1=(w1,I1)`` and ``c2=(w2,I2)``:
          1) Transport the simple directions of the second factor via ``h=w1^{-1}w2``
             by recognizing each conjugate ``h s_i h^{-1}`` (``i in I2``) as a
             **global** simple ``s_j``; collect the transported set ``J2``.
             If some conjugate is not a global simple, the product is 0.
          2) Require transversality **after transport**: ``I1 ∩ J2 = ∅``.
          3) Let ``U = I1 ∪ J2``. Check the fiber condition ``h ∈ W_U`` by
             testing that the right-minimal representative of ``h`` modulo ``U``
             is the identity. Otherwise the product is 0.
          4) The target is ``(w1, U)`` unless it is deleted by ``Delta``.
             Over ``ZZ`` attach the sign ``epsilon(I1, J2)``; over fields use 1.

        Returns either ``(sign, target_cell)`` or ``None`` (vanishing term).
        """
        (w1, I1) = c1
        (w2, I2) = c2
        I1s = tuple(sorted(I1))
        I2s = tuple(sorted(I2))

        # transport I2 into the frame of w1 via h = w1^{-1} w2
        h = w1.inverse() * w2
        J2 = self._conjugate_simple_subset_fullS(h, I2s)
        if J2 is None:
            return None

        # transversality after transport
        if set(I1s).intersection(J2):
            return None

        # union and fiber test
        U = tuple(sorted(set(I1s).union(J2)))
        wr = self._right_min_rep_mod(h, U)
        if wr != self.W.one():
            return None

        target = (w1, U)
        if target in self.Delta:
            return None

        sgn = self.epsilon(I1s, J2) if base_ring is ZZ else 1
        return (sgn, target)

    def product(self, v1, k1, v2, k2, base_ring=ZZ):
        r"""
        Chain-level transported intersection product.

        The product is computed by expanding on the basis of ``C_{k1}`` and
        ``C_{k2}``, applying ``_product_on_basis`` pairwise, and collecting
        into ``C_{k1+k2}``. Any target basis element not present among the
        allowed cells (i.e. deleted by the ideal) is silently skipped.

        INPUT
        -----
        v1 : vector (or list coerced to a vector)
            Coordinates in the basis of ``C_{k1}``.
        k1 : int
            Degree of ``v1``.
        v2 : vector (or list coerced to a vector)
            Coordinates in the basis of ``C_{k2}``.
        k2 : int
            Degree of ``v2``.
        base_ring : ``ZZ`` or a field (e.g. ``GF(2)``)
            Coefficient ring.

        OUTPUT
        ------
        A vector in the basis of ``C_{k1+k2}`` over ``base_ring``.
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
                    continue
                out[t] += a * b * (sgn if base_ring is ZZ else 1)
        return vector(base_ring, out)

    # ------------------------- public self-test -------------------------------

    def _basis_vec(self, k, idx, base_ring):
        """Return the unit vector ``e_idx`` in ``C_k`` over ``base_ring``."""
        n = len(self.by_k.get(k, []))
        v = [base_ring(0)] * n
        v[idx] = base_ring(1)
        return vector(base_ring, v)

    def _apply_boundary(self, mats, k, v):
        """Apply ``∂_k : C_k -> C_{k+1}`` to a column vector ``v`` in ``C_k``."""
        M = mats.get(k)
        if M is None:
            return vector(v.parent().base_ring(), [])
        return M * v

    def self_test(self, trials=8, verbose=True):
        r"""
        Internal checks for the transported intersection product:

        - (1) ``d^2 = 0`` over ``ZZ`` and ``GF(2)``.
        - (2) Leibniz rule over ``ZZ`` and ``GF(2)``:
              ``d(x·y) = (dx)·y + (-1)^{|x|} x·(dy)``.

        INPUT
        -----
        trials : int
            Number of random basis-pair trials for the Leibniz check.
        verbose : bool
            Print a short summary on success.

        OUTPUT
        ------
        ``True`` if all checks pass; raises ``AssertionError`` otherwise.
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

        # (2) Leibniz rule on a few random pairs
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
            print("[ParabolicDGA.self_test] All checks passed (d^2=0, Leibniz) using the transported intersection product.")
        return True
