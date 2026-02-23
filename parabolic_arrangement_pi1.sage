# -*- coding: utf-8 -*-
r"""
Fundamental Group of Real Parabolic Arrangement Complements
==========================================================

This module computes the fundamental group \pi_1(M_A) of real parabolic
arrangement complements using the Covering Theory for Complexes of Groups
developed by Bridson and Haefliger (Ch. III.C).

It builds upon the `ParabolicArrangementCohomology` class, using its
efficient Coxeter minimization and canonical cell representations to
compute orbits, stabilizers, and the monodromy kernel under the action
of a subgroup H.
"""

from sage.all import *
from sage.graphs.digraph import DiGraph
from sage.groups.free_group import FreeGroup

class ParabolicArrangementPi1(object):
    r"""
    Bridson & Haefliger Complex of Groups calculator for Parabolic Arrangements.
    """

    def __init__(self, dga, H_gens=None):
        r"""
        INPUT:
        - ``dga`` -- An instance of ParabolicArrangementCohomology.
        - ``H_gens`` -- List of generators for the subgroup H. If None, uses
                        the trivial subgroup.
        """
        self.dga = dga
        self.W = dga.W
        
        # 1. Initialize Subgroup H
        if H_gens is None:
            self.H = self.W.subgroup([self.W.one()])
        else:
            self.H = self.W.subgroup(H_gens)
            
        # Coerce back to W so elements have WeylGroup methods (e.g. .length())
        self.H_elements = [self.W(h) for h in self.H]
        
        # Core B&H Data Structures
        self.orbits = {}         # canonical_cell -> orbit_id
        self.reps = {}           # orbit_id -> canonical_cell
        self.stabilizers = {}    # orbit_id -> list of elements in H
        
        self.Y = DiGraph()       # Quotient Scwol 1-skeleton
        self.tree_edges = set()  # Edges in the maximal tree T
        self.h_a = {}            # Transport elements: edge -> h_a in H
        
        self._build_complex_of_groups()

    def _act(self, h, cell):
        r"""
        Acts with h \in H on a canonical cell (w, J) and returns the 
        canonicalized result using the DGA's fast minimize method.
        """
        w, J = cell
        hw = h * w
        return self.dga._canonize((hw, J))

    def _build_complex_of_groups(self):
        r"""
        Executes Steps 1 and 2 of the Bridson & Haefliger algorithm.
        """
        print(f"--- Building Complex of Groups ---")
        print(f"Subgroup H size: {len(self.H_elements)}")
        
        # Step 1: Compute Orbits and Local Stabilizers (Vertices of Scwol)
        orbit_counter = 0
        for cell in self.dga.cells:
            if cell not in self.orbits:
                orbit_id = f"v_{orbit_counter}"
                self.reps[orbit_id] = cell
                self.Y.add_vertex(orbit_id)
                
                # Compute Orbit and Stabilizer
                stab = []
                for h in self.H_elements:
                    target_cell = self._act(h, cell)
                    self.orbits[target_cell] = orbit_id
                    if target_cell == cell:
                        stab.append(h)
                
                self.stabilizers[orbit_id] = stab
                orbit_counter += 1
                
        print(f"Quotient Space Y has {orbit_counter} objects (orbits).")

        # Step 2: Build Quotient Scwol 1-skeleton (Edges of Scwol)
        # In a scwol, arrows go from Vertex to Edge (i.e., cell to its face)
        # Using the DGA grading, degree k goes to degree k+1
        
        edge_counter = 0
        for k in range(self.dga.max_grade):
            for cell in self.dga.by_grade.get(k, []):
                # Only explore arrows originating from official representatives
                source_orbit = self.orbits[cell]
                if self.reps[source_orbit] != cell:
                    continue
                    
                w, J = cell
                J_set = set(J)
                
                for s_idx in self.dga.S:
                    if s_idx not in J_set:
                        J_new = tuple(sorted(J_set | {s_idx}))
                        w_new = self.dga._minimize_in_coset(w, J_new)
                        target_cell = (w_new, J_new)
                        
                        if target_cell in self.dga.map_cell_to_idx:
                            target_orbit = self.orbits[target_cell]
                            
                            # Add edge to the quotient graph Y
                            edge_id = f"e_{edge_counter}"
                            self.Y.add_edge(source_orbit, target_orbit, edge_id)
                            edge_counter += 1
                            
                            # Step 3: Compute Transport Element h_a
                            # We need h_a * target_cell = reps[target_orbit]
                            official_target = self.reps[target_orbit]
                            
                            found_ha = False
                            for h in self.H_elements:
                                if self._act(h, target_cell) == official_target:
                                    # Store the directed edge and its h_a
                                    self.h_a[(source_orbit, target_orbit, edge_id)] = h
                                    found_ha = True
                                    break
                                    
                            if not found_ha:
                                raise ValueError(f"CRITICAL: Transport element not found for {target_cell} -> {official_target}")

    # ------------------------------------------------------------------
    # Step 1: Theorem-based analysis (W-invariant case)
    # ------------------------------------------------------------------

    def removed_rank2_strata(self):
        r"""
        Returns the set of pairs (s,t) (sorted tuples of S-indices) such that
        the rank-2 stratum sigma_{st} is removed (i.e. all its cosets are in Delta).

        A stratum J=(s,t) is removed iff NO cell (w, J) survives in dga.cells.
        """
        present = set()
        for _, J in self.dga.cells:
            if len(J) == 2:
                present.add(J)

        S = self.dga.S
        removed = set()
        for i in range(len(S)):
            for j in range(i + 1, len(S)):
                J = tuple(sorted([S[i], S[j]]))
                if J not in present:
                    removed.add(J)
        return removed

    def modified_coxeter_data(self):
        r"""
        Returns a dict {(s,t): m'_{st}} where:
          m'_{st} = m_{st}  if sigma_{st} is NOT removed,
          m'_{st} = Infinity if sigma_{st} IS removed.

        This encodes the Coxeter matrix of the orbifold group G_A.
        """
        sref = self.dga._sref
        S = self.dga.S
        removed = self.removed_rank2_strata()
        data = {}
        for i in range(len(S)):
            for j in range(i + 1, len(S)):
                s, t = S[i], S[j]
                J = tuple(sorted([s, t]))
                if J in removed:
                    data[(s, t)] = Infinity
                else:
                    data[(s, t)] = (sref[s] * sref[t]).order()
        return data

    def orbifold_group_presentation(self, verbose=False):
        r"""
        Builds the orbifold group G_A as a Sage finitely presented group.

        By Theorem (explicit presentation):
          Generators: {s_i | i in S}
          Relations :  s_i^2 = 1  (involutions)
                       (s_i s_j)^{m_{ij}} = 1  if sigma_{ij} is NOT removed

        Returns (G_A, idx_map) where idx_map maps S-index -> 0-based position.
        """
        S = self.dga.S
        n = len(S)
        idx_map = {S[i]: i for i in range(n)}

        gen_names = ['s{}'.format(s) for s in S]
        F = FreeGroup(gen_names)
        fgens = F.generators()

        relations = []

        # Involution relations
        for i in range(n):
            relations.append(fgens[i] ** 2)

        # Braid relations for retained strata
        m_data = self.modified_coxeter_data()
        for (s, t), m in m_data.items():
            if m != Infinity and m >= 2:
                i, j = idx_map[s], idx_map[t]
                relations.append((fgens[i] * fgens[j]) ** m)

        G_A = F / relations

        if verbose:
            print("G_A generators:", gen_names)
            print("Removed strata:", self.removed_rank2_strata())
            print("Modified Coxeter data:", m_data)

        return G_A, idx_map

    def compute_fundamental_group(self, verbose=True):
        r"""
        Computes the fundamental group using the B&H Tree and Monodromy.
        """
        # 1. Extract Maximal Tree T in the quotient Y
        # We use Sage's breadth-first search tree
        undirected_Y = self.Y.to_undirected()
        tree_graph = undirected_Y.min_spanning_tree(algorithm='Kruskal')
        
        # Store tree edges ignoring direction
        for u, v, label in tree_graph:
            self.tree_edges.add((u, v, label))
            self.tree_edges.add((v, u, label))
            
        # 2. Identify Generators
        # Every edge NOT in the tree gives a generator
        generators = []
        monodromy_map = {} # Maps string generator to element in H
        
        if verbose:
            print(f"--- Computing Presentation ---")
        
        gen_idx = 1
        for u, v, label in self.Y.edges(sort=True):
            if (u, v, label) not in self.tree_edges:
                gen_name = f"z_{gen_idx}"
                generators.append(gen_name)
                # The monodromy evaluates this generator to its transport element
                monodromy_map[gen_name] = self.h_a[(u, v, label)]
                
                if verbose:
                    print(f"Generator {gen_name} from edge {u} -> {v} (h_a length: {monodromy_map[gen_name].length()})")
                gen_idx += 1
                
        # 3. Handle Local Groups (Torsion)
        # Note: If H acts freely (all stabilizers trivial), this part is empty.
        # For non-free actions, we would add the stabilizer generators and 
        # the psi_a and g_{a,b} relations here.
        
        free_grp = FreeGroup(generators)
        
        if verbose:
            print(f"Orbifold Group has {len(generators)} free generators.")
            
        # 4. Compute Kernel of the Monodromy
        # We define a hom from the free group to H
        # Since H is a subgroup of W, we need a permutation or matrix representation
        # GAP handles this beautifully under the hood in Sage.
        
        # (Aca entraría el código para pasar el homomorfismo a GAP/Sage y pedir el Kernel)
        # Para acciones libres (donde G_sigma = {e}):
        # El núcleo de \rho: F_n \to H es el grupo fundamental.
        
        return free_grp, monodromy_map

    # ------------------------------------------------------------------
    # Step 2: π₁ via Orbifold Extension Theorem
    # ------------------------------------------------------------------

    def compute_pi1(self, verbose=True):
        r"""
        Computes π₁(K_A) = ker(ρ: G_A → W) via coset enumeration.

        Algorithm (Bridson-Haefliger Thm + Reidemeister-Schreier):
          1. Build the orbifold group G_A (modified Coxeter group).
          2. Build the regular permutation representation of W on itself.
          3. Define ρ: G_A → W by s_i ↦ left-multiplication by s_i.
          4. Compute ker(ρ) via GAP's coset enumeration + IsomorphismFpGroup.

        Returns (K_fp, phi_gap) where:
          K_fp  -- GAP FP group isomorphic to π₁(K_A)
          phi_gap -- the GAP homomorphism G_A → W
        """
        from sage.libs.gap.libgap import libgap

        S    = self.dga.S
        sref = self.dga._sref
        n    = len(S)

        if verbose:
            print("=" * 55)
            print("Computing π₁(K_A) via Orbifold Extension Theorem")
            print("=" * 55)

        # --- Build G_A ---
        G_A, idx_map = self.orbifold_group_presentation(verbose=verbose)
        G_A_gap  = G_A.gap()
        G_A_gens = list(G_A_gap.GeneratorsOfGroup())

        if verbose:
            removed = self.removed_rank2_strata()
            print(f"Removed rank-2 strata : {removed}")
            print(f"G_A has {n} generators, {len(G_A.relations())} relations")

        # --- Regular permutation representation of W ---
        # W acts on itself by left multiplication; gives Sym(|W|) rep.
        W_elts   = list(self.W)
        w_to_idx = {w: i + 1 for i, w in enumerate(W_elts)}  # 1-indexed for GAP

        gap_perms = []
        for s_idx in S:
            s    = sref[s_idx]
            perm = [w_to_idx[s * w] for w in W_elts]
            gap_perms.append(libgap.PermList(perm))

        W_gap = libgap.Group(gap_perms)

        if verbose:
            print(f"|W| = {len(W_elts)}, permutation degree = {len(W_elts)}")

        # --- Homomorphism rho: G_A -> W ---
        phi = libgap.GroupHomomorphismByImages(G_A_gap, W_gap, G_A_gens, gap_perms)
        if phi == libgap.fail:
            raise ValueError("rho: G_A -> W is not well-defined. Check that the modified Coxeter relations are consistent with W.")

        if verbose:
            print("Monodromy ρ: G_A → W defined successfully.")

        # --- Kernel = π₁(K_A) ---
        K_gap = phi.Kernel()

        if verbose:
            print(f"Kernel computed. Running IsomorphismFpGroup...")

        K_iso = K_gap.IsomorphismFpGroup()
        K_fp  = K_iso.Range()

        if verbose:
            K_gens = list(K_fp.GeneratorsOfGroup())
            K_rels = list(K_fp.RelatorsOfFpGroup())
            print("-" * 55)
            print(f"π₁(K_A) generators ({len(K_gens)}):")
            for g in K_gens:
                print(f"  {g}")
            print(f"π₁(K_A) relations ({len(K_rels)}):")
            for r in K_rels:
                print(f"  {r} = 1")
            print("=" * 55)

        return K_fp, phi

    def pi1_abelianization(self, verbose=True):
        r"""
        Returns the abelian invariants of π₁(K_A).

        Uses the exact sequence  π₁(K_A) → G_A → W  and the
        abelianization of the kernel (= H_1 of the arrangement complement).
        """
        K_fp, _ = self.compute_pi1(verbose=verbose)
        ab = list(K_fp.AbelianInvariants())
        if verbose:
            print(f"H₁(K_A; Z) abelian invariants: {ab}")
        return ab