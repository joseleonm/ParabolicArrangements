"""
ParabolicDGA - Implementación Completa con Producto Copa Corregido

Esta es una versión lista para usar de la clase ParabolicDGA con el
producto copa zonotopal correctamente implementado.

Basado en ProductosCopa9Nov.ipynb pero con correcciones.

USO:
    load("parabolic_dga_notebook.sage")

    W, P, Plist = build_W_P("A", 4)
    Delta = ideal_non_k_equal_A(W, Plist, k=3)
    D = ParabolicDGA(W, P, Plist, Delta, vertex_key='perm', ascii_mode='perm')

    # Calcular cohomología
    H1 = D.cohomology_basis(1, ring=GF(2))

    # Producto copa (YA CORREGIDO)
    u, v = H1[0], H1[1]
    cup = D.cup_on_cochains(u, v, 1, 1, ring=GF(2))
"""

import itertools
from collections import defaultdict, deque

# ============================================================================
# FUNCIONES AUXILIARES PARA CONSTRUCCIÓN DE W, P, Plist
# ============================================================================

def build_W_P(tipo, n):
    """
    Construye el grupo de Coxeter W y el conjunto de subgrupos parabólicos.

    INPUT:
    - tipo: 'A', 'B', 'C', 'D'
    - n: rango

    OUTPUT:
    - W: CoxeterGroup
    - P: conjunto de subgrupos parabólicos (como conjuntos de índices)
    - Plist: lista de P
    """
    W = CoxeterGroup([tipo, n], implementation="permutation")

    # Índices de generadores simples (1 a n para tipo A)
    if tipo == 'A':
        indices = list(range(1, n+1))
    else:
        indices = list(range(1, n+1))

    # Todos los subconjuntos (subgrupos parabólicos)
    P = set()
    for k in range(n+1):
        for subset in itertools.combinations(indices, k):
            P.add(frozenset(subset))

    Plist = sorted([tuple(sorted(p)) for p in P])

    return W, P, Plist


def ideal_non_k_equal_A(W, Plist, k):
    """
    Ideal parabólico para tipo A: arreglo de k-igualdades.

    Un subconjunto I está en Delta si |I| >= k.

    INPUT:
    - W: grupo de Coxeter
    - Plist: lista de subgrupos parabólicos
    - k: parámetro de k-igualdad

    OUTPUT:
    - Delta: conjunto de subconjuntos en el ideal
    """
    Delta = set()
    for I in Plist:
        if len(I) >= k:
            Delta.add(frozenset(I))

    return Delta


# ============================================================================
# CLASE PRINCIPAL: ParabolicDGA
# ============================================================================

class ParabolicDGA:
    """
    Álgebra Diferencial Graduada para complejos de Salvetti parabólicos.

    Implementa el producto copa zonotopal CORRECTO que suma sobre
    TODAS las particiones.
    """

    def __init__(self, W, P, Plist, Delta, vertex_key='perm', ascii_mode='perm'):
        """
        INPUT:
        - W: CoxeterGroup
        - P: conjunto de subgrupos parabólicos
        - Plist: lista de P
        - Delta: ideal de subgrupos (building set)
        - vertex_key: 'perm' o 'reduced_word'
        - ascii_mode: 'perm' o 'reduced_word'
        """
        self.W = W
        self.P = P
        self.Plist = Plist
        self.Delta = Delta
        self.vertex_key = vertex_key
        self.ascii_mode = ascii_mode

        # Referencia a generadores simples
        self.sref = {i: W.simple_reflection(i) for i in W.index_set()}

        # Construir complejo celular
        self._build_cells()

        # Caché para estructura de producto copa
        self._cup_structure_cache = {}

    def _build_cells(self):
        """Construye las celdas del complejo de Salvetti parabólico."""
        self.cells = []
        self.by_k = defaultdict(list)

        # Celdas: (w, I) donde w ∈ W, I ∉ Delta
        for w in self.W:
            for I in self.Plist:
                if frozenset(I) not in self.Delta:
                    cell = (w, I)
                    k = len(I)
                    self.cells.append(cell)
                    self.by_k[k].append(cell)

        # Ordenar por grado
        for k in self.by_k:
            self.by_k[k].sort()

        # Crear mapeos inversos
        self._cell_to_index = {cell: i for i, cell in enumerate(self.cells)}
        self._by_k_index = {}
        for k in self.by_k:
            self._by_k_index[k] = {cell: i for i, cell in enumerate(self.by_k[k])}

    def _right_min_rep_mod(self, w, I):
        """
        Representante minimal a derecha de w en W/W_I.

        Encuentra el elemento de menor longitud en el coset wW_I.
        """
        if len(I) == 0:
            return w

        current = w
        improved = True

        while improved:
            improved = False
            for i in I:
                s_i = self.sref[i]
                candidate = current * s_i
                if candidate.length() < current.length():
                    current = candidate
                    improved = True

        return current

    def _get_w0_J(self, J):
        """
        Calcula el elemento más largo del subgrupo parabólico W_J.

        Para tipo A con J = {i}, es simplemente s_i.
        Para J más grande, busca el elemento de mayor longitud.
        """
        if len(J) == 0:
            return self.W.one()

        if len(J) == 1:
            return self.sref[J[0]]

        # Para J más grande, generar el parabólico y buscar máximo
        gens = [self.sref[i] for i in J]

        # BFS para encontrar el elemento más largo
        queue = deque([self.W.one()])
        visited = {self.W.one()}
        max_element = self.W.one()
        max_length = 0

        while queue:
            w = queue.popleft()
            if w.length() > max_length:
                max_length = w.length()
                max_element = w

            for s in gens:
                new_w = w * s
                if new_w not in visited and new_w.length() > w.length():
                    visited.add(new_w)
                    queue.append(new_w)

        return max_element

    # ========================================================================
    # BOUNDARY OPERATOR
    # ========================================================================

    def boundary_matrix(self, k, ring=GF(2)):
        """
        Matriz del operador de borde ∂_k: C_k → C_{k-1}.

        ∂(w, I ∪ {i}) = Σ_{j ∈ I ∪ {i}} (-1)^{pos(j)} (w·s_j, I ∪ {i} \ {j})
        """
        if k == 0:
            return matrix(ring, 0, len(self.by_k[0]))

        cells_k = self.by_k.get(k, [])
        cells_km1 = self.by_k.get(k-1, [])

        if not cells_k or not cells_km1:
            return matrix(ring, len(cells_km1), len(cells_k))

        # Mapeo inverso
        index_km1 = self._by_k_index[k-1]

        M = matrix(ring, len(cells_km1), len(cells_k))

        for j, (w, I) in enumerate(cells_k):
            I_set = set(I)

            for pos, i in enumerate(I):
                # Face: I \ {i}
                I_face = tuple(sorted(I_set - {i}))

                # w·s_i
                w_face = w * self.sref[i]

                # Minimizar
                w_face_min = self._right_min_rep_mod(w_face, I_face)

                cell_face = (w_face_min, I_face)

                if cell_face in index_km1:
                    idx = index_km1[cell_face]

                    # Signo
                    if ring.characteristic() == 2:
                        M[idx, j] += 1
                    else:
                        sign = (-1)**pos
                        M[idx, j] += sign

        return M

    def coboundary_matrix(self, k, ring=GF(2)):
        """Matriz del operador de coborde δ^k: C^k → C^{k+1}."""
        return self.boundary_matrix(k+1, ring=ring).transpose()

    # ========================================================================
    # COHOMOLOGÍA
    # ========================================================================

    def cohomology_basis(self, k, ring=GF(2)):
        """
        Calcula base de H^k(X; ring).

        H^k = ker(δ^k) / im(δ^{k-1})
        """
        delta_k = self.coboundary_matrix(k, ring=ring)
        delta_km1 = self.coboundary_matrix(k-1, ring=ring)

        # Kernel de δ^k
        ker = delta_k.right_kernel()

        # Imagen de δ^{k-1}
        im = delta_km1.column_space()

        # Cocientes: H^k = ker / im
        # Encontrar base de ker que proyecta a base de H^k

        # Método: encontrar complemento de im en ker
        V = ker

        # Base de H^k: levantamiento de base de ker/im
        # Usamos que ker y im son subespacios vectoriales

        # Dimensión de H^k
        dim_Hk = ker.dimension() - im.dimension()

        if dim_Hk == 0:
            return []

        # Encontrar base de ker que no está en im
        ker_basis = ker.basis()
        im_basis = im.basis()

        # Extender im_basis a base de ker
        # Los vectores adicionales forman base de H^k

        # Método simple: tomar base de ker, eliminar componente en im
        Hk_basis = []

        for v in ker_basis:
            # Verificar si v está en im
            # Si no, agregarlo a la base de H^k

            # Crear matriz con im_basis + v
            test_space = (im + span([v], ring))

            if test_space.dimension() > im.dimension():
                # v no está en im, es independiente
                Hk_basis.append(v)
                im = test_space

                if len(Hk_basis) == dim_Hk:
                    break

        return Hk_basis

    # ========================================================================
    # PRODUCTO COPA - VERSIÓN CORREGIDA
    # ========================================================================

    def _precompute_cup_structure(self, p, q):
        """
        Precomputa la estructura del producto copa C^p × C^q → C^{p+q}.

        CORRECCIÓN: Suma sobre TODAS las particiones I = J ⊔ K.

        Retorna lista de (idx_target, idx_front, idx_back, sign)
        """
        m = p + q
        cells_p = self.by_k.get(p, [])
        cells_q = self.by_k.get(q, [])
        cells_m = self.by_k.get(m, [])

        if not cells_m:
            return []

        # Mapeos inversos
        map_p = self._by_k_index[p]
        map_q = self._by_k_index[q]

        structure = []

        # Para cada celda objetivo (w, I) en C^m
        for idx_target, (w, I) in enumerate(cells_m):

            # Sumar sobre TODAS las particiones I = J ⊔ K con |J|=p, |K|=q
            for partition_indices in itertools.combinations(range(len(I)), p):
                J_list = []
                K_list = []
                shuffle_inversions = 0
                k_count = 0

                # Particionar I en J y K, contando inversiones
                for pos in range(len(I)):
                    element = I[pos]
                    if pos in partition_indices:
                        J_list.append(element)
                        shuffle_inversions += k_count
                    else:
                        K_list.append(element)
                        k_count += 1

                J = tuple(J_list)
                K = tuple(K_list)

                # Cara front: (w, J)
                w_front = self._right_min_rep_mod(w, J)
                front_cell = (w_front, J)

                if front_cell not in map_p:
                    continue

                idx_front = map_p[front_cell]

                # Cara back: (w·w₀_J, K)
                w0_J = self._get_w0_J(J)
                w_shifted = w * w0_J
                w_back = self._right_min_rep_mod(w_shifted, K)
                back_cell = (w_back, K)

                if back_cell not in map_q:
                    continue

                idx_back = map_q[back_cell]

                # Signo de shuffle
                sign = -1 if (shuffle_inversions % 2 == 1) else 1

                structure.append((idx_target, idx_front, idx_back, sign))

        return structure

    def cup_on_cochains(self, u, v, p, q, ring=GF(2), use_cache=True):
        """
        Producto copa: u ∈ C^p, v ∈ C^q → u∪v ∈ C^{p+q}

        Formula correcta (diagonal zonotopal):
        (u∪v)(w,I) = Σ_{I=J⊔K, |J|=p, |K|=q} ε(J,K) · u(w,J) · v(w·w₀_J, K)

        CORRECCIÓN: Suma sobre TODAS las particiones (no solo una).
        """
        if not hasattr(u, "parent"):
            u = vector(ring, u)
        if not hasattr(v, "parent"):
            v = vector(ring, v)

        m = p + q
        cells_m = self.by_k.get(m, [])

        if not cells_m:
            return vector(ring, [])

        # Precomputar o usar caché
        if use_cache:
            cache_key = (p, q)

            if cache_key not in self._cup_structure_cache:
                self._cup_structure_cache[cache_key] = self._precompute_cup_structure(p, q)

            structure = self._cup_structure_cache[cache_key]
        else:
            structure = self._precompute_cup_structure(p, q)

        # Calcular producto
        result = vector(ring, len(cells_m))
        char2 = (ring.characteristic() == 2)

        for (idx_target, idx_front, idx_back, sign) in structure:
            term = u[idx_front] * v[idx_back]

            if term == 0:
                continue

            if char2:
                result[idx_target] += term
            else:
                if sign == -1:
                    result[idx_target] -= term
                else:
                    result[idx_target] += term

        return result

    # ========================================================================
    # VERIFICACIÓN DE LEIBNIZ
    # ========================================================================

    def verify_leibniz(self, p, q, ring=GF(2), verbose=False):
        """
        Verifica la regla de Leibniz: δ(u∪v) = (δu)∪v + (-1)^p u∪(δv)

        para todos los pares de generadores en C^p × C^q.
        """
        cells_p = self.by_k.get(p, [])
        cells_q = self.by_k.get(q, [])

        if not cells_p or not cells_q:
            if verbose:
                print(f"Leibniz({p},{q}): No hay celdas para verificar")
            return True

        delta_p = self.coboundary_matrix(p, ring=ring)
        delta_q = self.coboundary_matrix(q, ring=ring)
        delta_m = self.coboundary_matrix(p+q, ring=ring)

        # Verificar para generadores canónicos
        failures = 0
        total = 0

        for i in range(len(cells_p)):
            for j in range(len(cells_q)):
                # Generadores
                u = vector(ring, len(cells_p))
                v = vector(ring, len(cells_q))
                u[i] = 1
                v[j] = 1

                # Lado izquierdo: δ(u∪v)
                cup_uv = self.cup_on_cochains(u, v, p, q, ring=ring)
                lhs = delta_m * cup_uv

                # Lado derecho: (δu)∪v + (-1)^p u∪(δv)
                du = delta_p * u
                dv = delta_q * v

                term1 = self.cup_on_cochains(du, v, p+1, q, ring=ring)
                term2 = self.cup_on_cochains(u, dv, p, q+1, ring=ring)

                if ring.characteristic() == 2 or p % 2 == 0:
                    rhs = term1 + term2
                else:
                    rhs = term1 - term2

                # Comparar
                if lhs != rhs:
                    failures += 1
                    if verbose:
                        print(f"  Falla en ({i},{j}): lhs={lhs}, rhs={rhs}")

                total += 1

        if verbose:
            if failures == 0:
                print(f"Leibniz({p},{q}): ✓ Verificado para {total} pares")
            else:
                print(f"Leibniz({p},{q}): ✗ {failures}/{total} fallos")

        return failures == 0

    # ========================================================================
    # UTILIDADES
    # ========================================================================

    def info(self):
        """Imprime información sobre el complejo."""
        print("=" * 80)
        print("ParabolicDGA - Información")
        print("=" * 80)
        print(f"Grupo: {self.W.cartan_type()}")
        print(f"Ideal Delta: {len(self.Delta)} subgrupos")
        print()
        print("Números de Betti:")
        for k in sorted(self.by_k.keys()):
            print(f"  C^{k}: {len(self.by_k[k])} celdas")
        print("=" * 80)


# ============================================================================
# EJEMPLO DE USO
# ============================================================================

def ejemplo_A4_k3():
    """Ejemplo: A_4, k=3"""
    print("=" * 80)
    print("EJEMPLO: A_4, k=3")
    print("=" * 80)
    print()

    # Construir sistema
    W, P, Plist = build_W_P("A", 4)
    Delta = ideal_non_k_equal_A(W, Plist, k=3)
    D = ParabolicDGA(W, P, Plist, Delta, vertex_key='perm', ascii_mode='perm')

    # Info
    D.info()
    print()

    # Cohomología
    print("Calculando cohomología...")
    H1 = D.cohomology_basis(1, ring=GF(2))
    H2 = D.cohomology_basis(2, ring=GF(2))

    print(f"dim H^1 = {len(H1)}")
    print(f"dim H^2 = {len(H2)}")
    print()

    # Producto copa
    if len(H1) >= 2:
        print("Calculando productos copa en H^1 × H^1...")
        u, v = H1[0], H1[1]
        cup = D.cup_on_cochains(u, v, 1, 1, ring=GF(2))
        print(f"u∪v: soporte de tamaño {len(cup.nonzero_positions())}")
        print()

    # Verificar Leibniz
    print("Verificando regla de Leibniz...")
    leibniz_ok = D.verify_leibniz(1, 1, ring=GF(2), verbose=True)
    print()

    if leibniz_ok:
        print("✓ Implementación correcta: Leibniz satisfecho")
    else:
        print("✗ Error: Leibniz no satisfecho")

    print()
    print("=" * 80)

    return D


# ============================================================================
# CARGAR AUTOMÁTICAMENTE
# ============================================================================

print("Cargado: parabolic_dga_notebook.sage")
print()
print("Clases disponibles:")
print("  - ParabolicDGA")
print()
print("Funciones auxiliares:")
print("  - build_W_P(tipo, n)")
print("  - ideal_non_k_equal_A(W, Plist, k)")
print("  - ejemplo_A4_k3()")
print()
print("Ejemplo rápido:")
print("  D = ejemplo_A4_k3()")
print()
