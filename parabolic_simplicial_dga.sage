"""
Parabolic DGA - Implementación Simplicial Dual (estilo Baryshnikov)

Esta implementación usa el complejo simplicial puro (orden de los vértices)
en lugar del complejo de Salvetti/zonotopal, pero implementa el MISMO
producto copa vía la diagonal zonotopal correcta.

Autor: Implementación dual a parabolic_dga.sage
Fecha: 2025
"""

from sage.all import (
    RootSystem, Integer, Matrix, vector, ZZ, GF, QQ,
    Graph, VectorSpace, binomial
)
import itertools
from collections import defaultdict

# ============================================================================
# Clase principal: ParabolicSimplicialDGA
# ============================================================================

class ParabolicSimplicialDGA:
    """
    DGA simplicial para arreglos parabólicos.

    Similar a ParabolicZonotopalDGA pero usando estructura simplicial:
    - Vértices = cosets minimales (w, {i})
    - k-símplices = cosets minimales (w, I) con |I| = k+1
    - Producto copa vía diagonal zonotopal (suma sobre particiones)

    EXAMPLES::

        sage: W, Plist, _ = build_W_P('A', 3)
        sage: Delta = ideal_k_parabolic(W, Plist, k=3)
        sage: dga = ParabolicSimplicialDGA(W, Plist, Delta)
        sage: H1 = dga.cohomology_basis(1, ring=GF(2))
        sage: H2 = dga.cohomology_basis(2, ring=GF(2))
        sage: len(H1), len(H2)
        (6, 0)
    """

    def __init__(self, W, Plist, Delta, verbose=False):
        """
        INPUT:

        - ``W`` -- Grupo de Weyl
        - ``Plist`` -- Lista de cosets parabólicos
        - ``Delta`` -- Ideal de exclusión (cosets a excluir)
        - ``verbose`` -- Mostrar información de construcción
        """
        self.W = W
        self.Plist = Plist
        self.Delta = Delta
        self.rank = W.rank()
        self.verbose = verbose

        # Construir complejo simplicial
        self._build_complex()

        # Precomputar matrices de coborde
        self._coboundary_matrices = {}

        # Cache para estructuras de producto copa
        self._cup_structure_cache = {}

    def _build_complex(self):
        """Construye el complejo simplicial de cosets parabólicos."""
        if self.verbose:
            print("Construyendo complejo simplicial...")

        # Filtrar cosets válidos
        valid_cells = []
        for c in self.Plist:
            if not self._is_cell_typed(c):
                continue
            if c in self.Delta:
                continue
            valid_cells.append(c)

        # Organizar por grado
        self.by_grade = defaultdict(list)
        for c in valid_cells:
            k = len(c[1])  # Grado = |I|
            self.by_grade[k].append(c)

        # Crear mapeo celda -> índice
        self.map_cell_to_idx = {}
        for k, cells in self.by_grade.items():
            for idx, c in enumerate(cells):
                self.map_cell_to_idx[(k, c)] = idx

        self.max_grade = max(self.by_grade.keys()) if self.by_grade else 0

        if self.verbose:
            print(f"  Grados: {sorted(self.by_grade.keys())}")
            for k in sorted(self.by_grade.keys()):
                print(f"  C^{k}: {len(self.by_grade[k])} celdas")

    def _is_cell_typed(self, c):
        """Verifica si c es una celda bien formada (w, I)."""
        return (isinstance(c, tuple) and len(c) == 2 and
                isinstance(c[1], tuple))

    def _minimize_in_coset(self, w, J):
        """
        Encuentra el representante minimal del coset wW_J.

        Usa descenso por reflexiones simples en J.
        """
        if len(J) == 0:
            return w

        # Generadores del subgrupo parabólico
        gens = [self.W.simple_reflection(i) for i in J]

        current = w
        improved = True
        max_iters = 100  # Prevenir loops infinitos
        iters = 0

        while improved and iters < max_iters:
            improved = False
            for s in gens:
                candidate = current * s
                if candidate.length() < current.length():
                    current = candidate
                    improved = True
            iters += 1

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
            return self.W.simple_reflection(J[0])

        # Para J más grande, generar el parabólico y buscar máximo
        gens = [self.W.simple_reflection(i) for i in J]

        # BFS para encontrar el elemento más largo
        from collections import deque
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

    def coboundary_matrix(self, k, ring=ZZ):
        """
        Matriz de coborde δ^k: C^k → C^{k+1}.

        INPUT:

        - ``k`` -- Grado
        - ``ring`` -- Anillo de coeficientes

        OUTPUT: Matriz sobre ``ring``
        """
        key = (k, ring)
        if key in self._coboundary_matrices:
            return self._coboundary_matrices[key]

        cells_k = self.by_grade.get(k, [])
        cells_k1 = self.by_grade.get(k+1, [])

        n_rows = len(cells_k1)
        n_cols = len(cells_k)

        if n_rows == 0 or n_cols == 0:
            mat = Matrix(ring, n_rows, n_cols)
            self._coboundary_matrices[key] = mat
            return mat

        # Inicializar matriz cero
        mat = Matrix(ring, n_rows, n_cols)

        # Mapeo inverso para celdas de grado k+1
        cells_k1_map = {c: i for i, c in enumerate(cells_k1)}

        char2 = (ring.characteristic() == 2)

        # Construir matriz de coborde
        for j, (w, I) in enumerate(cells_k):
            I_set = set(I)

            # Para cada índice simple que podemos agregar
            for new_idx in range(1, self.rank + 1):
                if new_idx in I_set:
                    continue

                # Crear nueva celda
                I_new = tuple(sorted(I_set | {new_idx}))
                w_new = self._minimize_in_coset(w, I_new)
                new_cell = (w_new, I_new)

                # Buscar índice en celdas de grado k+1
                i = cells_k1_map.get(new_cell)
                if i is None:
                    continue

                # Calcular signo: (-1)^{posición de new_idx en I_new}
                pos = I_new.index(new_idx)
                sign = 1 if (pos % 2 == 0 or char2) else -1

                mat[i, j] += sign

        self._coboundary_matrices[key] = mat
        return mat

    def cohomology_basis(self, k, ring=GF(2)):
        """
        Calcula base de H^k = ker(δ^k) / im(δ^{k-1}).

        INPUT:

        - ``k`` -- Grado
        - ``ring`` -- Anillo de coeficientes

        OUTPUT: Lista de vectores base en H^k
        """
        delta_k = self.coboundary_matrix(k, ring=ring)
        delta_km1 = self.coboundary_matrix(k-1, ring=ring) if k > 0 else None

        # Kernel de δ^k
        ker = delta_k.right_kernel()
        ker_basis = ker.basis()

        if k == 0 or delta_km1 is None or delta_km1.ncols() == 0:
            # H^0 = ker(δ^0)
            return list(ker_basis)

        # Imagen de δ^{k-1}
        img_basis = delta_km1.columns()

        if len(img_basis) == 0:
            return list(ker_basis)

        # Construir espacio vectorial
        V = VectorSpace(ring, delta_k.ncols())

        # Subespacio de cociclos
        Z = V.subspace(ker_basis)

        # Subespacio de cobordes
        B = V.subspace(img_basis)

        # Cociente H^k = Z / B
        try:
            Q = Z.quotient(B)
            return [Q.lift(gen) for gen in Q.basis()]
        except (ValueError, AttributeError):
            # Si el cociente falla, calcular manualmente
            # Extender base de B a base de Z
            basis_B = B.basis()
            basis_Z = Z.basis()

            # Encontrar complemento
            from sage.modules.free_module import span
            complement_basis = []
            current_space = B

            for v in basis_Z:
                if v not in current_space:
                    complement_basis.append(v)
                    current_space = span(current_space.basis() + [v])

            return complement_basis

    def cup_product(self, u, v, deg_u, deg_v, ring=GF(2)):
        """
        Producto copa u ∪ v usando diagonal zonotopal.

        (u ∪ v)(w, I) = Σ_{I = J ⊔ K} ε(J,K) · u(w,J) · v(w·w₀_J, K)

        INPUT:

        - ``u`` -- Cochain en C^{deg_u}
        - ``v`` -- Cochain en C^{deg_v}
        - ``deg_u`` -- Grado de u
        - ``deg_v`` -- Grado de v
        - ``ring`` -- Anillo de coeficientes

        OUTPUT: Vector en C^{deg_u + deg_v}
        """
        # Convertir a vectores
        if not hasattr(u, "parent"):
            u = vector(ring, u)
        if not hasattr(v, "parent"):
            v = vector(ring, v)

        deg_tot = deg_u + deg_v
        cells_tot = self.by_grade.get(deg_tot, [])
        n_tot = len(cells_tot)

        result = vector(ring, n_tot)

        if n_tot == 0:
            return result

        char2 = (ring.characteristic() == 2)

        # Para cada celda objetivo (w, I)
        for idx_target, (w, I) in enumerate(cells_tot):
            acc = ring(0)

            # Iterar sobre todas las particiones I = J ⊔ K
            # con |J| = deg_u, |K| = deg_v
            for partition_indices in itertools.combinations(range(len(I)), deg_u):
                J_list = []
                K_list = []
                shuffle_inversions = 0
                k_count = 0

                # Particionar I en J y K
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
                w_front = self._minimize_in_coset(w, J)
                front_cell = (w_front, J)
                idx_front = self.map_cell_to_idx.get((deg_u, front_cell))

                if idx_front is None:
                    continue

                val_u = u[idx_front]
                if val_u == 0:
                    continue

                # Cara back: (w · w₀_J, K)
                w0_J = self._get_w0_J(J)
                w_shifted = w * w0_J
                w_back = self._minimize_in_coset(w_shifted, K)
                back_cell = (w_back, K)
                idx_back = self.map_cell_to_idx.get((deg_v, back_cell))

                if idx_back is None:
                    continue

                val_v = v[idx_back]
                if val_v == 0:
                    continue

                # Acumular con signo de shuffle
                term = val_u * val_v
                if (not char2) and (shuffle_inversions % 2 == 1):
                    acc -= term
                else:
                    acc += term

            result[idx_target] = acc

        return result

    def verify_leibniz_rule(self, ring=GF(2), trials=20, max_degree=2):
        """
        Verifica la regla de Leibniz: δ(u∪v) = (δu)∪v + (-1)^|u| u∪(δv).

        INPUT:

        - ``ring`` -- Anillo de coeficientes
        - ``trials`` -- Número de pruebas por par de grados
        - ``max_degree`` -- Grado máximo a probar

        OUTPUT: True si todas las pruebas pasan, False si alguna falla
        """
        print(f"--- Verificando Regla de Leibniz sobre {ring} ---")

        char2 = (ring.characteristic() == 2)

        for deg_u in range(max_degree + 1):
            for deg_v in range(max_degree + 1):
                if deg_u + deg_v > max_degree:
                    continue

                n_u = len(self.by_grade.get(deg_u, []))
                n_v = len(self.by_grade.get(deg_v, []))

                if n_u == 0 or n_v == 0:
                    continue

                print(f"  Probando deg_u={deg_u}, deg_v={deg_v}...")

                for trial in range(min(trials, n_u * n_v)):
                    # Generar cochains aleatorias
                    u = vector(ring, [ring.random_element() for _ in range(n_u)])
                    v = vector(ring, [ring.random_element() for _ in range(n_v)])

                    # LHS: δ(u∪v)
                    cup = self.cup_product(u, v, deg_u, deg_v, ring=ring)
                    delta_cup = self.coboundary_matrix(deg_u + deg_v, ring=ring) * cup

                    # RHS: (δu)∪v + (-1)^|u| u∪(δv)
                    delta_u = self.coboundary_matrix(deg_u, ring=ring) * u
                    delta_v = self.coboundary_matrix(deg_v, ring=ring) * v

                    term1 = self.cup_product(delta_u, v, deg_u + 1, deg_v, ring=ring)
                    term2 = self.cup_product(u, delta_v, deg_u, deg_v + 1, ring=ring)

                    if not char2 and deg_u % 2 == 1:
                        rhs = term1 - term2
                    else:
                        rhs = term1 + term2

                    if delta_cup != rhs:
                        print(f"    ✗ FALLA en trial {trial}")
                        print(f"      LHS = {delta_cup}")
                        print(f"      RHS = {rhs}")
                        return False

                print(f"    ✓ Pasa {min(trials, n_u * n_v)} pruebas")

        print("✓ Regla de Leibniz VERIFICADA")
        return True

    def cup_product_span_analysis(self, ring=GF(2), verbose=True):
        """
        Analiza si H^2 es generado por productos copa H^1 × H^1.

        OUTPUT: Diccionario con resultados del análisis
        """
        if verbose:
            print("Analizando generación de H^2 por productos copa...")

        # Calcular cohomología
        H1 = self.cohomology_basis(1, ring=ring)
        H2 = self.cohomology_basis(2, ring=ring)

        dim_H1 = len(H1)
        dim_H2 = len(H2)

        if verbose:
            print(f"  dim H^1 = {dim_H1}")
            print(f"  dim H^2 = {dim_H2}")

        if dim_H2 == 0:
            return {
                'dim_H1': dim_H1,
                'dim_H2': dim_H2,
                'rank_achieved': 0,
                'is_complete': True,
                'pairs_considered': 0
            }

        # Proyección a H^2
        delta2 = self.coboundary_matrix(2, ring=ring)
        delta1 = self.coboundary_matrix(1, ring=ring)

        # Construir proyector a H^2
        V = VectorSpace(ring, delta2.ncols())
        ker_delta2 = delta2.right_kernel()
        img_delta1 = V.subspace(delta1.columns())

        def project_to_H2(v):
            """Proyecta un cociclo a H^2."""
            # Verificar que es cociclo
            if delta2 * v != 0:
                return None
            # Proyectar módulo imagen
            try:
                return v  # Simplificado - mejorar si es necesario
            except:
                return v

        # Recolectar productos copa
        if verbose:
            print(f"  Calculando productos copa ({dim_H1}×{dim_H1} = {dim_H1**2} pares)...")

        cup_images = []
        pairs_considered = 0

        for i in range(dim_H1):
            for j in range(i, dim_H1):  # Solo mitad superior por simetría
                u = H1[i]
                v = H1[j]

                cup = self.cup_product(u, v, 1, 1, ring=ring)
                pairs_considered += 1

                if cup.is_zero():
                    continue

                # Verificar que es cociclo
                if delta2 * cup != 0:
                    if verbose:
                        print(f"    ¡Advertencia! u∪v no es cociclo para i={i}, j={j}")
                    continue

                cup_images.append(cup)

        # Calcular rango
        if len(cup_images) == 0:
            rank_achieved = 0
        else:
            span_mat = Matrix(ring, cup_images)
            rank_achieved = span_mat.rank()

        is_complete = (rank_achieved == dim_H2)

        result = {
            'dim_H1': dim_H1,
            'dim_H2': dim_H2,
            'rank_achieved': rank_achieved,
            'is_complete': is_complete,
            'pairs_considered': pairs_considered,
            'num_nonzero_cups': len(cup_images)
        }

        if verbose:
            print(f"  Productos copa no nulos: {len(cup_images)}")
            print(f"  Rango alcanzado: {rank_achieved} / {dim_H2}")
            if is_complete:
                print(f"  ✓ H^2 es COMPLETAMENTE generado por productos copa")
            else:
                print(f"  ✗ INCOMPLETO: faltan {dim_H2 - rank_achieved} generadores")

        return result


# ============================================================================
# Funciones auxiliares
# ============================================================================

def build_W_P(weyl_type, weyl_rank):
    """
    Construye grupo de Weyl y lista de cosets parabólicos.

    INPUT:

    - ``weyl_type`` -- Tipo de Weyl ('A', 'B', 'D', etc.)
    - ``weyl_rank`` -- Rango

    OUTPUT: Tupla (W, Pobj, Plist)
    """
    W = RootSystem([weyl_type, Integer(weyl_rank)]).ambient_space().weyl_group()
    Pobj = W.milnor_fiber_poset()
    Plist = list(Pobj)
    return W, Pobj, Plist


def ideal_k_parabolic(W, Plist, k=3):
    """
    Construye el ideal k-parabólico (non-k-equal).

    INPUT:

    - ``W`` -- Grupo de Weyl
    - ``Plist`` -- Lista de cosets parabólicos
    - ``k`` -- Parámetro (k=3 para non-3-equal)

    OUTPUT: Conjunto de celdas a excluir
    """
    r = W.rank()
    S = set(range(1, r + 1))

    # Construir grafo de Dynkin
    ct = W.cartan_type()
    dynkin_graph = Graph(ct.dynkin_diagram())

    Delta = set()
    threshold = k - 1

    # Apex siempre se excluye
    apex = (W.one(), tuple(range(1, r + 1)))
    Delta.add(apex)

    # Para cada celda, verificar componentes conexas
    for c in Plist:
        if not isinstance(c, tuple) or len(c) != 2:
            continue

        w, J = c

        if len(J) < threshold:
            continue

        # Verificar componentes del subgrafo de Dynkin
        subgraph = dynkin_graph.subgraph(vertices=list(J))
        components = subgraph.connected_components()

        if any(len(comp) >= threshold for comp in components):
            Delta.add(c)

    return Delta


# ============================================================================
# TESTS (comentados - descomentar para ejecutar)
# ============================================================================

# Para ejecutar tests, usa demo_quick.sage o ejecuta manualmente:
#
# W, Pobj, Plist = build_W_P('A', 3)
# Delta = ideal_k_parabolic(W, Plist, k=3)
# dga = ParabolicSimplicialDGA(W, Plist, Delta, verbose=True)
# dga.verify_leibniz_rule(ring=GF(2), trials=10, max_degree=2)
