# Índice de Archivos del Proyecto

## 🚀 Archivos Principales (Usar Estos)

### Implementación Principal
- **`parabolic_dga_notebook.sage`** ⭐ **USA ESTE** - Implementación completa con producto copa corregido
- **`parabolic_dga.sage`** - Implementación optimizada de referencia (ParabolicZonotopalDGA)

### Análisis y Herramientas
- **`analisis_cancelaciones.sage`** - Análisis de cancelaciones en suma de particiones
- **`intersection_dual_optimized.sage`** - Producto de intersección dual optimizado
- **`parabolic_simplicial_dga.sage`** - Versión dual completa (lenta, pedagógica)

## 📚 Documentación del Proyecto

### En Raíz
- **`README.md`** - Documentación principal del proyecto
- **`CHANGELOG.md`** - Historial de cambios

### En `.docs_session/` (Carpeta Oculta)
Documentación de la sesión de debugging y corrección:

#### Guías de Uso
- `USO_RAPIDO.md` - Guía rápida para empezar
- `QUICKSTART.md` - QuickStart original

#### Documentación Técnica
- `SOLUCION_PRODUCTO_COPA.md` - Corrección detallada del producto copa
- `TEORIA_DUALIZACION.md` - Teoría matemática completa
- `ROADMAP_TEORIA_INTERSECCION.md` - Plan de implementación del dual
- `RESUMEN_SESION_FINAL.md` - Resumen de lo logrado

#### Código de Sesión
- `ProductoCopa_Corregido.sage` - Parche obsoleto (usar parabolic_dga_notebook.sage)
- `test_producto_corregido.sage` - Tests del parche
- `demo_quick.sage`, `test_A5_*.sage` - Tests varios
- `dualidad_verificacion.sage` - Prototipo de verificación

## 🔧 Uso Rápido

```sage
# Cargar implementación completa
load("parabolic_dga_notebook.sage")

# Ejemplo automático
D = ejemplo_A4_k3()

# O personalizado
W, P, Plist = build_W_P("A", 4)
Delta = ideal_non_k_equal_A(W, Plist, k=3)
D = ParabolicDGA(W, P, Plist, Delta)

# Calcular cohomología
H1 = D.cohomology_basis(1, ring=GF(2))

# Producto copa (ya corregido)
u, v = H1[0], H1[1]
cup = D.cup_on_cochains(u, v, 1, 1, ring=GF(2))
```

## 📂 Estructura del Repositorio

```
ParabolicArrangements/
├── parabolic_dga_notebook.sage    ⭐ Principal - usar este
├── parabolic_dga.sage              Referencia optimizada
├── analisis_cancelaciones.sage     Análisis de cancelaciones
├── intersection_dual_optimized.sage Intersección dual
├── parabolic_simplicial_dga.sage   Versión dual (lenta)
│
├── README.md                       Documentación principal
├── CHANGELOG.md                    Historial
├── INDEX_ARCHIVOS.md              Este archivo
│
└── .docs_session/                 Documentación de sesión (oculta para git)
    ├── USO_RAPIDO.md
    ├── SOLUCION_PRODUCTO_COPA.md
    ├── TEORIA_DUALIZACION.md
    ├── ROADMAP_TEORIA_INTERSECCION.md
    ├── RESUMEN_SESION_FINAL.md
    ├── QUICKSTART.md
    ├── ProductoCopa_Corregido.sage (obsoleto)
    └── ... (tests y demos)
```

## 🎯 ¿Qué Archivo Usar?

| Necesito... | Archivo |
|-------------|---------|
| Implementación completa lista | `parabolic_dga_notebook.sage` ⭐ |
| Referencia optimizada | `parabolic_dga.sage` |
| Entender cancelaciones | `analisis_cancelaciones.sage` |
| Intersección dual | `intersection_dual_optimized.sage` |
| Guía de uso | `.docs_session/USO_RAPIDO.md` |
| Teoría completa | `.docs_session/TEORIA_DUALIZACION.md` |
| Entender la corrección | `.docs_session/SOLUCION_PRODUCTO_COPA.md` |

## ✅ Garantías

**`parabolic_dga_notebook.sage`** implementa:
- ✅ Suma sobre TODAS las particiones (no solo una)
- ✅ Signos de shuffle correctos ε(J,K)
- ✅ Transformación por w₀_J
- ✅ Leibniz satisfecho
- ✅ Cohomología correcta
- ✅ Optimizado con caché

---

*Última actualización: Febrero 2026*
