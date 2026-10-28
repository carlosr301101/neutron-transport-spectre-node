"""Camino A: problema de autovalores matricial ``A_N a = (1/nu) a`` (ec. 8).

La matriz es la de los apuntes (pág. 3), que usan el ansatz
``psi = a e^{+sigma_t x/nu}``:

    (A_N)_{mn} = c*w_n / (2*mu_m)          si m != n
    (A_N)_{mm} = (c*w_m - 2) / (2*mu_m)

Los autovalores de ``A_N`` son ``lambda = 1/nu``. Es un camino independiente del
de la relación de dispersión y sirve para contrastarlo.

La matriz no depende de ``sigma_t``: al dividir la ecuación de transporte por
``sigma_t``, el cociente ``mu_m/sigma_t`` multiplica a ``d/dx`` y el factor
``sigma_t`` del exponencial lo cancela.
"""

from __future__ import annotations

import cmath
import math

from .medium import Medium
from .quadrature import Quadrature


# ------------------------------------------------------------------- matriz
def build_matrix(c: float, quad: Quadrature) -> list[list[complex]]:
    """Construye ``A_N`` para un ``c`` y una cuadratura dados."""
    n = quad.n
    matrix: list[list[complex]] = []
    for m in range(n):
        mu = quad.mu[m]
        row: list[complex] = []
        for j in range(n):
            if j == m:
                row.append(complex((c * quad.w[m] - 2.0) / (2.0 * mu)))
            else:
                row.append(complex(c * quad.w[j] / (2.0 * mu)))
        matrix.append(row)
    return matrix


# ------------------------------------------------------- álgebra numérica
def _hessenberg(matrix: list[list[complex]]) -> list[list[complex]]:
    """Reduce a forma de Hessenberg superior con Householder complejas."""
    n = len(matrix)
    h = [row[:] for row in matrix]
    for k in range(n - 2):
        x = [h[i][k] for i in range(k + 1, n)]
        norm_x = math.sqrt(sum(abs(v) ** 2 for v in x))
        if norm_x == 0.0:
            continue
        phase = x[0] / abs(x[0]) if abs(x[0]) > 0.0 else 1.0 + 0.0j
        v = list(x)
        v[0] += phase * norm_x
        norm_v = math.sqrt(sum(abs(t) ** 2 for t in v))
        if norm_v == 0.0:
            continue
        v = [t / norm_v for t in v]
        length = len(v)

        # h <- (I - 2 v v^H) h   (izquierda)
        for j in range(k, n):
            s = sum(v[i].conjugate() * h[k + 1 + i][j] for i in range(length))
            for i in range(length):
                h[k + 1 + i][j] -= 2.0 * v[i] * s
        # h <- h (I - 2 v v^H)   (derecha)
        for i in range(n):
            s = sum(h[i][k + 1 + j] * v[j] for j in range(length))
            for j in range(length):
                h[i][k + 1 + j] -= 2.0 * s * v[j].conjugate()
    return h


def _qr_step(h: list[list[complex]], shift: complex) -> list[list[complex]]:
    """Un paso de iteración QR con desplazamiento sobre una matriz de Hessenberg."""
    n = len(h)
    a = [row[:] for row in h]
    for i in range(n):
        a[i][i] -= shift

    rotations: list[tuple[complex, complex]] = []
    for k in range(n - 1):
        top, bottom = a[k][k], a[k + 1][k]
        abs_top = abs(top)
        if abs(bottom) == 0.0:
            rotations.append((1.0 + 0.0j, 0.0j))
            continue
        radius = math.hypot(abs_top, abs(bottom))
        if abs_top == 0.0:
            cosine = 0.0 + 0.0j
            sine = bottom.conjugate() / abs(bottom)
        else:
            cosine = complex(abs_top / radius)
            sine = (top / abs_top) * bottom.conjugate() / radius
        rotations.append((cosine, sine))
        for j in range(k, n):
            r1, r2 = a[k][j], a[k + 1][j]
            a[k][j] = cosine * r1 + sine * r2
            a[k + 1][j] = -sine.conjugate() * r1 + cosine * r2

    # R * Q  (con Q = G^H, aplicado por la derecha)
    for k in range(n - 1):
        cosine, sine = rotations[k]
        for i in range(n):
            r1, r2 = a[i][k], a[i][k + 1]
            a[i][k] = cosine * r1 + sine.conjugate() * r2
            a[i][k + 1] = -sine * r1 + cosine * r2

    for i in range(n):
        a[i][i] += shift
    return a


def _wilkinson_shift(block: list[list[complex]]) -> complex:
    """Autovalor de la submatriz 2x2 final más cercano a la esquina inferior."""
    m = len(block)
    if m == 1:
        return block[0][0]
    a, b = block[m - 2][m - 2], block[m - 2][m - 1]
    c, d = block[m - 1][m - 2], block[m - 1][m - 1]
    trace = 0.5 * (a + d)
    disc = cmath.sqrt(trace * trace - (a * d - b * c))
    first, second = trace + disc, trace - disc
    return first if abs(first - d) <= abs(second - d) else second


def eigenvalues(matrix: list[list[complex]], tol: float = 1e-13, max_iter: int = 3000) -> list[complex]:
    """Autovalores de una matriz compleja por Hessenberg + QR con desplazamiento."""
    n = len(matrix)
    if n == 0:
        return []
    h = _hessenberg(matrix)
    found: list[complex] = []
    size = n
    iterations = 0
    while size > 1:
        subdiagonal = abs(h[size - 1][size - 2])
        scale = abs(h[size - 2][size - 2]) + abs(h[size - 1][size - 1]) + 1e-300
        if subdiagonal <= tol * scale:
            found.append(h[size - 1][size - 1])
            size -= 1
            iterations = 0
            continue
        block = [row[:size] for row in h[:size]]
        updated = _qr_step(block, _wilkinson_shift(block))
        for i in range(size):
            h[i][:size] = updated[i]
        iterations += 1
        if iterations > max_iter:
            raise RuntimeError(
                f"la iteración QR no convergió para el bloque de tamaño {size} "
                f"(residuo {subdiagonal:.3e})"
            )
    found.append(h[0][0])
    return found


# ---------------------------------------------------------------------- api
def _clean(z: complex, tol: float = 1e-10) -> complex:
    """Elimina la parte imaginaria de ruido numérico (``1e-16``) de un autovalor."""
    scale = max(1.0, abs(z))
    if abs(z.imag) <= tol * scale:
        return complex(z.real, 0.0)
    if abs(z.real) <= tol * scale:
        return complex(0.0, z.imag)
    return z


def path_a(medium: Medium, quad: Quadrature) -> list[complex]:
    """Autovalores ``nu`` del problema matricial, en unidades de los documentos.

    La matriz da ``nu_A`` con el ansatz ``e^{+sigma_t x/nu}``; los documentos
    usan ``e^{-x/nu}``, de modo que ``nu = nu_A/sigma_t`` (el signo lo absorbe la
    simetría ``+nu / -nu`` del espectro).
    """
    from .spectrum import canonical_order

    matrix = build_matrix(medium.c, quad)
    lambdas = eigenvalues(matrix)
    nus = [_clean(1.0 / (lam * medium.sigma_t)) for lam in lambdas]
    return canonical_order(nus)


def matrix_residual(medium: Medium, quad: Quadrature, spectrum) -> float:
    """Residuo relativo máximo de ``A_N a = a/(sigma_t*nu)``, con ``a = v(-nu)``.

    Es un contraste entre caminos: la autofunción del camino B (evaluada en
    ``-nu``) debe ser autovector de la matriz del camino A. Se normaliza por
    ``max|a|`` para que no dependa de la escala de la matriz.
    """
    matrix = build_matrix(medium.c, quad)
    n = quad.n
    modes = spectrum.modes
    worst = 0.0
    for j in range(0, len(modes), 2):
        for mode, partner in ((modes[j], modes[j + 1]), (modes[j + 1], modes[j])):
            vec = partner.v
            eigenvalue = 1.0 / (medium.sigma_t * mode.nu)
            residual = 0.0
            for m in range(n):
                accumulated = sum(matrix[m][k] * vec[k] for k in range(n))
                residual = max(residual, abs(accumulated - eigenvalue * vec[m]))
            scale = max(abs(z) for z in vec) or 1.0
            worst = max(worst, residual / scale)
    return worst
