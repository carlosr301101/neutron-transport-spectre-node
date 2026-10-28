"""Camino B: relación de dispersión, polinomio característico y autofunciones.

Se resuelve el problema espectral de los apuntes (ecs. 9 y 10)

    1 = (c*nu/2) * sum_n w_n / (mu_n + nu)          (dispersión)

o, en la convención de los documentos (``psi = v e^{-x/nu}``, ``a = sigma_t*nu``),

    sum_m w_m / (1 - mu_m/a) = 2/c

y las autofunciones angulares, normalizadas a flujo escalar modal unitario,

    v_m(nu) = (c/2) / (1 - mu_m/a).

Para ``c > 1`` la fundamental deja de ser real: se resuelve el par complejo
``a = +-i*omega`` con la ecuación real

    sum_m w_m / (1 + mu_m^2/omega^2) = 2/c.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from fractions import Fraction

from .medium import Medium
from .quadrature import Quadrature

#: Tolerancias internas del buscador de raíces.
RTOL = 1e-15
MAX_ITER = 400


# --------------------------------------------------------------------- utils
def dispersion(a: complex, quad: Quadrature, c: float) -> complex:
    """``sum_m w_m a/(a - mu_m) - 2/c`` (real si ``a`` es real)."""
    total = 0j
    for mu, w in zip(quad.mu, quad.w):
        total += w * a / (a - mu)
    return total - 2.0 / c


def _bisect(f, lo: float, hi: float) -> float:
    """Bisección con cambio de signo garantizado en ``[lo, hi]``."""
    f_lo = f(lo)
    f_hi = f(hi)
    if f_lo == 0.0:
        return lo
    if f_hi == 0.0:
        return hi
    if (f_lo > 0.0) == (f_hi > 0.0):
        raise ValueError("la bisección requiere un cambio de signo en el intervalo")
    for _ in range(MAX_ITER):
        mid = 0.5 * (lo + hi)
        f_mid = f(mid)
        if f_mid == 0.0 or (hi - lo) <= RTOL * max(1.0, abs(mid)):
            return mid
        if (f_mid > 0.0) == (f_lo > 0.0):
            lo, f_lo = mid, f_mid
        else:
            hi, f_hi = mid, f_mid
    return 0.5 * (lo + hi)


# ---------------------------------------------------------------- autovalores
def real_roots_positive(c: float, quad: Quadrature) -> list[float]:
    """Raíces reales ``a > 0`` de la dispersión, ordenadas de mayor a menor.

    Para ``0 < c < 1`` hay ``N/2`` raíces: la fundamental en ``(mu_1, inf)`` y
    una transitoria en cada intervalo ``(mu_k, mu_{k-1})``. Para ``c >= 1`` la
    fundamental no es real y solo quedan las ``N/2 - 1`` transitorias.
    """
    if c <= 0.0:
        raise ValueError(f"c debe ser > 0 para tener espectro (recibido {c})")

    f = lambda a: float(dispersion(a, quad, c).real)  # noqa: E731
    roots: list[float] = []

    if c < 1.0:
        lo = quad.mu_half[0] * (1.0 + 1e-12)
        hi = max(quad.mu_half[0] * 2.0, 1.0)
        guard = 0
        while f(hi) > 0.0:
            hi *= 2.0
            guard += 1
            if guard > 400:  # pragma: no cover - c extremadamente cercano a 1
                raise RuntimeError("no se pudo encerrar la raíz fundamental")
        roots.append(_bisect(f, lo, hi))

    for k in range(1, quad.k):
        lo = quad.mu_half[k] * (1.0 + 1e-12)
        hi = quad.mu_half[k - 1] * (1.0 - 1e-12)
        roots.append(_bisect(f, lo, hi))

    return roots


def complex_fundamental(c: float, quad: Quadrature) -> float:
    """``omega > 0`` de la fundamental compleja ``a = +-i*omega`` (solo ``c > 1``)."""
    if c <= 1.0:
        raise ValueError("la fundamental compleja requiere c > 1")

    def f(omega: float) -> float:
        total = 0.0
        for mu, w in zip(quad.mu, quad.w):
            total += w / (1.0 + (mu / omega) ** 2)
        return total - 2.0 / c

    lo = max(quad.mu_half[0], 1e-3) * 1e-6
    hi = max(quad.mu_half[0], 1e-3) * 2.0
    while f(hi) < 0.0:
        hi *= 2.0
    return _bisect(f, lo, hi)


def eigenvector(a: complex, quad: Quadrature, c: float) -> list[complex]:
    """Autofunción angular ``v_m = (c/2)/(1 - mu_m/a)`` (compleja si ``a`` lo es)."""
    return [(c / 2.0) / (1.0 - mu / a) for mu in quad.mu]


# ------------------------------------------------------- polinomio exacto
def _convolve(p: list[float], q: list[float]) -> list[float]:
    out = [0.0] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        for j, b in enumerate(q):
            out[i + j] += a * b
    return out


def _synthetic_div(p: list[float], root: float) -> list[float]:
    """Divide ``p`` (descendente) por ``(y - root)`` y devuelve el cociente."""
    out = [p[0]]
    for coef in p[1:-1]:
        out.append(coef + out[-1] * root)
    return out


def characteristic_polynomial(quad: Quadrature, c: float) -> list[float]:
    """Coeficientes (descendentes) del polinomio característico en ``y = a^2``.

    Se obtiene de ``sum_k w_k y/(y - mu_k^2) = 1/c`` multiplicando por
    ``prod_k (y - mu_k^2)``. El grado es ``N/2``.
    """
    if c <= 0.0:
        raise ValueError(f"c debe ser > 0 (recibido {c})")

    q = [1.0]
    for mu in quad.mu_half:
        q = _convolve(q, [1.0, -mu * mu])

    acc = [0.0] * (quad.k + 1)
    for mu, w in zip(quad.mu_half, quad.w_half):
        quotient = _synthetic_div(q, mu * mu)
        # en orden descendente, multiplicar por y equivale a añadir un cero al final
        term = quotient + [0.0]
        for i, coef in enumerate(term):
            acc[i] += w * coef
    for i, coef in enumerate(q):
        acc[i] -= coef / c
    return acc


def _snap(value: float, max_den: int = 100_000, tol: float = 1e-12) -> Fraction | None:
    """Aproxima ``value`` por una fracción simple, si está cerca."""
    if abs(value) < 1e-14:
        return Fraction(0)
    frac = Fraction(value).limit_denominator(max_den)
    if abs(float(frac) - value) <= tol * max(1.0, abs(value)):
        return frac
    return None


def exact_coefficients(coeffs: list[float]) -> list[Fraction] | None:
    """Normaliza los coeficientes a enteros primos entre sí, si son racionales."""
    snapped = [_snap(x) for x in coeffs]
    if any(f is None for f in snapped):
        return None
    lcm = 1
    for frac in snapped:
        lcm = lcm * frac.denominator // math.gcd(lcm, frac.denominator)
    ints = [int(frac * lcm) for frac in snapped]
    common = 0
    for value in ints:
        common = math.gcd(common, abs(value))
    if common > 1:
        ints = [value // common for value in ints]
    if ints[0] < 0:
        ints = [-value for value in ints]
    return [Fraction(value) for value in ints]


# -------------------------------------------------------------------- modelo
@dataclass
class Mode:
    """Un modo espectral: autovalor ``nu`` y autofunción angular ``v``."""

    index: int
    nu: complex
    a: complex
    v: list[complex]
    kind: str

    @property
    def is_complex(self) -> bool:
        return abs(self.nu.imag) > 1e-10

    @property
    def decay(self) -> float | None:
        """``1/nu`` en cm^-1 para modos reales; ``None`` si es complejo."""
        if self.is_complex or self.nu.real == 0.0:
            return None
        return 1.0 / self.nu.real


@dataclass
class Spectrum:
    """Espectro completo del problema (2*N/2 modos)."""

    c: float
    sigma_t: float
    quad: Quadrature
    modes: list[Mode]
    poly: list[float]
    poly_exact: list[Fraction] | None = None
    notes: list[str] = field(default_factory=list)

    @property
    def nu(self) -> list[complex]:
        return [mode.nu for mode in self.modes]

    @property
    def fundamental(self) -> Mode:
        return self.modes[0]

    def by_kind(self, kind: str) -> list[Mode]:
        return [mode for mode in self.modes if mode.kind == kind]


def solve_spectrum(medium: Medium, quad: Quadrature) -> Spectrum:
    """Autovalores y autofunciones por el camino de la dispersión."""
    c = medium.c
    if c <= 0.0:
        raise ValueError(f"c = {c:.6g} <= 0: el medio no tiene espectro físico")

    notes: list[str] = []
    pairs: list[tuple[complex, str]] = []

    if abs(c - 1.0) < 1e-12:
        notes.append("c = 1: no hay modo fundamental finito (nu_1 -> infinito)")
        for a in real_roots_positive(c, quad):
            pairs.append((complex(a, 0.0), "transitorio"))
    elif c < 1.0:
        roots = real_roots_positive(c, quad)
        pairs.append((complex(roots[0], 0.0), "fundamental"))
        for a in roots[1:]:
            pairs.append((complex(a, 0.0), "transitorio"))
    else:
        omega = complex_fundamental(c, quad)
        pairs.append((complex(0.0, omega), "fundamental compleja"))
        for a in real_roots_positive(c, quad):
            pairs.append((complex(a, 0.0), "transitorio"))

    modes: list[Mode] = []
    index = 0
    for a, kind in pairs:
        for signed in (a, -a):
            index += 1
            modes.append(
                Mode(
                    index=index,
                    nu=signed / medium.sigma_t,
                    a=signed,
                    v=eigenvector(signed, quad, c),
                    kind=kind,
                )
            )

    coeffs = characteristic_polynomial(quad, c)
    return Spectrum(
        c=c,
        sigma_t=medium.sigma_t,
        quad=quad,
        modes=modes,
        poly=coeffs,
        poly_exact=exact_coefficients(coeffs),
        notes=notes,
    )


# ----------------------------------------------------------------- chequeos
def checks(spec: Spectrum) -> dict[str, float]:
    """Invariantes numéricos del espectro (todos deberían ser ~0)."""
    c = spec.c
    quad = spec.quad
    norm = 0.0
    disp = 0.0
    for mode in spec.modes:
        total = sum(w * v for w, v in zip(quad.w, mode.v))
        norm = max(norm, abs(total - 1.0))
        disp = max(disp, abs(dispersion(mode.a, quad, c)))

    mirror = 0.0
    n = quad.n
    for i in range(0, len(spec.modes), 2):
        plus = spec.modes[i].v
        minus = spec.modes[i + 1].v
        for m in range(n):
            mirror = max(mirror, abs(minus[m] - plus[n - 1 - m]))

    return {"normalization": norm, "dispersion": disp, "mirror": mirror}


def canonical_order(nus: list[complex]) -> list[complex]:
    """Ordena el espectro: mayor ``|nu|`` primero y, en cada par, ``+nu`` antes que ``-nu``."""
    def key(z: complex):
        magnitude = -round(abs(z), 9)
        positive_first = 0 if (z.real > 0.0 or (z.real == 0.0 and z.imag >= 0.0)) else 1
        return (magnitude, positive_first)

    return sorted(nus, key=key)


def continuum_nu(c: float, sigma_t: float = 1.0) -> float | None:
    """Autovalor fundamental del transporte continuo (referencia teórica).

    Resuelve ``(c*a/2) * ln((a+1)/(a-1)) = 1`` con ``a = sigma_t*nu``. Devuelve
    ``None`` si ``c >= 1`` (no hay raíz real en ese caso).
    """
    if c <= 0.0 or c >= 1.0:
        return None

    def f(a: float) -> float:
        return 0.5 * c * a * math.log((a + 1.0) / (a - 1.0)) - 1.0

    lo, hi = 1.0 + 1e-12, 2.0
    while f(hi) > 0.0:
        hi *= 2.0
    return _bisect(f, lo, hi) / sigma_t
