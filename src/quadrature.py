"""Cuadratura angular de Gauss-Legendre para S_N.

El archivo ``cuadraturas.py`` de la cátedra lista solo la mitad positiva de los
nodos (``mu_m``) con sus pesos (``omega_m``), para N = 2, 4, 6 y 8. Aquí se
espeja esa mitad para construir el conjunto completo, de modo que

* ``mu`` queda ordenado de izquierda a derecha (ascendente, con signo), igual
  que las tablas de ``resolucion_ejercicio_s2.md`` / ``_s4.md``;
* ``w`` es simétrico (``w(+mu) = w(-mu)``) y suma 2.

Si el N pedido no está en el archivo, los nodos se generan por Newton sobre el
polinomio de Legendre P_N.
"""

from __future__ import annotations

import importlib.util
import math
import pathlib

#: Ruta por defecto al archivo de cuadraturas de la cátedra.
DEFAULT_TABLE = pathlib.Path(__file__).resolve().parent.parent / "cuadraturas.py"


class Quadrature:
    """Conjunto completo de direcciones y pesos de una cuadratura S_N."""

    __slots__ = ("n", "k", "mu", "w", "mu_half", "w_half", "source", "warnings")

    def __init__(
        self,
        mu_half: list[float],
        w_half: list[float],
        source: str = "generated",
        warnings: list[str] | None = None,
    ) -> None:
        if len(mu_half) != len(w_half):
            raise ValueError("mu_half y w_half deben tener la misma longitud")
        if not mu_half:
            raise ValueError("la cuadratura no puede estar vacía")
        if any(m <= 0.0 for m in mu_half):
            raise ValueError("mu_half debe contener solo nodos positivos")
        if any(w <= 0.0 for w in w_half):
            raise ValueError("los pesos deben ser positivos")

        half_mu = list(mu_half)
        half_w = list(w_half)
        # media cuadratura: descendente (como en cuadraturas.py)
        order = sorted(range(len(half_mu)), key=lambda i: -half_mu[i])
        self.mu_half = [half_mu[i] for i in order]
        self.w_half = [half_w[i] for i in order]
        self.k = len(self.mu_half)
        self.n = 2 * self.k
        # conjunto completo: ascendente y simétrico
        self.mu = [-m for m in self.mu_half] + self.mu_half[::-1]
        self.w = list(self.w_half) + self.w_half[::-1]
        self.source = source
        self.warnings = list(warnings or [])

    # ------------------------------------------------------------------
    @property
    def total_weight(self) -> float:
        """Suma de todos los pesos (debe ser 2)."""
        return sum(self.w)

    @property
    def half_weight(self) -> float:
        """Suma de la mitad positiva (debe ser 1)."""
        return sum(self.w_half)

    def positive_index(self) -> list[int]:
        """Índices m con ``mu_m > 0`` (direcciones entrantes por la izquierda)."""
        return [m for m, mu in enumerate(self.mu) if mu > 0.0]

    def negative_index(self) -> list[int]:
        """Índices m con ``mu_m < 0`` (direcciones entrantes por la derecha)."""
        return [m for m, mu in enumerate(self.mu) if mu < 0.0]

    def rows(self) -> list[tuple[int, float, float]]:
        """``[(m, mu_m, w_m), ...]`` con m empezando en 1."""
        return [(i + 1, mu, w) for i, (mu, w) in enumerate(zip(self.mu, self.w))]

    def __repr__(self) -> str:  # pragma: no cover - ayuda de depuración
        return f"Quadrature(n={self.n}, source={self.source!r})"


# ---------------------------------------------------------------------- api
def gauss_legendre(n: int) -> tuple[list[float], list[float]]:
    """Nodos (mitad positiva, descendente) y pesos de Gauss-Legendre en [-1, 1].

    Usa Newton sobre ``P_n(x)`` con la estimación inicial clásica. El resultado
    se simetriza por construcción (solo se calculan ``n/2`` nodos), lo que
    garantiza ``w(+mu) = w(-mu)`` exacto y que el polinomio característico
    salga con coeficientes limpios (p. ej. ``63 a^4 - 733 a^2 + 180`` para N=4).
    """
    if n < 2 or n % 2 != 0:
        raise ValueError(f"S_N requiere N par y >= 2 (recibido {n})")

    half = n // 2
    mus: list[float] = []
    ws: list[float] = []
    for i in range(1, half + 1):
        x = math.cos(math.pi * (i - 0.25) / (n + 0.5))
        derivative = 1.0
        for _ in range(100):
            p_prev, p_curr = 1.0, x
            for degree in range(2, n + 1):
                p_prev, p_curr = p_curr, (
                    (2 * degree - 1) * x * p_curr - (degree - 1) * p_prev
                ) / degree
            derivative = n * (x * p_curr - p_prev) / (x * x - 1.0)
            step = -p_curr / derivative
            x += step
            if abs(step) < 1e-16:
                break
        mus.append(abs(x))
        ws.append(2.0 / ((1.0 - x * x) * derivative * derivative))

    order = sorted(range(half), key=lambda i: -mus[i])
    return [mus[i] for i in order], [ws[i] for i in order]


def load_table(n: int, path: pathlib.Path | str = DEFAULT_TABLE) -> tuple[list[float], list[float]]:
    """Lee la mitad positiva de ``cuadraturas.py`` para el orden ``n``."""
    path = pathlib.Path(path)
    if not path.exists():
        raise FileNotFoundError(f"no se encontró el archivo de cuadraturas: {path}")
    spec = importlib.util.spec_from_file_location("_sn_cuadraturas", path)
    if spec is None or spec.loader is None:  # pragma: no cover - defensivo
        raise ImportError(f"no se pudo cargar {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    data = getattr(module, "DATA", None)
    if not isinstance(data, dict):
        raise ValueError(f"{path} no expone un diccionario DATA")

    mus: list[float] = []
    ws: list[float] = []
    for nn, _m, mu, w in zip(data["N"], data["m"], data["mu_m"], data["omega_m"]):
        if int(nn) == n:
            mus.append(float(mu))
            ws.append(float(w))
    if not mus:
        raise ValueError(f"cuadraturas.py no contiene N={n} (disponibles: {sorted(set(data['N']))})")
    return mus, ws


def _discrepancies(mus: list[float], ws: list[float], ref_mu: list[float], ref_w: list[float]) -> list[str]:
    """Diferencias entre una cuadratura dada y la generada (referencia)."""
    notes: list[str] = []
    if len(mus) != len(ref_mu):
        notes.append(f"la tabla tiene {len(mus)} nodos y la generada {len(ref_mu)}")
        return notes
    for i, (mu, w) in enumerate(zip(mus, ws)):
        if abs(mu - ref_mu[i]) > 1e-13:
            notes.append(f"nodo {i + 1}: tabla {mu:.16f} vs generado {ref_mu[i]:.16f}")
        if abs(w - ref_w[i]) > 1e-13:
            notes.append(
                f"peso {i + 1}: tabla {w:.16f} vs generado {ref_w[i]:.16f} "
                f"(diferencia {abs(w - ref_w[i]):.2e})"
            )
    return notes


def build_quadrature(
    n: int,
    source: str = "auto",
    path: pathlib.Path | str = DEFAULT_TABLE,
) -> Quadrature:
    """Construye la cuadratura ``S_n``.

    Args:
        n: orden (par, >= 2).
        source: ``"auto"`` (usa la tabla si coincide con la generada; si no,
            genera y avisa), ``"table"`` (obliga a leer el archivo) o
            ``"generated"`` (Newton sobre P_N).

    Notes:
        La tabla de la cátedra tiene un peso de N=6 con un dígito cambiado
        (``0.4679139345866910`` en lugar de ``0.4679139345726910``). Con
        ``source="auto"`` se detecta y se usa la cuadratura generada, dejando
        el aviso en ``Quadrature.warnings``.
    """
    if n < 2 or n % 2 != 0:
        raise ValueError(f"S_N requiere N par y >= 2 (recibido {n})")
    if source not in {"auto", "table", "generated"}:
        raise ValueError(f"source debe ser auto|table|generated (recibido {source!r})")

    ref_mu, ref_w = gauss_legendre(n)
    ref = Quadrature(ref_mu, ref_w, source="generated")
    if abs(ref.half_weight - 1.0) > 1e-12:
        raise ValueError(
            f"los pesos generados para N={n} suman {ref.half_weight:.12f}, se esperaba 1"
        )

    if source == "generated":
        return ref

    try:
        mus, ws = load_table(n, path)
    except (FileNotFoundError, ValueError):
        if source == "table":
            raise
        return ref

    table = Quadrature(mus, ws, source="table")
    notes = _discrepancies(table.mu_half, table.w_half, ref.mu_half, ref.w_half)
    if abs(table.half_weight - 1.0) > 1e-8:
        notes.append(
            f"los pesos de la tabla suman {table.half_weight:.12f}, se esperaba 1"
        )

    if not notes:
        return table
    if source == "table":
        table.warnings = notes + [
            "se usa la tabla tal como está; 'auto' preferiría la cuadratura generada"
        ]
        return table

    ref.warnings = [f"tabla de cuadraturas.py inconsistente para N={n}, se usa la generada"] + notes
    return ref
