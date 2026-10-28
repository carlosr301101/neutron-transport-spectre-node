"""Medio homogéneo 1D y cálculo de la constante ``c``.

Las dos ramas están tomadas de ``calculo-c-metodo-espectral.md``:

* no multiplicativo:  ``c = sigma_s0 / sigma_t``
* multiplicativo:     ``c = sigma_s0/sigma_t + 4*pi*nu_bar*sigma_f/(k_eff*sigma_t)``

``nu_bar`` es el número medio de neutrones por fisión. Se llama así (y no
``nu``) para no confundirlo con ``nu``, el autovalor espectral: el mismo
símbolo tiene dos significados en la literatura y aquí se separan a propósito.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

#: Factor del término de fisión usado en la cátedra (``4*math.pi``).
FISSION_FACTOR = 4.0 * math.pi


@dataclass
class Medium:
    """Medio homogéneo de una sola energía.

    Attributes:
        sigma_t: sección eficaz macroscópica total (cm^-1), debe ser > 0.
        sigma_s0: sección eficaz de dispersión isotrópica (cm^-1).
        sigma_f: sección eficaz de fisión (cm^-1).
        nu_bar: neutrones producidos por fisión.
        k_eff: factor de multiplicación efectivo; obligatorio si hay fisión.
        fission_factor: peso del término de fisión (por defecto ``4*pi``).
    """

    sigma_t: float
    sigma_s0: float = 0.0
    sigma_f: float = 0.0
    nu_bar: float = 0.0
    k_eff: float | None = None
    fission_factor: float = FISSION_FACTOR

    def __post_init__(self) -> None:
        if self.sigma_t <= 0.0:
            raise ValueError(f"sigma_t debe ser > 0 (recibido {self.sigma_t})")
        for name in ("sigma_s0", "sigma_f", "nu_bar"):
            value = getattr(self, name)
            if value < 0.0:
                raise ValueError(f"{name} no puede ser negativo (recibido {value})")
        if self.sigma_f > 0.0:
            if self.k_eff is None:
                raise ValueError(
                    "con sigma_f > 0 hace falta k_eff (ver calculo-c-metodo-espectral.md)"
                )
            if self.k_eff <= 0.0:
                raise ValueError(f"k_eff debe ser > 0 (recibido {self.k_eff})")
        elif self.k_eff is not None and self.k_eff <= 0.0:
            raise ValueError(f"k_eff debe ser > 0 (recibido {self.k_eff})")

    # ------------------------------------------------------------------ c
    @property
    def has_fission(self) -> bool:
        """True si el medio es multiplicativo."""
        return self.sigma_f > 0.0

    @property
    def c(self) -> float:
        """Constante de segundo vuelo, según la rama que corresponda."""
        base = self.sigma_s0 / self.sigma_t
        if not self.has_fission:
            return base
        return base + (
            self.fission_factor * self.nu_bar * self.sigma_f / (self.k_eff * self.sigma_t)
        )

    # --------------------------------------------------------------- utils
    @classmethod
    def from_c(cls, c: float, sigma_t: float = 1.0) -> "Medium":
        """Medio no multiplicativo con un ``c`` impuesto (``sigma_s0 = c*sigma_t``)."""
        return cls(sigma_t=sigma_t, sigma_s0=c * sigma_t)

    def warnings(self) -> list[str]:
        """Advertencias físicas que el reporte muestra tal cual."""
        out: list[str] = []
        if self.has_fission and self.nu_bar == 0.0:
            out.append("sigma_f > 0 pero nu_bar = 0: la fisión no aporta neutrones")
        c = self.c
        if c <= 0.0:
            out.append("c <= 0: no hay espectro físico (el medio solo absorbe)")
        elif abs(c - 1.0) < 1e-12:
            out.append("c = 1: sin atenuación fundamental (nu_1 -> infinito)")
        elif c > 1.0:
            out.append("c > 1: medio supercrítico, la fundamental pasa a ser compleja (+-i*omega)")
        return out

    def describe(self) -> str:
        """Explica de dónde sale ``c``, con la sustitución numérica."""
        lines: list[str] = []
        if self.has_fission:
            lines.append("Medio multiplicativo")
            lines.append(
                f"  sigma_t = {self.sigma_t:.4f}   sigma_s0 = {self.sigma_s0:.4f}"
                f"   sigma_f = {self.sigma_f:.4f}"
            )
            lines.append(
                f"  nu_bar  = {self.nu_bar:.4f}   k_eff   = {self.k_eff:.4f}"
                f"   factor  = {self.fission_factor:.4f}  (4*pi = {FISSION_FACTOR:.4f})"
            )
            lines.append("  c = sigma_s0/sigma_t + factor*nu_bar*sigma_f/(k_eff*sigma_t)")
            lines.append(
                f"    = {self.sigma_s0:.4f}/{self.sigma_t:.4f} + "
                f"{self.fission_factor:.4f}*{self.nu_bar:.4f}*{self.sigma_f:.4f}"
                f"/({self.k_eff:.4f}*{self.sigma_t:.4f})"
            )
            lines.append(f"    = {self.c:.6f}")
        else:
            lines.append("Medio homogéneo no multiplicativo")
            lines.append(f"  sigma_t = {self.sigma_t:.4f}   sigma_s0 = {self.sigma_s0:.4f}")
            lines.append(
                f"  c = sigma_s0/sigma_t = {self.sigma_s0:.4f}/{self.sigma_t:.4f}"
                f" = {self.c:.6f}"
            )
        for warning in self.warnings():
            lines.append(f"  [!] {warning}")
        return "\n".join(lines)

    def to_dict(self) -> dict:
        """Representación serializable (para ``--json``)."""
        return {
            "sigma_t": self.sigma_t,
            "sigma_s0": self.sigma_s0,
            "sigma_f": self.sigma_f,
            "nu_bar": self.nu_bar,
            "k_eff": self.k_eff,
            "fission_factor": self.fission_factor,
            "has_fission": self.has_fission,
            "c": self.c,
            "warnings": self.warnings(),
        }
