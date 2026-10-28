"""Reportes del análisis espectral (texto plano o Markdown)."""

from __future__ import annotations

from fractions import Fraction

from .medium import Medium
from .quadrature import Quadrature
from .spectrum import Spectrum, checks, continuum_nu


# ---------------------------------------------------------------- formateo
def fmt_num(value: complex | float, digits: int = 7) -> str:
    """Formatea un real o complejo sin ruido numérico."""
    z = complex(value)
    scale = max(1.0, abs(z))
    if abs(z.imag) <= 1e-12 * scale:
        return f"{z.real:.{digits}f}"
    if abs(z.real) <= 1e-12 * scale:
        sign = "+" if z.imag >= 0 else "-"
        return f"{sign}i*{abs(z.imag):.{digits}f}"
    sign = "+" if z.imag >= 0 else "-"
    return f"{z.real:.{digits}f} {sign} i*{abs(z.imag):.{digits}f}"


def _table(headers: list[str], rows: list[list[str]], markdown: bool) -> str:
    if markdown:
        lines = ["| " + " | ".join(headers) + " |",
                 "|" + "|".join([":---:"] * len(headers)) + "|"]
        lines += ["| " + " | ".join(str(cell) for cell in row) + " |" for row in rows]
        return "\n".join(lines)
    widths = []
    for i, header in enumerate(headers):
        widths.append(max(len(header), *(len(str(row[i])) for row in rows)) if rows else len(header))
    sep = "  "
    lines = [sep.join(h.rjust(w) for h, w in zip(headers, widths))]
    lines.append(sep.join("-" * w for w in widths))
    lines += [sep.join(str(cell).rjust(w) for cell, w in zip(row, widths)) for row in rows]
    return "\n".join(lines)


# ----------------------------------------------------------------- bloques
def quadrature_table(quad: Quadrature, markdown: bool = True) -> str:
    rows = []
    for m, mu, w in quad.rows():
        sense = "izquierda" if mu < 0 else "derecha"
        rows.append([str(m), f"{mu:+.7f}", f"{w:.7f}", sense])
    return _table(["m", "mu_m", "w_m", "sentido"], rows, markdown)


def polynomial_equation(spec: Spectrum) -> str:
    """Ecuación característica en ``a = sigma_t*nu`` (``a^2 = y``)."""
    coeffs = spec.poly_exact if spec.poly_exact is not None else spec.poly
    degree = len(coeffs) - 1
    terms: list[str] = []
    for i, coef in enumerate(coeffs):
        exponent = 2 * (degree - i)
        value = float(coef) if isinstance(coef, Fraction) else coef
        if abs(value) < 1e-14:
            continue
        magnitude = abs(value)
        if isinstance(coef, Fraction) and coef.denominator == 1:
            body = f"{int(magnitude)}"
        else:
            body = f"{magnitude:.6g}"
        if exponent == 0:
            text = body
        elif exponent == 1:
            text = body if body == "1" else f"{body} a"
        else:
            text = f"{body} a^{exponent}" if body != "1" else f"a^{exponent}"
        if not terms:
            terms.append(("-" if value < 0 else "") + text)
        else:
            terms.append((" - " if value < 0 else " + ") + text)
    return "".join(terms) + " = 0"


def eigenvalue_table(spec: Spectrum, markdown: bool = True) -> str:
    rows = []
    for mode in spec.modes:
        decay = mode.decay
        decay_text = fmt_num(decay) if decay is not None else fmt_num(1.0 / mode.nu)
        rows.append([str(mode.index), fmt_num(mode.nu), decay_text, mode.kind])
    return _table(["j", "nu (cm)", "1/nu (1/cm)", "tipo"], rows, markdown)


def eigenvector_table(spec: Spectrum, markdown: bool = True) -> str:
    n = spec.quad.n
    rows = []
    for m in range(n):
        row = [str(m + 1), f"{spec.quad.mu[m]:+.7f}"]
        row += [fmt_num(mode.v[m]) for mode in spec.modes]
        rows.append(row)
    headers = ["m", "mu_m"] + [f"j={mode.index}" for mode in spec.modes]
    return _table(headers, rows, markdown)


def paths_table(spec: Spectrum, nus_a: list[complex], markdown: bool = True) -> str:
    rows = []
    for mode, other in zip(spec.modes, nus_a):
        rows.append([str(mode.index), fmt_num(other), fmt_num(mode.nu), f"{abs(other - mode.nu):.2e}"])
    return _table(["j", "nu (camino A)", "nu (camino B)", "|A-B|"], rows, markdown)


def checks_table(spec: Spectrum, residual_a: float, markdown: bool = True) -> str:
    values = checks(spec)
    rows = [
        ["sum_m w_m v_m - 1", f"{values['normalization']:.2e}", "~0"],
        ["residuo de dispersión", f"{values['dispersion']:.2e}", "~0"],
        ["simetría v(-nu)", f"{values['mirror']:.2e}", "~0"],
        ["residuo A_N (autovector)", f"{residual_a:.2e}", "camino A"],
    ]
    return _table(["chequeo", "valor", "referencia"], rows, markdown)


def convergence_table(medium: Medium, quads: list[Quadrature], markdown: bool = True) -> str:
    reference = continuum_nu(medium.c, medium.sigma_t)
    rows = []
    for quad in quads:
        spectrum = _solve_cached(medium, quad)
        nu1 = spectrum.modes[0].nu
        if abs(nu1.imag) > 1e-10:
            rows.append([str(quad.n), fmt_num(nu1), "-", "-"])
            continue
        if reference is None:
            error = "-"
        else:
            error = f"{abs(nu1.real - reference) / abs(reference):.2e}"
        rows.append([str(quad.n), fmt_num(nu1), fmt_num(1.0 / nu1.real), error])
    return _table(
        ["N", "nu_1 (cm)", "1/nu_1 (1/cm)", "error rel. vs continuo"], rows, markdown
    )


def _solve_cached(medium: Medium, quad: Quadrature) -> Spectrum:
    # import local para evitar ciclo en tiempo de importación
    from .spectrum import solve_spectrum

    return solve_spectrum(medium, quad)


# ----------------------------------------------------------------- reporte
def full_report(
    medium: Medium,
    quad: Quadrature,
    spec: Spectrum,
    nus_a: list[complex],
    residual_a: float,
    markdown: bool = True,
    include_medium: bool = True,
) -> str:
    """Informe completo de un orden S_N."""
    parts: list[str] = []
    title = f"S_{quad.n}  (cuadratura: {quad.source}, suma de pesos = {quad.total_weight:.12f})"
    parts.append(f"## {title}" if markdown else f"=== {title} ===")
    parts.append("")
    if include_medium:
        parts.append(medium.describe())
        parts.append("")
    if markdown:
        parts.append("### Cuadratura angular")
    else:
        parts.append("-- Cuadratura angular")
    parts.append(quadrature_table(quad, markdown))
    for warning in quad.warnings:
        parts.append(f"    [!] {warning}")
    parts.append("")
    if markdown:
        parts.append("### Ecuación característica")
    else:
        parts.append("-- Ecuación característica")
    parts.append(f"    {polynomial_equation(spec)}")
    if spec.poly_exact is None:
        parts.append("    (coeficientes en punto flotante: no hay forma entera simple)")
    parts.append("")
    for note in spec.notes:
        parts.append(f"    [!] {note}")
    if spec.notes:
        parts.append("")
    if markdown:
        parts.append("### Autovalores y autofunciones")
    else:
        parts.append("-- Autovalores")
    parts.append(eigenvalue_table(spec, markdown))
    parts.append("")
    if not markdown:
        parts.append("-- Autofunciones (columnas = modos; j=1,2 -> +nu1,-nu1; ...)")
    parts.append(eigenvector_table(spec, markdown))
    parts.append("")
    if markdown:
        parts.append("### Contraste de los dos caminos")
    else:
        parts.append("-- Contraste de los dos caminos")
    parts.append(paths_table(spec, nus_a, markdown))
    parts.append("")
    if markdown:
        parts.append("### Chequeos numéricos")
    else:
        parts.append("-- Chequeos numéricos")
    parts.append(checks_table(spec, residual_a, markdown))
    return "\n".join(parts)
