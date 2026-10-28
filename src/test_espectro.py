"""Regresión del análisis espectral contra los documentos S2 y S4.

Se ejecuta con ``python codigo/test_espectro.py`` desde ``metodos-espectral/``
(sin pytest ni dependencias externas).
"""

from __future__ import annotations

import io
import math
import pathlib
import sys
from contextlib import redirect_stdout

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))

from codigo import eigen  # noqa: E402
from codigo.__main__ import main  # noqa: E402
from codigo.medium import FISSION_FACTOR, Medium  # noqa: E402
from codigo.quadrature import build_quadrature, gauss_legendre  # noqa: E402
from codigo.spectrum import (  # noqa: E402
    checks,
    continuum_nu,
    solve_spectrum,
)

PASSED = 0


def check(name: str, condition: bool, detail: str = "") -> None:
    global PASSED
    if not condition:
        raise AssertionError(f"FALLA: {name} {detail}")
    PASSED += 1
    print(f"  ok  {name}" + (f"  ({detail})" if detail else ""))


def close(a: float, b: float, tol: float) -> bool:
    return abs(a - b) <= tol


def test_quadrature() -> None:
    print("cuadratura")
    for n in (2, 4, 8):
        table = build_quadrature(n, source="table")
        generated = build_quadrature(n, source="generated")
        check(
            f"N={n}: tabla == generado",
            max(abs(a - b) for a, b in zip(table.mu, generated.mu)) < 1e-13
            and max(abs(a - b) for a, b in zip(table.w, generated.w)) < 1e-13,
            f"max dmu={max(abs(a-b) for a,b in zip(table.mu, generated.mu)):.1e}",
        )
        check(f"N={n}: suma de pesos = 2", close(table.total_weight, 2.0, 1e-12))

    # N=6: la tabla de la cátedra trae un peso con un dígito cambiado
    table6 = build_quadrature(6, source="table")
    auto6 = build_quadrature(6, source="auto")
    check(
        "N=6 tabla: peso 3 desplazado 1.4e-11",
        abs(table6.w_half[2] - 0.4679139345726910) < 2e-11,
        f"w3={table6.w_half[2]:.16f}",
    )
    check("N=6 auto: detecta y usa la generada", auto6.source == "generated")
    check("N=6 auto: deja aviso", any("inconsistente" in w for w in auto6.warnings))

    check("N=4: nodos simétricos", close(build_quadrature(4).mu[0], -build_quadrature(4).mu[3], 0.0))
    mus, ws = gauss_legendre(10)
    check("N=10 generado: 5 nodos y suma 1", len(mus) == 5 and close(sum(ws), 1.0, 1e-12))


def test_s2() -> None:
    print("S2 (resolucion_ejercicio_s2.md)")
    medium = Medium.from_c(0.97)
    quad = build_quadrature(2)
    spec = solve_spectrum(medium, quad)

    check("nu = +-10/3", close(spec.modes[0].nu.real, 10 / 3, 1e-12))
    check("nu2 = -10/3", close(spec.modes[1].nu.real, -10 / 3, 1e-12))
    check("v1(+) = 0.4133975", close(spec.modes[0].v[0].real, 0.4133975, 1e-7))
    check("v2(+) = 0.5866025", close(spec.modes[0].v[1].real, 0.5866025, 1e-7))
    check("1/nu = 0.3", close(spec.modes[0].decay, 0.3, 1e-12))
    check("polinomio 9a^2 - 100", spec.poly_exact == [c for c in (9, -100)], str(spec.poly_exact))

    auto = eigen.path_a(medium, quad)
    check(
        "camino A coincide",
        max(abs(a - b) for a, b in zip(auto, spec.nu)) < 1e-12,
    )


def test_s4() -> None:
    print("S4 (resolucion_ejercicio_s4.md)")
    medium = Medium.from_c(0.97)
    quad = build_quadrature(4)
    spec = solve_spectrum(medium, quad)

    check(
        "polinomio 63a^4 - 733a^2 + 180",
        spec.poly_exact == [c for c in (63, -733, 180)],
        str(spec.poly_exact),
    )
    check("nu1 = 3.3740096", close(spec.modes[0].nu.real, 3.3740096, 1e-7))
    check("nu2 = 0.5009792", close(spec.modes[2].nu.real, 0.5009792, 1e-7))
    check("1/nu1 = 0.2963833", close(spec.modes[0].decay, 0.2963833, 1e-7))
    check("1/nu2 = 1.9960910", close(spec.modes[2].decay, 1.9960910, 1e-7))

    esperado = [
        [0.3863845, 0.6512046, 0.1783805, -0.6746358],
        [0.4406028, 0.5393472, 0.2889256, 1.5091784],
        [0.5393472, 0.4406028, 1.5091784, 0.2889256],
        [0.6512046, 0.3863845, -0.6746358, 0.1783805],
    ]
    error = max(
        abs(spec.modes[j].v[m].real - esperado[m][j])
        for m in range(4)
        for j in range(4)
    )
    check("tabla de autovectores del documento", error < 1e-7, f"max err={error:.1e}")

    auto = eigen.path_a(medium, quad)
    error_a = max(abs(a - b) for a, b in zip(auto, spec.nu))
    check("camino A coincide", error_a < 1e-12, f"max|A-B|={error_a:.1e}")
    residual = eigen.matrix_residual(medium, quad, spec)
    check("residuo A_N", residual < 1e-12, f"{residual:.1e}")


def test_convergencia() -> None:
    print("convergencia")
    medium = Medium.from_c(0.97)
    esperado = {
        6: 3.3740313741,
        8: 3.3740313858,
        16: 3.3740313858,
    }
    for n, value in esperado.items():
        spec = solve_spectrum(medium, build_quadrature(n))
        check(f"N={n}: nu1 = {value}", close(spec.modes[0].nu.real, value, 1e-9), f"{spec.modes[0].nu.real:.10f}")
    reference = continuum_nu(0.97)
    check("continuo nu1 = 3.3740313858", close(reference, 3.3740313858, 1e-9), f"{reference:.10f}")
    check(
        "error relativo S4 ~ 6.5e-6",
        close(abs(3.374009559245682 - reference) / reference, 6.47e-6, 1e-8),
    )


def test_invariantes() -> None:
    print("invariantes")
    medium = Medium.from_c(0.97)
    for n in (2, 4, 6, 8, 12, 16):
        spec = solve_spectrum(medium, build_quadrature(n))
        values = checks(spec)
        check(
            f"N={n}: norm/disp/mirror < 1e-10",
            max(values.values()) < 1e-10,
            f"max={max(values.values()):.1e}",
        )
        check(f"N={n}: 2N/2 modos", len(spec.modes) == n)


def test_multiplicativo() -> None:
    print("medio multiplicativo")
    medium = Medium(
        sigma_t=1.0, sigma_s0=0.5, sigma_f=0.02, nu_bar=2.4, k_eff=1.0
    )
    check("c = 1.1031858", close(medium.c, 1.1031857894892403, 1e-12), f"c={medium.c:.7f}")
    check("factor por defecto = 4*pi", close(medium.fission_factor, FISSION_FACTOR, 0.0))

    quad = build_quadrature(2)
    spec = solve_spectrum(medium, quad)
    omega_expected = 0.5773502691896257 / math.sqrt(medium.c - 1.0)
    check("S2: fundamental compleja +i*omega", abs(spec.modes[0].nu.imag - omega_expected) < 1e-10)
    check("S2: parte real nula", abs(spec.modes[0].nu.real) < 1e-12)
    check("S2: par conjugado", abs(spec.modes[1].nu + spec.modes[0].nu) < 1e-12)

    quad4 = build_quadrature(4)
    spec4 = solve_spectrum(medium, quad4)
    check("S4: fundamental compleja y un transitorio", len(spec4.modes) == 4)
    check("S4: transitorio real", abs(spec4.modes[2].nu.imag) < 1e-12)
    auto = eigen.path_a(medium, quad4)
    error = max(abs(a - b) for a, b in zip(auto, spec4.nu))
    check("camino A coincide con c>1", error < 1e-10, f"max|A-B|={error:.1e}")

    sin_fission = Medium.from_c(medium.c)
    check(
        "misma c, mismo espectro (sin fisión explícita)",
        max(
            abs(a - b)
            for a, b in zip(sin_fission.c and solve_spectrum(sin_fission, quad).nu, spec.nu)
        )
        < 1e-12,
    )

    # equivalencia de las dos ramas: sigma_s0 elegido para reproducir c = 0.97
    sigma_s0 = 0.97 - FISSION_FACTOR * 2.4 * 0.02
    equivalent = Medium(sigma_t=1.0, sigma_s0=sigma_s0, sigma_f=0.02, nu_bar=2.4, k_eff=1.0)
    check("rama multiplicativa reproduce c = 0.97", close(equivalent.c, 0.97, 1e-12))
    spec_eq = solve_spectrum(equivalent, build_quadrature(4))
    check(
        "mismo espectro que el caso de los documentos",
        max(abs(a - b) for a, b in zip(spec_eq.nu, build_quadrature(4) and solve_spectrum(Medium.from_c(0.97), build_quadrature(4)).nu))
        < 1e-12,
    )


def test_c_unidad() -> None:
    print("caso c = 1")
    medium = Medium.from_c(1.0)
    spec = solve_spectrum(medium, build_quadrature(4))
    check("sin fundamental finita", len(spec.modes) == 2)
    check("solo transitorios", all(mode.kind == "transitorio" for mode in spec.modes))
    check("aviso en el medio", any("c = 1" in w for w in medium.warnings()))


def test_errores() -> None:
    print("validaciones")
    try:
        Medium(sigma_t=1.0, sigma_f=0.02)
    except ValueError:
        check("sigma_f sin k_eff falla", True)
    else:
        check("sigma_f sin k_eff falla", False)
    try:
        Medium(sigma_t=0.0)
    except ValueError:
        check("sigma_t = 0 falla", True)
    else:
        check("sigma_t = 0 falla", False)
    try:
        build_quadrature(10, source="table")
    except ValueError:
        check("tabla sin N=10 falla", True)
    else:
        check("tabla sin N=10 falla", False)


def test_cli() -> None:
    print("CLI")
    buffer = io.StringIO()
    with redirect_stdout(buffer):
        code = main(["--n", "4", "--no-medium"])
    output = buffer.getvalue()
    check("exit code 0", code == 0)
    check("ecuación característica", "63 a^4 - 733 a^2 + 180 = 0" in output)
    check("nu1 en la salida", "3.3740096" in output)
    check("reporta los dos caminos", "camino A" in output)


def main_test() -> int:
    print("Regresión del análisis espectral S_N\n")
    test_quadrature()
    test_s2()
    test_s4()
    test_convergencia()
    test_invariantes()
    test_multiplicativo()
    test_c_unidad()
    test_errores()
    test_cli()
    print(f"\n{PASSED} verificaciones OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main_test())
