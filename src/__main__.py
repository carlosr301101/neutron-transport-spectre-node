"""CLI del análisis espectral S_N.

Ejemplos
--------
Caso de los documentos (sigma_t = 1, sigma_s0 = 0.97, sin fisión)::

    python -m codigo --n 2,4

Medio multiplicativo (c = 1.1032, fundamental compleja)::

    python -m codigo --n 4 --sigma-t 1 --sigma-s 0.5 \
        --sigma-f 0.02 --nu-bar 2.4 --k-eff 1

Barrido de convergencia::

    python -m codigo --n 2,4,6,8,16 --compare --format md --out reporte.md
"""

from __future__ import annotations

import argparse
import csv
import json
import pathlib
import sys

from . import eigen
from .medium import FISSION_FACTOR, Medium
from .quadrature import DEFAULT_TABLE, build_quadrature
from .report import convergence_table, full_report
from .spectrum import checks, solve_spectrum


def _parse_orders(text: str) -> list[int]:
    orders: list[int] = []
    for chunk in text.split(","):
        chunk = chunk.strip()
        if not chunk:
            continue
        try:
            value = int(chunk)
        except ValueError as exc:
            raise ValueError(f"orden inválido: {chunk!r}") from exc
        if value < 2 or value % 2 != 0:
            raise ValueError(f"S_N requiere N par y >= 2 (recibido {value})")
        if value not in orders:
            orders.append(value)
    if not orders:
        raise ValueError("hay que indicar al menos un orden N")
    return orders


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="python -m codigo",
        description="Autovalores y autofunciones del problema espectral S_N en placa 1D "
        "por dos caminos (matricial y relación de dispersión).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--n", default="2,4", help="órdenes S_N separados por coma (pares)")
    parser.add_argument("--sigma-t", type=float, default=1.0, dest="sigma_t",
                        help="sección eficaz total (cm^-1)")
    parser.add_argument("--sigma-s", type=float, default=0.97, dest="sigma_s0",
                        help="sección eficaz de dispersión isotrópica (cm^-1)")
    parser.add_argument("--sigma-f", type=float, default=0.0, dest="sigma_f",
                        help="sección eficaz de fisión (cm^-1); activa la rama multiplicativa")
    parser.add_argument("--nu-bar", type=float, default=0.0, dest="nu_bar",
                        help="neutrones por fisión")
    parser.add_argument("--k-eff", type=float, default=None, dest="k_eff",
                        help="factor de multiplicación (obligatorio si --sigma-f > 0)")
    parser.add_argument("--fission-factor", type=float, default=FISSION_FACTOR,
                        dest="fission_factor",
                        help=f"peso del término de fisión (por defecto 4*pi = {FISSION_FACTOR:.6f})")
    parser.add_argument("--c", type=float, default=None,
                        help="impone c directamente (medio no multiplicativo; ignora sigma_s y fisión)")
    parser.add_argument("--quadrature", choices=["auto", "table", "generated"], default="auto",
                        help="origen de los nodos/pesos")
    parser.add_argument("--quadrature-file", default=str(DEFAULT_TABLE),
                        help="ruta a cuadraturas.py")
    parser.add_argument("--format", choices=["text", "md"], default="text",
                        help="formato de las tablas")
    parser.add_argument("--out", default=None, help="archivo de salida (por defecto stdout)")
    parser.add_argument("--json", default=None, help="exporta el espectro completo a JSON")
    parser.add_argument("--csv", default=None, help="exporta los autovalores a CSV")
    parser.add_argument("--compare", action="store_true",
                        help="incluye la tabla de convergencia de nu_1 (se fuerza si hay más de un N)")
    parser.add_argument("--no-medium", action="store_true",
                        help="omite el bloque del medio en cada orden")
    return parser


def _make_medium(args: argparse.Namespace) -> Medium:
    if args.c is not None:
        return Medium.from_c(args.c, sigma_t=args.sigma_t)
    return Medium(
        sigma_t=args.sigma_t,
        sigma_s0=args.sigma_s0,
        sigma_f=args.sigma_f,
        nu_bar=args.nu_bar,
        k_eff=args.k_eff,
        fission_factor=args.fission_factor,
    )


def _jsonable(z: complex):
    return [z.real, z.imag]


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    markdown = args.format == "md"

    try:
        medium = _make_medium(args)
        orders = _parse_orders(args.n)
        quads = [
            build_quadrature(order, source=args.quadrature, path=args.quadrature_file)
            for order in orders
        ]

        chunks: list[str] = []
        payload: dict = {"medium": medium.to_dict(), "orders": []}

        for index, (quad, order) in enumerate(zip(quads, orders)):
            spectrum = solve_spectrum(medium, quad)
            nus_a = eigen.path_a(medium, quad)
            residual = eigen.matrix_residual(medium, quad, spectrum)

            chunks.append(
                full_report(
                    medium,
                    quad,
                    spectrum,
                    nus_a,
                    residual,
                    markdown=markdown,
                    include_medium=not (args.no_medium or index > 0),
                )
            )
            payload["orders"].append(
                {
                    "n": quad.n,
                    "quadrature_source": quad.source,
                    "quadrature": {
                        "mu": quad.mu,
                        "w": quad.w,
                        "total_weight": quad.total_weight,
                    },
                    "characteristic_polynomial": [float(c) for c in spectrum.poly],
                    "characteristic_polynomial_exact": (
                        [int(c) for c in spectrum.poly_exact]
                        if spectrum.poly_exact is not None
                        else None
                    ),
                    "modes": [
                        {
                            "j": mode.index,
                            "nu": _jsonable(mode.nu),
                            "lambda_1_over_nu": _jsonable(1.0 / mode.nu),
                            "kind": mode.kind,
                            "v": [_jsonable(v) for v in mode.v],
                        }
                        for mode in spectrum.modes
                    ],
                    "path_a": [_jsonable(z) for z in nus_a],
                    "checks": dict(residual_matrix=residual, **checks(spectrum)),
                }
            )

        if args.compare or len(quads) > 1:
            title = "Convergencia de la fundamental"
            chunks.append(f"## {title}" if markdown else f"=== {title} ===")
            chunks.append(convergence_table(medium, quads, markdown))

        text = "\n\n".join(chunks) + "\n"
        if args.out:
            pathlib.Path(args.out).write_text(text, encoding="utf-8")
            print(f"reporte escrito en {args.out}")
        else:
            print(text)

        if args.json:
            pathlib.Path(args.json).write_text(
                json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8"
            )
            print(f"JSON escrito en {args.json}")

        if args.csv:
            with open(args.csv, "w", newline="", encoding="utf-8") as handle:
                writer = csv.writer(handle)
                writer.writerow(["n", "c", "j", "nu_real", "nu_imag", "decay", "kind"])
                for entry in payload["orders"]:
                    for mode in entry["modes"]:
                        nu_re, nu_im = mode["nu"]
                        decay = 1.0 / nu_re if abs(nu_im) < 1e-12 and nu_re != 0.0 else ""
                        writer.writerow(
                            [entry["n"], medium.c, mode["j"], nu_re, nu_im, decay, mode["kind"]]
                        )
            print(f"CSV escrito en {args.csv}")

    except (ValueError, RuntimeError, FileNotFoundError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2

    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
