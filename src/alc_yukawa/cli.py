from __future__ import annotations

import argparse
from pathlib import Path
import sys

from .pipeline import run_fit
from .charm_study import run_charm_study
from .global_study import run_global_study


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="alc-yukawa",
        description="ALC-inspired Yukawa ladder test using PDG public masses.",
    )
    p.add_argument(
        "--config",
        type=str,
        required=True,
        help="Path to preregistered YAML config (e.g. config/preregistered.yaml)",
    )
    sub = p.add_subparsers(dest="cmd", required=True)

    fit = sub.add_parser("fit", help="Run the ladder fit and write artifacts.")
    fit.add_argument(
        "--output-dir",
        type=str,
        default=None,
        help="Override output dir (otherwise uses config).",
    )

    charm = sub.add_parser("charm-study", help="Run the charm robustness study.")
    charm.add_argument(
        "--output-dir",
        type=str,
        default=None,
        help="Override output dir (otherwise uses config).",
    )

    glob = sub.add_parser("global-study", help="Run the global common-scale alpha-ladder test.")
    glob.add_argument(
        "--output-dir",
        type=str,
        default=None,
        help="Override output dir (otherwise uses config).",
    )

    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    cfg_path = Path(args.config)
    if not cfg_path.exists():
        print(f"Config not found: {cfg_path}", file=sys.stderr)
        return 2

    if args.cmd == "fit":
        out_dir = Path(args.output_dir) if args.output_dir else None
        run_fit(cfg_path=cfg_path, output_dir_override=out_dir)
        return 0

    if args.cmd == "charm-study":
        out_dir = Path(args.output_dir) if args.output_dir else None
        run_charm_study(cfg_path=cfg_path, output_dir_override=out_dir)
        return 0

    if args.cmd == "global-study":
        out_dir = Path(args.output_dir) if args.output_dir else None
        run_global_study(cfg_path=cfg_path, output_dir_override=out_dir)
        return 0

    return 1
