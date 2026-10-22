"""Thin wrapper to launch the Snakemake workflow from the ``pseudodefense`` CLI."""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path


def snakemake_command(
    cores: int,
    config: Sequence[str] = (),
    configfile: Path | None = None,
    dry_run: bool = False,
    targets: Sequence[str] = (),
) -> list[str]:
    """Build the ``snakemake`` command line for a pipeline run.

    Args:
        cores: CPU cores Snakemake may use.
        config: ``KEY=VALUE`` overrides of ``config/config.yaml``.
        configfile: alternative config file (replaces ``config/config.yaml``).
        dry_run: only print the plan (``-n``).
        targets: rules or files to build (default: rule ``all``).

    """
    bad = [c for c in config if "=" not in c]
    if bad:
        raise ValueError(f"config overrides must be KEY=VALUE, got: {bad}")
    if cores < 1:
        raise ValueError("cores must be >= 1")
    cmd = ["snakemake", "--cores", str(cores), "--rerun-incomplete"]
    if configfile is not None:
        cmd += ["--configfile", str(configfile)]
    if config:
        cmd += ["--config", *config]
    if dry_run:
        cmd.append("-n")
    cmd += list(targets)
    return cmd
