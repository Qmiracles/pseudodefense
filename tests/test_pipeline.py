from pathlib import Path

import pytest

from pseudodefense import pipeline


def test_default_command():
    assert pipeline.snakemake_command(4) == ["snakemake", "--cores", "4", "--rerun-incomplete"]


def test_overrides_dry_run_and_targets():
    cmd = pipeline.snakemake_command(
        8,
        config=["outdir=results/pilot", "n_genomes=10"],
        configfile=Path("my.yaml"),
        dry_run=True,
        targets=["report"],
    )
    assert cmd == [
        "snakemake", "--cores", "8", "--rerun-incomplete", "--configfile", "my.yaml",
        "--config", "outdir=results/pilot", "n_genomes=10", "-n", "report",
    ]  # fmt: skip


@pytest.mark.parametrize("bad", [{"config": ["outdir"]}, {"cores": 0}])
def test_invalid_arguments(bad):
    kwargs = {"cores": 2, **bad}
    with pytest.raises(ValueError):
        pipeline.snakemake_command(**kwargs)
