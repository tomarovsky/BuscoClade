# Unit tests

Pytest unit tests for the pure-logic helpers in `workflow/scripts/`. These test
the Python scripts directly (parsing, coordinate math, sequence transforms) —
they do **not** run Snakemake or any external tool.

## Running

Run from the repository root, inside the `buscoclade_main` conda env:

```bash
conda activate buscoclade_main
pytest tests/unit/
```

Useful variants:

```bash
pytest tests/unit/ -v                                  # one line per test
pytest tests/unit/test_quastcore.py                    # a single file
pytest tests/unit/ -k apply_snps                        # tests matching a keyword
```

## How imports work

`workflow/scripts/` is not an installed package — the scripts are executed
directly by the pipeline. `conftest.py` puts that directory on `sys.path` so a
test can do `import busco_reconstruct_common` (etc.) by module name. This is why
tests must be run from the repo root.

## Layout

- `conftest.py` — adds `workflow/scripts/` to the import path.
- `test_<script>.py` — one file per script under test.

Tests that need heavier dependencies guard the import with
`pytest.importorskip` (`pysam` for the reconstruction scripts, `ete3` for
`astral_metrics`), so they skip rather than error when a dependency is absent.

Pure-plotting scripts (`draw_phylotrees*`, `busco_histogram`,
`unique_ids_plot`) and the vendored third-party `vcf2phylip.py` are intentionally
not covered here.
