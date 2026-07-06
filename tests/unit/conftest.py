"""Make workflow/scripts importable so the shared reconstruction helpers can be
unit-tested directly (they are executed as scripts in the pipeline, so they are
not an installed package)."""
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "workflow" / "scripts"))
