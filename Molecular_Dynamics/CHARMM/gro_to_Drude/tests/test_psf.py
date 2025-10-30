import pytest

pytest.importorskip("parmed", reason="parmed not installed")

# The original file contained a manual example that required the `parmed`
# package.  Without the dependency available, importing the module would fail
# during test collection.  The lines below are kept for reference but the test
# is skipped when the dependency is missing.

import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from grocharmm.gro_to_psf import PSFeditor

editor = PSFeditor("data/tld1_generated.psf")
editor.load_inp(
    topol_file="data/topol.top",
    crd_file="data/tld1_final.crd",
    output_psf="data/tld1_generated.psf",
)

editor.read_lines()
editor.update_psf_segments()
editor.update_psf_resids()
editor.write_lines("data/output.psf")

print("Final PSF saved to test/output.psf")
