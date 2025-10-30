import sys
import types
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
import tempfile

# Ensure gro_to_psf can be imported even if parmed is missing
try:
    import parmed  # noqa: F401
except ImportError:  # pragma: no cover
    sys.modules['parmed'] = types.ModuleType('parmed')

from grocharmm.gro_to_psf import PSFeditor


def _make_line(atom_id: int, resid: int, resname: str, atomname: str) -> str:
    """Return a single PSF atom line with predictable spacing."""
    return (
        f"{atom_id:>10} SEG1     {resid:<9}{resname:<6}  {atomname:<6}A     0.000000  1.0000\n"
    )


def _write_sample_psf() -> str:
    """Create a temporary PSF file with a minimal !NATOM section."""
    lines = [
        _make_line(1, 1, "DOPE", "CA"),
        _make_line(2, 2, "SOD", "NA"),
        _make_line(3, 3, "TIP3", "O"),
        _make_line(4, 3, "TIP3", "H1"),
        _make_line(5, 3, "TIP3", "H2"),
        _make_line(6, 3, "TIP3", "H3"),
        _make_line(7, 4, "GLY", "N"),
        _make_line(8, 5, "SOD", "NA"),
    ]
    content = "PSF\n\n   1 !NTITLE\n\n       8 !NATOM\n" + "".join(lines) + "!BOND\n"
    tmp = tempfile.NamedTemporaryFile("w+", delete=False)
    tmp.write(content)
    tmp.close()
    return tmp.name


def test_update_psf_segments_basic():
    path = _write_sample_psf()
    editor = PSFeditor(path)
    editor.read_lines()
    editor.update_psf_segments()

    natom_index = next(i for i, l in enumerate(editor.lines) if "!NATOM" in l)
    atoms = editor.lines[natom_index + 1 : natom_index + 9]
    segments = [line[11:16].strip() for line in atoms]

    assert segments == [
        "MEMB",
        "IONS",
        "TIP3",
        "TIP3",
        "TIP3",
        "TIP3",
        "PROA",
        "IONS",
    ]


def test_update_psf_resids_tip3_handling():
    path = _write_sample_psf()
    editor = PSFeditor(path)
    editor.read_lines()
    editor.update_psf_segments()
    editor.update_psf_resids()

    natom_index = next(i for i, l in enumerate(editor.lines) if "!NATOM" in l)
    atoms = editor.lines[natom_index + 1 : natom_index + 9]
    resids = [int(line[20:29].strip()) for line in atoms]

    assert resids == [1, 1, 1, 1, 1, 2, 1, 1]
