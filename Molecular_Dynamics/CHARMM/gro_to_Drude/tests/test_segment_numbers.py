import os
import sys
import tempfile

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from grocharmm.gro_to_crd import CRDeditor


BASE_LINE = (
    "         2         1  MET       HY1            38.9799995422       68.9799957275"
    "      115.1799926758  PROA      1               0.0000000000"
)


def _make_line(atom_id: int, resid: int, segment: str) -> str:
    """Return a single CRD line with the given values inserted."""
    line = BASE_LINE
    line = f"{atom_id:>10}" + line[10:]
    line = line[:10] + f"{resid:>10}" + line[20:]
    line = line[:102] + f"{segment:<6}" + line[108:]
    return line


def test_update_segment_numbers_basic():
    """CRDeditor.update_segment_numbers should renumber residues per segment."""

    lines = [
        _make_line(1, 1, "SEGA"),
        _make_line(2, 2, "SEGA"),
        _make_line(3, 3, "SEGB"),
        _make_line(4, 4, "SEGB"),
    ]

    sample = "*\n" + "\n".join(lines) + "\n"

    with tempfile.NamedTemporaryFile("w+", delete=False) as tmp:
        tmp.write(sample)
        tmp_path = tmp.name

    editor = CRDeditor(tmp_path)
    editor.read()
    editor.update_segment_numbers()

    expected_numbers = [1, 2, 1, 2]
    expected_lines = [
        "*\n",
    ] + [
        line[:112] + f"{num:<16}" + line[128:] + "\n"
        for line, num in zip(lines, expected_numbers)
    ]

    assert editor.lines == expected_lines
