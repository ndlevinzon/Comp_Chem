import pytest

pytest.skip("legacy functional tests - skip during automated test run", allow_module_level=True)

import tempfile
import sys, os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from grocharmm.gro_to_crd import replace_resnames, update_segment, update_segment_number

def test_replace_resnames_basic():
    input_content = "TRI DOP POP SWM\nATOM 1 TRI 1\nATOM 2 POP 2"
    expected_output = "TRIO DOPE POPC SWM4\nATOM 1 TRIO 1\nATOM 2 POPC 2"
    with tempfile.NamedTemporaryFile("w+", delete=False) as temp_in:
        temp_in.write(input_content)
        temp_in.seek(0)
        input_path = temp_in.name
    with tempfile.NamedTemporaryFile("r", delete=False) as temp_out:
        output_path = temp_out.name
    replace_resnames(input_path, output_path)
    with open(output_path, 'r') as f:
        result = f.read()

    os.unlink(input_path)
    os.unlink(output_path)
    assert result == expected_output



def test_update_segment():
    input_content = (
        "* Comment line\n"
        "         1         1  TIP3      O      0.000   0.000   0.000  SYSTEM    1     0.0000000000\n"
        "         2         2  TRIO      N      0.000   0.000   0.000  SYSTEM    1     0.0000000000\n"
        "         3         3  ALA       C      0.000   0.000   0.000  SYSTEM    1     0.0000000000\n"
        "         4         4  TIP3      O      0.000   0.000   0.000  MEMB      1     0.0000000000\n"
    )
    expected_output = (
        "* Comment line\n"
        "         1         1  TIP3      O      0.000   0.000   0.000  TIP3      1     0.0000000000\n"
        "         2         2  TRIO      N      0.000   0.000   0.000  MEMB      1     0.0000000000\n"
        "         3         3  ALA       C      0.000   0.000   0.000  PROA      1     0.0000000000\n"
        "         4         4  TIP3      O      0.000   0.000   0.000  MEMB      1     0.0000000000\n"
    )

    with tempfile.NamedTemporaryFile(mode="w+", delete=False) as temp_in:
        temp_in.write(input_content)
        temp_in.seek(0)
        input_path = temp_in.name
    with tempfile.NamedTemporaryFile(mode="r+", delete=False) as temp_out:
        output_path = temp_out.name
    update_segment(input_path, output_path)
    with open(output_path, "r") as f:
        output = f.read()
    os.remove(input_path)
    os.remove(output_path)
    assert output == expected_output


