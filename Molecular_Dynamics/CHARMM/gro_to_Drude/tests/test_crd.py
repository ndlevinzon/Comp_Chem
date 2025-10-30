import pytest

pytest.skip("legacy manual script - skip during automated test run", allow_module_level=True)

import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from grocharmm.gro_to_crd import CRDeditor

editor = CRDeditor('data/tld1.crd')
editor.read()
editor.update_segments()
editor.update_segment_numbers()
editor.write('data/tld1_updated.crd')

print("Output saved to test/tld1_updated.crd")

