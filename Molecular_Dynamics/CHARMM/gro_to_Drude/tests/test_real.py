import pytest

pytest.skip("legacy integration script - skip during test run", allow_module_level=True)

import shutil
import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from grocharmm.gro_to_crd import replace_resnames, update_segment, update_segment_number

def run_replace_real_crd_file():
    input_real = "../original_psf_crd/tld1.crd"
    temp_file = "../original_psf_crd/tld1_test_temp.crd"

    shutil.copy(input_real, temp_file)
    print(f"Copied to: {temp_file}")

    replace_resnames(temp_file)

    print(f"Done! Modified file saved at: {temp_file}")


def run_update_segment():
    input_real = "../original_psf_crd/tld1.crd"
    temp_file = "../original_psf_crd/tld1_seg_temp.crd"

    shutil.copy(input_real, temp_file)
    print(f"Copied to: {temp_file}")

    update_segment(input_real,temp_file)

    print(f"Done! Modified file saved at: {temp_file}")


def run_update_segment_number():
    input_real = "../original_psf_crd/tld1_seg_temp.crd"
    temp_file = "../original_psf_crd/tld1_renumberd_temp.crd"

    shutil.copy(input_real, temp_file)
    print(f"Copied to: {temp_file}")

    update_segment_number(input_real,temp_file)

    print(f"Done! Modified file saved at: {temp_file}")






if __name__ == "__main__":
    run_update_segment_number()
