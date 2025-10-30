# GRO to CHARMM-Drude File Editor (GUI)

A lightweight Python GUI for converting `.gro` files into CHARMM-compatible `.crd` and `.psf` formats. Designed for Drude and other CHARMM-based workflows.

## Features

- Convert `.gro` to `.crd` and `.psf` 
- Optional residue renaming (CRD mode only)
- Automatic segment name and residue number updates
- Custom segment rule editor

## Requirements

- Python 3.10 or later
- Parmed
- `tkinter` (bundled with Python on Windows; may need to be installed on Linux/macOS)

## Installation

Install the package and its command-line entry point:

```bash
pip install .
```

## Usage

After installation, launch the GUI with:

```bash
grocharmm-gui
```

### Steps

1. Select an input `.gro` file.
2. Choose an output path and format (`.crd` or `.psf`).
3. Select mode: `CRD` or `PSF`.
4. (Optional) Enable "Replace resnames" (CRD only).
5. Use **Preview Segments** to automatically populate segment rules.
6. Modify segment rules if needed (format: `RESNAME:SEGNAME`).
7. Click **Run** to generate the output file.

## Default Segment Rules

The following rules are used unless you override them:

```
TIP3:TIP3
SOD:IONS
CLA:IONS
DOPE:MEMB
POPC:MEMB
TRIO:MEMB
```

## Output

The output `.crd` or `.psf` file will be saved to the selected location. Existing files will be overwritten.



## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

© 2025 Jay Braun
