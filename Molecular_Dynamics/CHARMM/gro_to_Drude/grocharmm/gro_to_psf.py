import parmed as pmd
from grocharmm.utils import assign_segment


class PSFeditor:
    def __init__(self, input_file):
        self.input_file = input_file

    
    def read_lines(self):
        with open(self.input_file,'r') as f:
            self.lines = f.readlines()

    def write_lines(self, output_path): 
        with open(output_path, 'w') as f:
            f.writelines(self.lines)
        
    def load_inp(self, topol_file='topol.top', crd_file='tld1_final.crd', output_psf='out1.psf'):
        topology = pmd.load_file(topol_file)
        coordinates = pmd.load_file(crd_file)
        topology.coordinates = coordinates.coordinates
        topology.save(output_psf, overwrite=True)

    def update_psf_segments(self, segment_rules=None):

        """
        Update segment names in the atom section of a CHARMM .psf file.
        Only modifies lines between !NATOM and the next line starting with '!'.
        """
        # segment_rules = {
        #     'TIP3': 'TIP3',
        #     'SOD': 'IONS',
        #     'CLA': 'IONS',
        #     'DOPE': 'MEMB',
        #     'POPC': 'MEMB',
        #     'TRIO': 'MEMB',
        # }

        # def assign_segment(resname, segment_rules=None):
        #     return segment_rules.get(resname, 'PROA')[:5]

        updated_lines = []
        in_atom_section = False
        atom_lines_remaining = 0

        # with open(input_path, 'r') as f:
        for line in self.lines:
            # for line in f:
            if '!NATOM' in line:
                in_atom_section = True
                try:
                    atom_lines_remaining = int(line.strip().split()[0])
                except Exception:
                    raise ValueError("Failed to read atom count from !NATOM line.")
                updated_lines.append(line)
                continue

            if in_atom_section and atom_lines_remaining > 0:
                if len(line) >= 35 and line[11:16].strip():
                    resname = line[29:35].strip()
                    new_segment = assign_segment(resname, segment_rules)
                    updated_line = line[:11] + f"{new_segment:<5}" + line[16:]
                    updated_lines.append(updated_line)
                else:
                    updated_lines.append(line)
                atom_lines_remaining -= 1
            else:
                updated_lines.append(line)

        self.lines = updated_lines


    def update_psf_resids(self):
        """
        Update residue numbers (column 3: chars 21–29) in PSF atom section lines.
        - When the segment changes, counter resets (except PROA, which continues).
        - For TIP3, increment every 3 atoms.
        - Operates strictly between !NATOM and the next section header (!...).
        """
        updated_lines = []
        in_atom_section = False
        atom_lines_remaining = 0
        counter = 0
        tip3_atom_count = 0
        last_segment = None
        last_resid = None

        for line in self.lines:
            # Detect start of atom section
            if '!NATOM' in line:
                in_atom_section = True
                atom_lines_remaining = int(line.strip().split()[0])
                updated_lines.append(line)
                continue

            # Only update during atom section
            if in_atom_section and atom_lines_remaining > 0:
                if len(line) < 50 or not line[11:16].strip():
                    updated_lines.append(line)
                    atom_lines_remaining -= 1
                    continue

                atom_index = int(line[:10])
                segment = line[10:20]     # col 2
                resid_str = line[20:30]   # col 3
                resid = resid_str.strip()
                segment_clean = segment.strip()

                # Logic for updating counter
                if segment_clean != last_segment:
                    counter = 1
                    tip3_atom_count = 1 if segment_clean == "TIP3" else 0
                else:
                    if segment_clean == "TIP3":
                        tip3_atom_count += 1
                        if tip3_atom_count > 3:
                            counter += 1
                            tip3_atom_count = 1
                    elif resid != last_resid:
                        counter += 1

                last_segment = segment_clean
                last_resid = resid

                updated_line = (
                    f"{atom_index:>10}"
                    f"{segment:<10}"
                    f"{counter:<9}"      # proper width: chars 21–29
                    + line[29:]
                )
                updated_lines.append(updated_line)
                atom_lines_remaining -= 1
            else:
                updated_lines.append(line)

        self.lines = updated_lines

