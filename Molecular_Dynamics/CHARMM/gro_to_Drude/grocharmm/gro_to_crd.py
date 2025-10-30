from grocharmm.utils import assign_segment


class CRDeditor:
    def __init__(self, input_file):
        self.input_file = input_file
        self.lines = []

    def read(self):
        with open(self.input_file, 'r') as f:
            self.lines = f.readlines()

    def write(self, output_file):
        with open(output_file, 'w') as f:
            f.writelines(self.lines)

    def replace_resnames(self, replacements=None):
        if replacements is None:
            replacements = {
                'TRI': 'TRIO',
                'POP': 'POPC',
                'DOP': 'DOPE',
                'SWM': 'SWM4'
            }

        content = ''.join(self.lines)
        for old, new in replacements.items():
            content = content.replace(old, new)
        self.lines = content.splitlines(keepends=True)

    def update_segments(self, segment_rules=None):
        # segment_rules = {
        #     'TIP3': 'TIP3',
        #     'SOD': 'IONS',
        #     'CLA': 'IONS',
        #     'DOPE': 'MEMB',
        #     'POPC': 'MEMB',
        #     'TRIO': 'MEMB',
        # }

        # def assign_segment(resname):
        #     return segment_rules.get(resname, 'PROA')

        updated_lines = []
        for line in self.lines:
            if 'SYSTEM' not in line:
                updated_lines.append(line)
                continue

            resname = line[20:26].strip()
            new_segment = assign_segment(resname, segment_rules)
            updated_line = line.replace('SYSTEM', new_segment.ljust(6), 1)
            updated_lines.append(updated_line)

        self.lines = updated_lines

    def update_segment_numbers(self):
        current_segment = None
        current_resid = None
        counter = 0
        updated_lines = []

        for line in self.lines:
            if line.startswith("*") or "EXT" in line or not line.strip():
                updated_lines.append(line)
                continue

            segment = line[102:108].strip()
            resid = int(line[10:20].strip())

            if segment != current_segment:
                counter = 1
                current_segment = segment
                current_resid = resid
            elif resid != current_resid:
                counter += 1
                current_resid = resid

            updated_line = line[:112] + f"{counter:<16}" + line[128:]
            updated_lines.append(updated_line)

        self.lines = updated_lines
