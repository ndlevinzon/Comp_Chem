
DEFAULT_SEGMENT_RULES = {
    'TIP3': 'TIP3',
    'SOD': 'IONS',
    'CLA': 'IONS',
    'DOPE': 'MEMB',
    'POPC': 'MEMB',
    'TRIO': 'MEMB',
}

def assign_segment(resname, rules=None):
    rules = rules or DEFAULT_SEGMENT_RULES
    return rules.get(resname, 'PROA')[:5]


def extract_segs_crd(file_path):
    """Extract unique segment names from a CHARMM .crd file."""
    segments = set()
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith("*") or "EXT" in line or not line.strip():
                continue
            if len(line) >= 108:
                segment = line[102:108].strip()
                if segment:
                    segments.add(segment)
    return segments


def extract_resnames_crd(crd_file):
    resnames = set()
    with open(crd_file, 'r') as f:
        for line in f:
            if len(line) >= 26 and not line.startswith('*') and line.strip():
                resname = line[20:26].strip()
                if resname:
                    resnames.add(resname)
    return resnames
