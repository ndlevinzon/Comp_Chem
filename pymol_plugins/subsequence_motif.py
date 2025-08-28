from pymol import cmd
import re

AA3_TO_1 = {
 'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C',
 'GLN':'Q','GLU':'E','GLY':'G','HIS':'H','ILE':'I',
 'LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P',
 'SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V',
 'MSE':'M','SEC':'U','PYL':'O'
}

def _seq_and_reslist(sel):
    """Return (sequence_string, [(segi, chain, resi)]) for CA-ordered residues."""
    model = cmd.get_model(f"({sel}) and polymer.protein and name CA")
    seq, reskeys, seen = [], [], set()
    for a in model.atom:
        key = (a.segi, a.chain, a.resi)
        if key in seen: 
            continue
        seen.add(key)
        seq.append(AA3_TO_1.get(a.resn.upper(), 'X'))
        reskeys.append(key)
    return ''.join(seq), reskeys

def find_subseq(obj, subseq, chain=None, name_prefix="motif"):
    """Exact (non-regex) subsequence search, e.g. find_subseq 1abc, RGD, chain=A"""
    sel_base = obj if chain is None else f"{obj} and chain {chain}"
    seq, reskeys = _seq_and_reslist(sel_base)
    pat = subseq.upper()
    i, hits = 0, []
    while True:
        j = seq.find(pat, i)
        if j == -1: break
        hits.append((j, j+len(pat)))
        i = j+1
    if not hits:
        print(f"No hits for '{subseq}' in {sel_base}")
        return []
    out = []
    for n, (s, e) in enumerate(hits, 1):
        keys = reskeys[s:e]
        by_chain = {}
        for segi, ch, r in keys:
            by_chain.setdefault(ch, []).append(r)
        parts = [f"(chain {ch} and resi {'+'.join(rlist)})" for ch, rlist in by_chain.items()]
        selection = f"({obj}) and polymer.protein and (" + " or ".join(parts) + ")"
        name = f"{name_prefix}_{n}"
        cmd.select(name, selection)
        print(f"{name}: {selection}")
        out.append(name)
    return out

def find_subseq_regex(obj, regex, chain=None, name_prefix="motif"):
    """Regex search (e.g., N[^P][ST] for N-X-S/T). Use like: find_subseq_regex 1abc, N[^P][ST], chain=A"""
    sel_base = obj if chain is None else f"{obj} and chain {chain}"
    seq, reskeys = _seq_and_reslist(sel_base)
    pat = re.compile(regex.upper())
    out = []
    for n, m in enumerate(pat.finditer(seq), 1):
        s, e = m.span()
        keys = reskeys[s:e]
        by_chain = {}
        for segi, ch, r in keys:
            by_chain.setdefault(ch, []).append(r)
        parts = [f"(chain {ch} and resi {'+'.join(rlist)})" for ch, rlist in by_chain.items()]
        selection = f"({obj}) and polymer.protein and (" + " or ".join(parts) + ")"
        name = f"{name_prefix}_{n}"
        cmd.select(name, selection)
        print(f"{name}: {selection}")
        out.append(name)
    if not out:
        print(f"No regex hits for '{regex}' in {sel_base}")
    return out

cmd.extend("find_subseq", find_subseq)
cmd.extend("find_subseq_regex", find_subseq_regex)
