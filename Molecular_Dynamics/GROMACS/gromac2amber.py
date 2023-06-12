import parmed as pmd

# convert GROMACS topology to AMBER format


IN  = 'file.rst'
OUT = 'file.gro'

amb_top = pmd.load_file('file.top', xyz=IN)
amb_top.save(OUT)
