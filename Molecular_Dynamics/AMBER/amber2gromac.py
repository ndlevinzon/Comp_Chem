import parmed as pmd

# convert GROMACS topology to AMBER format

IN  = 'frame.gro'
OUT = 'frame.crd'
TOP = 'frame.top'

gmx_top = pmd.load_file('topol.top', xyz=IN)
gmx_top.save(TOP, format='amber')
gmx_top.save(OUT, format='rst7')
