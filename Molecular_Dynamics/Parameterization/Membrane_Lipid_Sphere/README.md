Membrane Lipid Sphere Setup:

1)	Prepare Protein Structure (MODELLER)
a.	Check if crystal structure as any open ends
b.	Connect sequence gaps if possible (no need to model large unresolved domains, add few residues to obtain continuous structure.
2)	Convert .pdb to GROMACS (gmx pdb2gmx)
a.	Make sure you have the force field folder (charmm36-nov2018.ff) in your working directory
b.	gmx pdb2gmx -f model.pdb -ignh
i.	select correct force field and water model (TIP3P)
ii.	you could also select the protonation state of individual residues (probably not necessary at this stage)
c.	This should generate conf.gro, topol.top  posre.itp
d.	Get conf.pdb
i.	gmx editconf -f conf.gro -o conf.pdb
e.	Copy topol.top to Protein-atomistic.itp:
 		In Protein-atomistic.itp delete header until  
		[moleculetype] 
		and everything below:
		;Include Position restraint file 
		#ifdef POSRES 
		#include "posre.itp"
		#endif

3)	Convert Protein structure to coarse-grained (CG) MARTINIv2.2 representation (martinize.py)
a.	Be sure you have the correct secondary structure assignment using DSSP
b.	./martinize.py -f conf.gro -dssp /anaconda3/bin/mkdssp -x prot-cg.pdb -o prot-cg.top -ff martini22 -elastic -ep 6 -ea 3 -ef 500 -eu 0.9 -v -p ALL -cys auto -nt
c.	This will generate a CG protein structure (prot-cg.pdb), corase-grained topology (prot-cg.top), and Protein_A.itp. 
d.	CHECK IF YOUR PROTEIN IS COMPLETE!!! Possibly the last residue is missing. If that is the case change your conf.gro: Delete OT2 from C-term and rename OT1 to O, fix number of atoms in .gro header. 


4)	Set up Membrane (insane.py)
a.	./insane.py -f prot-cg.pdb -pbc rectangular -l POPC -center -o cg-membrane.gro -sol W -d 6 -dz 4
b.	Check terminal output for number of lipids and water added.
c.	Inspect cg-membrane.gro with pymol. Is the protein oriented correctly? 
(gmx editconf -f prot-cg.pdb -rotate x y z -o rot-prot-cg.pdb, use rot-prot-cg.pdb for input in insane.py)
d.	Prepare new topol.top
i.	Use topol-cg-temp and adjust number of lipid molecules. Save as topol.top

e.	Copy martini force field files: 
i.	martini_v2.2.itp
ii.	martini_v2.0_lipids_all.itp



5)	Minimize CG System 
a.	Cp martini_new-rf_min.mdp and martini_v2.x_new-rf.mdp
b.	Generate index.ndx file
i.	Gmx make_ndx -f cg-membrane.gro -o index.ndx
ii.	Select Protein and POPC (1 | 13; q)
c.	gmx grompp -f martini_new-rf_min.mdp -c cg-membrane.gro -p topol.top -n index.ndx -o min.tpr -po min-out.mdp
i.	This will generate a MD input file (min.tpr)
d.	gmx mdrun -v -deffnm min
i.	This will run a simulation and generate various outputs such as min.edr, min.gro, min.log, min.trr



6)	 Simulate CG system with position restraints on protein
a.	gmx grompp -f martini_v2.x_new-rf.mdp -c min.gro -p topol.top -n index.ndx -r min.gro -o md.tpr -po md-out.mdp
i.	Generates md.tpr
b.	Gmx mdrun -v -deffnm md
i.	Generates md.cpt, md.edr, md.gro, md.log, md.xtc, md_prev.cpt


7)	 Backmapping to atomistic resolution (initram.sh, backward.py)
a.	mkdir Backmap, cd
b.	Make molecules whole across PBC
i.	gmx trjconv -f md.gro -s md.tpr -o md-mol.pdb -pbc mol (echo 0)
c.	Generate atomistic topology
i.	Cp topol-aa-temp
ii.	Adjust number of lipids and water, save topol-aa.top
d.	Cp -r charmm36-nov2018.ff, popc.itp, posre_popc.itp, posre.itp, Protein-atomistic.itp
e.	./initram.sh -f md-mol.pdb -p top-aa.top -keep -to charmm36 
i.	You may have to try multiple times! It’s only done if you get until 6-mdpr-0.002.gro

8)	 Fit crystal-structure/prepared compete structure on backmapped protein in PyMOL
a.	Make backmapped.gro whole across PBC
i.	gmx trjconv -f backmapped.gro -s 6-mdpr-0.002.tpr -pbc mol -o backmapped-mol.pdb (echo 0)
ii.	In pymol align conf.pdb onto backmapped-mol.pdb
iii.	Remove protein in backmapped-mol.pdb
iv.	set retain_order, 1
v.	save molecule backmapped-mol.pdb as backmapped-environment.pdb
vi.	save fitted conf.pdb as fitted.pdb
vii.	cat fitted.pdb backampped-environment.pdb > fitted-system.pdb
viii.	Remove any END or TERM in fitted-system.pdb between protein and environment
ix.	Add box dimensions to header in fitted-system.pdb from backmapped-mol.pdb
x.	gmx editconf -f fitted-system.pdb -o fitted-renumber.pdb


9)	 Run minimization with frozen protein
a.	Cp min_freeze.mdp
b.	Generate index_freeze.ndx
i.	gmx make_ndx -f fitted-system.pdb -o index_freeze.ndx
ii.	Group Membrane and SOL, name group Rest
c.	gmx grompp -f min_freeze.mdp -c fitted-renumber.pdb -n index_freeze.ndx -p backmapped.top -po min_freeze-out.mdp -o min_freeze.tpr 
d.	gmx mdrun -v -deffnm min_freeze
e.	gmx grompp -f min.mdp -c min_freeze.gro -n index.ndx -p backmapped.top -po min-out.mdp -o min.tpr
f.	gmx mdrun -v -deffnm min


10)	Select lipid spheres
a.	gmx trjconv -f min.gro -s min.tpr -pbc mol -o min-mol.pdb -p backmapped.top 
b.	run prepare.pml
