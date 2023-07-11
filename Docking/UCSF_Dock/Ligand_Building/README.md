# The DB2 File Format and Building Ligands
## The Brass Tacks

DB2 files encode & compress 3D molecule rotomers, sometimes referred to as conformations. DB2 files do not store the dihedral angles of rotomers, but instead the per-atom coordinates resulting from those angles. Atoms of different rotomers whose coordinates overlap are merged together in db2 files (how this is done will be detailed more later), thus allowing DOCK to avoid double-calculating per-atom energies.
## The nitty gritty

DB2 is a text format, with each row of information being printed on a new line (not to exceed 80 characters per line, due to old fortran nonsense). Mutiple DB2 entries can be stored in a db2 file- there is no limit on the number of entries that can be fit in a single db2 file. Thus concatenating any number of db2 files together still results in a valid db2 file.

There are seven(/eight) species of lines in a db2 entry- they will be listed in the order they appear.

* "M" lines serve a few functions. First is to define information about a db2 molecule entry- molecule name, smiles, properties, etc. The second is to create a boundary between distinct db2 entries

* "A" and "B" lines define atoms and bonds for a molecule. These are practically identical to atom/bond lines in .mol2 files- not much else to say here

* "X" lines define the set of all coordinates for the entry- these will be referenced later on. coordinates are guaranteed to be "distinct", in the sense that no coordinate is within a threshold distance (typically 0.001A) of any other coordinate in the set of all coordinates

* "R" lines define the coordinates of the rigid component. this is the basis component of the molecule that does not move between rotomers, usually a benzene ring or the like. the choice of which structure should be the rigid component is arbitrary, and in fact usually a separate db2 entry is created for each possible rigid component (usually each ring system) in a molecule

* "C" lines are "confs", not to be confused with conformations. We will get back to these.

* "S" lines or "set" lines define the rotomers/conformations. Each set entry may span multiple "S" lines to fit line-width requirements, and each set entry refers to a single rotomer/conformation. We will get back to these as well.

* "D" lines are "clusters", and are a feature of older db2 files. The cluster lines were meant to group together similar rotomers, though the idea was scrapped and newer db2 files no longer have them, thus "D" may now stand for "deprecated"
## What is the deal with sets and confs?

First off- the naming is super confusing. Confs don't describe conformations(rotomers), rather they describe a subset of a conformation. "Set" is an annoyingly general term for something that is so specific. These are the terms we are stuck with though! If I could turn back time, I would probably rename set->conformation and conf->subconformation

Sets are collections of confs, and confs are collections of atom xyz positions. Each set describes a conformation(rotomer) of the molecule. A single conf may be a part of one or more sets, thus confs can be used to define exactly where sets overlap. It is common to see confs comprised of a single atom- this may describe a situation where two rotomers happen to overlap at a particular atom. A conf may also be a part of just one set- this describes a situation where none of the atoms in the conf overlap with other sets, and is quite common.

## Example

Learning by example is best, so let's imagine the following (fake) rotomer set for a very simple molecule
```
+-------+
|   H   |
|   |   |
|O==C   |
+-------+
|       |
|       |
|O==C--H|
+-------+
```

Let us say the C atom lies at (0, 0). A somewhat abbreviated db2 for the set would look like this:
```
M fake
A 1 O
A 2 C
A 3 H
B 1 1 2 2
B 2 2 3 1
X 1 1 1 -1.000 +0.000
X 2 2 1 +0.000 +0.000
X 3 3 2 +0.000 +1.000
X 4 3 3 +1.000 +0.000
R 1 7 -1.000 +0.000
R 2 7 +0.000 +0.000
C 1 1 2
C 2 3 3
C 3 4 4
S 1 1 2 0 0 +0.000 +0.000
S 1 2 1 2
S 2 1 2 0 0 +0.000 +0.000
S 1 2 1 3
```
