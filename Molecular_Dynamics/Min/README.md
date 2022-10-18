# AmberMdPrep.sh
Wrapper Script for Preparing Explicitly Solvated Systems for MD with AMBER

# Usage
```
AmberMdPrep.sh <input options>
Input File Options:
  Required:"
    -p <file>            : Topology
    -c <file>            : Coordinates
    --temp <temp>        : Temperature
  Optional:"
    -i <file>            : File containing custom minimization steps
    --thermo <type>      : Thermostat: berendsen, langevin (default)
    --baro <type>        : Barostat: berendsen, montecarlo (default)
    --finalthermo <type> : Final stage thermostat: berendsen, langevin (default)
    --finalbaro <type>   : Final stage barostat: berendsen, montecarlo (default)
    --nsolute <#>        : Number of solute residues
    --type <type>        : Residues type {protein* | nucleic}; determines backbone mask
    --mask <mask>        : Additional mask to use for restraints during steps 1-8
    --ares <name>        : Residue name to add to heavy atom masks if present
    --pmask <mask>       : Restraint mask to use during \"production\" (steps 9 and above)
    --pwt <weight>       : Restraint weight to use for '--pmask'; required if '--pmask' specified
    --pref <file>        : Optional reference structure to use for '--pmask'
    --charmmwater        : If specified assume CHARMM water (i.e. 'TIP3')
    --cutoff <cut>       : If specified, override default cutoffs with <cut>
    --test               : Test only. Do not run
    --norestart          : Do standard Eq with no restarts
    --skipfinaleq        : If specified, skip final eq. (step 10)
    --nprocs <#>         : Number of CPU processes to use (default 4)
    -O                   : Overwrite existing files, otherwise skip
    --keyhelp            : Print help for recognized input file keys
    --statusfile <file>  : Status file for final density eq
  Environment vars:
    PROG_MIN             : Command for minimization steps
    PROG_MD              : Command for MD steps
```
