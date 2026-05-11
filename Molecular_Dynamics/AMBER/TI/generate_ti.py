from pathlib import Path
import shutil
import textwrap

# ============================================================
# USER SETTINGS
# ============================================================

PROJECT_DIR = Path("C:/Users/ndlev/OneDrive/Documents/Research/mg/mg_ff_NEW/york_analysis/ti/grotz_m")

# Input files
START_RST7 = PROJECT_DIR / "step9.rst7"

CHARGE_TOP = PROJECT_DIR / "mg_opc_grotz_m_fix_box.prmtop"
# Mg has charge +2, LJ intact

VDW_TOP_V0 = PROJECT_DIR / "mg_opc_grotz_m_fix_box_neutral.prmtop"
VDW_TOP_V1 = PROJECT_DIR / "mg_opc_grotz_m_fix_box_dummy.prmtop"

TI_MASK = ":MG"

# Lambda schedules
CHARGE_LAMBDAS = [
    0.000, 0.050, 0.100, 0.150, 0.200,
    0.300, 0.400, 0.500, 0.600, 0.700,
    0.800, 0.850, 0.900, 0.950, 1.000
]

VDW_LAMBDAS = [
    0.000, 0.025, 0.050, 0.075, 0.100,
    0.150, 0.200, 0.300, 0.400, 0.500,
    0.600, 0.700, 0.800, 0.900,
    0.950, 0.975, 1.000
]

# Simulation parameters
TEMP = 300.0
CUT = 9.0
DT = 0.002
PROD_NS = 2.0
EQUIL_PS = 250.0

NTWX = 5000
NTPR = 500
NTWR = 5000

MAXCYC_MIN = 5000

TI_MASK = ":MG"   # change if your Mg residue/atom name differs

# Softcore parameters for vdW leg
SCALPHA = 0.5
SCBETA = 12.0

# Optional SLURM script generation    # e.g. "pmemd.cuda", "pmemd.MPI", "sander"
SLURM_PARTITION = "lonepeak-gpu"
SLURM_ACCOUNT = "lonepeak-gpu"             # e.g. "my_account" or None
SLURM_TIME = "04:00:00"
SLURM_MEM = "4G"


# ============================================================
# INTERNAL FUNCTIONS
# ============================================================

def write_text_file(path, text, executable=False):
    """
    Write text files safely from Windows for Linux/HPC use.
    Forces UTF-8 and Unix LF line endings.
    """
    text = text.replace("\r\n", "\n").replace("\r", "\n")
    path.write_text(text, encoding="utf-8", newline="\n")

    if executable:
        path.chmod(0o755)

def steps_from_time(time_value, units, dt):
    """
    Convert ps or ns to MD steps.
    dt is in ps.
    """
    if units == "ps":
        total_ps = time_value
    elif units == "ns":
        total_ps = time_value * 1000.0
    else:
        raise ValueError("units must be 'ps' or 'ns'")

    return int(round(total_ps / dt))


def format_lambda(lam):
    return f"{lam:.3f}"


def copy_if_exists(src, dst):
    """
    Copy file if it exists. If not, create a placeholder note.
    """
    dst.parent.mkdir(parents=True, exist_ok=True)

    if src.exists():
        shutil.copy2(src, dst)
    else:
        warning = (
            f"WARNING: Source file not found:\n"
            f"{src}\n\n"
            f"Expected destination:\n"
            f"{dst}\n\n"
            f"Place the correct file here before running AMBER.\n"
        )
        write_text_file(dst.with_suffix(dst.suffix + ".MISSING.txt"), warning)


def write_min_mdin(path, title="Minimization"):
    text = f"""\
{title}
&cntrl
  imin=1,
  maxcyc={MAXCYC_MIN},
  ncyc={MAXCYC_MIN // 2},

  ntb=1,
  cut={CUT},

  ntpr=100,
/
"""
    write_text_file(path, text)


def write_equil_mdin_pair(wdir, lam, leg):
    nstlim = steps_from_time(EQUIL_PS, "ps", DT)

    if leg == "charge":

        base = f"""\
Equilibration
&cntrl
  imin=0,
  irest=0,
  ntx=1,

  ntb=2,
  ntp=1,
  cut={CUT},

  ntt=3,
  gamma_ln=2.0,
  temp0={TEMP},
  tempi={TEMP},

  ntc=2,
  ntf=2,

  dt={DT},
  nstlim={nstlim},

  ntpr={NTPR},
  ntwx={NTWX},
  ntwr={NTWR},

  icfe=1,
  clambda={lam:.3f},
  ifsc=0,
"""

        v0 = base + "  crgmask='',\n/\n"
        v1 = base + f"  crgmask='{TI_MASK}',\n/\n"

    elif leg == "vdw":

        base = f"""\
Equilibration
&cntrl
  imin=0,
  irest=0,
  ntx=1,

  ntb=2,
  ntp=1,
  cut={CUT},

  ntt=3,
  gamma_ln=2.0,
  temp0={TEMP},
  tempi={TEMP},

  ntc=2,
  ntf=1,

  dt={DT},
  nstlim={nstlim},

  ntpr={NTPR},
  ntwx={NTWX},
  ntwr={NTWR},

  icfe=1,
  clambda={lam:.3f},
  ifsc=1,
  scalpha={SCALPHA},
  scbeta={SCBETA},
"""

        v0 = base + f"  scmask='{TI_MASK}',\n  logdvdl=1,\n/\n"
        v1 = base + f"  scmask='{TI_MASK}',\n  logdvdl=1,\n/\n"

    write_text_file(wdir / "equil_v0.in", v0)
    write_text_file(wdir / "equil_v1.in", v1)


def write_prod_mdin_pair(wdir, lam, leg):
    nstlim = steps_from_time(PROD_NS, "ns", DT)

    if leg == "charge":

        base = f"""\
Production
&cntrl
  imin=0,
  irest=1,
  ntx=5,

  cut={CUT},

  ntt=3,
  gamma_ln=2.0,
  temp0={TEMP},

  ntc=2,
  ntf=2,

  dt={DT},
  nstlim={nstlim},

  ntpr={NTPR},
  ntwx={NTWX},
  ntwr={NTWR},

  icfe=1,
  clambda={lam:.3f},
  ifsc=0,
"""

        v0 = base + "  crgmask='',\n/\n"
        v1 = base + f"  crgmask='{TI_MASK}',\n/\n"

    elif leg == "vdw":

        base = f"""\
Production
&cntrl
  imin=0,
  irest=1,
  ntx=5,

  cut={CUT},

  ntt=3,
  gamma_ln=2.0,
  temp0={TEMP},

  ntc=2,
  ntf=1,

  dt={DT},
  nstlim={nstlim},

  ntpr={NTPR},
  ntwx={NTWX},
  ntwr={NTWR},

  icfe=1,
  clambda={lam:.3f},
  ifsc=1,
  scalpha={SCALPHA},
  scbeta={SCBETA},
"""

        v0 = base + f"  scmask='{TI_MASK}',\n  logdvdl=1,\n/\n"
        v1 = base + f"  scmask='{TI_MASK}',\n  logdvdl=1,\n/\n"

    write_text_file(wdir / "prod_v0.in", v0)
    write_text_file(wdir / "prod_v1.in", v1)

def write_groupfiles(wdir, leg):

    if leg == "charge":
        p0 = "prmtop"
        p1 = "prmtop"

        write_text_file(wdir / "group_min", f"""\
-O -i min.in -o min_v0.out -p {p0} -c start.rst7 -r min_v0.rst7 -inf min_v0.info
-O -i min.in -o min_v1.out -p {p1} -c start.rst7 -r min_v1.rst7 -inf min_v1.info
""")

        write_text_file(wdir / "group_equi", f"""\
-O -i equil_v0.in -o equil_v0.out -p {p0} -c min_v0.rst7 -r equil_v0.rst7 -x equil_v0.nc -inf equil_v0.info
-O -i equil_v1.in -o equil_v1.out -p {p1} -c min_v1.rst7 -r equil_v1.rst7 -x equil_v1.nc -inf equil_v1.info
""")

    elif leg == "vdw":
        p0 = "prmtop_v0"
        p1 = "prmtop_v1"

        # No separate minimization-derived coordinates for v0/v1.
        # Both groups start from identical start.rst7.
        write_text_file(wdir / "group_min", f"""\
-O -i min.in -o min_v0.out -p {p0} -c start.rst7 -r min_v0.rst7 -inf min_v0.info
-O -i min.in -o min_v1.out -p {p1} -c start.rst7 -r min_v1.rst7 -inf min_v1.info
""")

        write_text_file(wdir / "group_equi", f"""\
-O -i equil_v0.in -o equil_v0.out -p {p0} -c start.rst7 -r equil_v0.rst7 -x equil_v0.nc -inf equil_v0.info
-O -i equil_v1.in -o equil_v1.out -p {p1} -c start.rst7 -r equil_v1.rst7 -x equil_v1.nc -inf equil_v1.info
""")

    else:
        raise ValueError("leg must be charge or vdw")

    if leg == "charge":
        prod_c0 = "equil_v0.rst7"
        prod_c1 = "equil_v1.rst7"
        p0 = "prmtop"
        p1 = "prmtop"
    else:
        prod_c0 = "equil_v0.rst7"
        prod_c1 = "equil_v1.rst7"
        p0 = "prmtop_v0"
        p1 = "prmtop_v1"

    write_text_file(wdir / "group_prod", f"""\
-O -i prod_v0.in -o prod_v0.out -p {p0} -c {prod_c0} -r prod_v0.rst7 -x prod_v0.nc -inf prod_v0.info
-O -i prod_v1.in -o prod_v1.out -p {p1} -c {prod_c1} -r prod_v1.rst7 -x prod_v1.nc -inf prod_v1.info
""")


def write_run_script(path):
    text = """\
#!/bin/bash
set -e
set -o pipefail

export OMPI_MCA_pml=ob1
export OMPI_MCA_btl=^openib
export AMBERHOME="/uufs/chpc.utah.edu/sys/installdir/r8/amber/24"
ENGINE="${AMBERHOME}/bin/sander.MPI"
MPI_NP="${MPI_NP:-96}"

echo "Running in $(pwd)"
echo "ENGINE = ${ENGINE}"
echo "MPI_NP = ${MPI_NP}"

mpirun -np ${MPI_NP} "${ENGINE}" -ng 2 -groupfile group_min
mpirun -np ${MPI_NP} "${ENGINE}" -ng 2 -groupfile group_equi
mpirun -np ${MPI_NP} "${ENGINE}" -ng 2 -groupfile group_prod
"""
    write_text_file(path, text, executable=True)


def make_leg(leg_name, lambdas, topology_v0, topology_v1=None):
    leg_dir = PROJECT_DIR / leg_name
    leg_dir.mkdir(parents=True, exist_ok=True)

    for lam in lambdas:
        lam_str = format_lambda(lam)
        wdir = leg_dir / f"lambda_{lam_str}"
        wdir.mkdir(parents=True, exist_ok=True)

        copy_if_exists(START_RST7, wdir / "start.rst7")

        leg_type = "charge" if "charge" in leg_name else "vdw"

        if leg_type == "charge":
            copy_if_exists(topology_v0, wdir / "prmtop")
        else:
            copy_if_exists(topology_v0, wdir / "prmtop_v0")
            copy_if_exists(topology_v1, wdir / "prmtop_v1")

        write_min_mdin(
            wdir / "min.in",
            title=f"Minimization for {leg_type} leg, lambda={lam_str}"
        )
        write_equil_mdin_pair(wdir, lam, leg_type)
        write_prod_mdin_pair(wdir, lam, leg_type)
        write_groupfiles(wdir, leg_type)
        write_run_script(wdir / "run_window.sh")


def write_manifest():
    manifest = f"""\
Mg2+ OPC TI solvation setup

Project directory:
  {PROJECT_DIR}

Input restart:
  {START_RST7}

Charge leg:
  A topology: {CHARGE_TOP}
  lambdas: {CHARGE_LAMBDAS}

vdW leg:
  V0 topology (neutral LJ-on): {VDW_TOP_V0}
  V1 topology (dummy LJ-off): {VDW_TOP_V1}
  lambdas: {VDW_LAMBDAS}

Simulation parameters:
  TEMP = {TEMP}
  CUT = {CUT}
  DT = {DT} ps
  EQUIL_PS = {EQUIL_PS}
  PROD_NS = {PROD_NS}
  SCALPHA = {SCALPHA}
  SCBETA = {SCBETA}

Important:
  This script creates templates and directory structure.
  You must confirm AMBER TI topology/mask syntax for your AMBER engine/version.
  Mg2+ charge-removal free energies require finite-size/convention corrections.
"""
    (PROJECT_DIR / "README_TI_SETUP.txt").write_text(manifest)


# ============================================================
# MAIN SETUP
# ============================================================

PROJECT_DIR.mkdir(parents=True, exist_ok=True)

make_leg("01_charge_off", CHARGE_LAMBDAS, CHARGE_TOP)

make_leg("02_vdw_off", VDW_LAMBDAS, VDW_TOP_V0, VDW_TOP_V1)

write_manifest()

print(f"Created TI setup in: {PROJECT_DIR.resolve()}")
print("Review README_TI_SETUP.txt before running.")
