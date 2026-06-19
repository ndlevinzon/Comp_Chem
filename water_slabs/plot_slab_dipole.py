import MDAnalysis as mda
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ============================================================
# USER SETTINGS
# ============================================================

topology = "../../build/opc_slab.prmtop"
trajectory = "warm_prod_centered_slab.nc"

z_min = -100.0
z_max = 100.0
dz = 0.5

water_resname = "WAT"
oxygen_name = "O"
hydrogen_selector = "name H*"

# ============================================================
# LOAD
# ============================================================

u = mda.Universe(topology, trajectory)

z_bins = np.arange(z_min, z_max + dz, dz)
z_centers = 0.5 * (z_bins[:-1] + z_bins[1:])
n_bins = len(z_centers)

# ============================================================
# BUILD WATER ATOM INDEX ARRAYS ONCE
# ============================================================

waters = u.select_atoms(f"resname {water_resname}").residues

O_indices = []
H1_indices = []
H2_indices = []

for wat in waters:
    O_atoms = wat.atoms.select_atoms(f"name {oxygen_name}")
    H_atoms = wat.atoms.select_atoms(hydrogen_selector)

    if len(O_atoms) == 1 and len(H_atoms) == 2:
        O_indices.append(O_atoms.indices[0])
        H1_indices.append(H_atoms.indices[0])
        H2_indices.append(H_atoms.indices[1])

O_indices = np.array(O_indices, dtype=int)
H1_indices = np.array(H1_indices, dtype=int)
H2_indices = np.array(H2_indices, dtype=int)

print(f"Waters analyzed: {len(O_indices)}")

# ============================================================
# ACCUMULATORS
# ============================================================

sum_cos = np.zeros(n_bins)
sum_p2 = np.zeros(n_bins)
counts = np.zeros(n_bins)

# ============================================================
# LOOP OVER FRAMES, VECTORIZED OVER WATERS
# ============================================================

for ts in u.trajectory:

    pos = u.atoms.positions

    O = pos[O_indices]
    H1 = pos[H1_indices]
    H2 = pos[H2_indices]

    H_mid = 0.5 * (H1 + H2)

    mu = H_mid - O
    mu_norm = np.linalg.norm(mu, axis=1)

    valid = mu_norm > 0.0

    z = O[:, 2]
    cos_theta = np.zeros_like(z)
    cos_theta[valid] = mu[valid, 2] / mu_norm[valid]

    p2 = 0.5 * (3.0 * cos_theta**2 - 1.0)

    # Bin by z
    bin_idx = np.digitize(z, z_bins) - 1
    valid_bins = (bin_idx >= 0) & (bin_idx < n_bins)

    bin_idx = bin_idx[valid_bins]
    cos_theta = cos_theta[valid_bins]
    p2 = p2[valid_bins]

    sum_cos += np.bincount(bin_idx, weights=cos_theta, minlength=n_bins)
    sum_p2 += np.bincount(bin_idx, weights=p2, minlength=n_bins)
    counts += np.bincount(bin_idx, minlength=n_bins)

# ============================================================
# FINAL AVERAGES
# ============================================================

mean_cos = np.divide(sum_cos, counts, out=np.full_like(sum_cos, np.nan), where=counts > 0)
mean_p2 = np.divide(sum_p2, counts, out=np.full_like(sum_p2, np.nan), where=counts > 0)

profile = pd.DataFrame({
    "z": z_centers,
    "mean_cos_theta": mean_cos,
    "mean_P2": mean_p2,
    "count": counts,
})

profile.to_csv("water_orientation_z_profile.csv", index=False)

# ============================================================
# PLOTS
# ============================================================

plt.figure(figsize=(6, 8))
plt.plot(profile["mean_cos_theta"], profile["z"], linewidth=2)
plt.axvline(0, color="black", linestyle="--", alpha=0.5)
plt.axhline(0, color="black", linestyle="--", alpha=0.5)
plt.xlabel(r"$\langle \cos\theta \rangle$")
plt.ylabel("Z coordinate (Å)")
plt.title("Water Dipole Orientation Across Slab")
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig("water_cos_theta_z.png", dpi=300)
plt.show()

plt.figure(figsize=(6, 8))
plt.plot(profile["mean_P2"], profile["z"], linewidth=2)
plt.axvline(0, color="black", linestyle="--", alpha=0.5)
plt.axhline(0, color="black", linestyle="--", alpha=0.5)
plt.xlabel(r"$\langle P_2(\cos\theta) \rangle$")
plt.ylabel("Z coordinate (Å)")
plt.title("Water Orientational Order Across Slab")
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig("water_P2_z.png", dpi=300)
plt.show()
