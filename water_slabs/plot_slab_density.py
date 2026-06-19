import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import savgol_filter

# ============================================================
# USER SETTINGS
# ============================================================

density_file = "cool_water_O_number_density_z.dat"

Lx = 194.567
Ly = 188.781
dz = 0.25

# If cpptraj output is accumulated over frames, set this correctly.
# If already averaged, leave as 1.
n_frames = 1010

MW_WATER = 18.01528
NA = 6.02214076e23

# ============================================================
# LOAD CPPTRAJ OXYGEN NUMBER OUTPUT
# ============================================================

df = pd.read_csv(
    density_file,
    comment="#",
    sep=r"\s+",
    names=["z", "O_count_per_bin", "sd_O_count_per_bin"]
)

# ============================================================
# CONVERT O COUNTS TO WATER DENSITY
# ============================================================

bin_volume_A3 = Lx * Ly * dz

area_A2 = Lx * Ly

df["number_density_A3"] = df["O_count_per_bin"] / area_A2

df["density_g_cm3"] = (
    df["number_density_A3"] * MW_WATER / NA * 1.0e24
)

df["sd_number_density_A3"] = (
    df["sd_O_count_per_bin"] / n_frames / bin_volume_A3
)

df["sd_g_cm3"] = (
    df["sd_number_density_A3"] * MW_WATER / NA * 1.0e24
)

df.to_csv("warm_water_O_density_z_converted_g_cm3.csv", index=False)

# ============================================================
# AUTOMATIC PLATEAU / INTERFACE DETECTION
# ============================================================

z = df["z"].to_numpy()
rho = df["density_g_cm3"].to_numpy()

window = 51
polyorder = 2

if window >= len(rho):
    window = len(rho) - 1 if len(rho) % 2 == 0 else len(rho)

if window % 2 == 0:
    window -= 1

rho_smooth = savgol_filter(rho, window_length=window, polyorder=polyorder)

# Prevent smoothing from producing negative density near vapor/interface regions
rho_smooth = np.clip(rho_smooth, 0.0, None)

bulk_density = np.percentile(rho_smooth, 95)
plateau_threshold = 0.90 * bulk_density

bulk_density = np.percentile(rho_smooth, 95)
plateau_threshold = 0.90 * bulk_density

plateau_mask = rho_smooth >= plateau_threshold
indices = np.where(plateau_mask)[0]

if len(indices) == 0:
    raise RuntimeError("No plateau detected. Try lowering plateau_threshold.")

segments = np.split(indices, np.where(np.diff(indices) > 1)[0] + 1)
largest_segment = max(segments, key=len)

z_lower_bulk = z[largest_segment[0]]
z_upper_bulk = z[largest_segment[-1]]

vapor_threshold = 0.05 * bulk_density

water_mask = rho_smooth >= vapor_threshold
water_indices = np.where(water_mask)[0]

if len(water_indices) == 0:
    raise RuntimeError("No water slab detected. Try lowering vapor_threshold.")

z_lower_water = z[water_indices[0]]
z_upper_water = z[water_indices[-1]]

lower_interface = (z_lower_water, z_lower_bulk)
upper_interface = (z_upper_bulk, z_upper_water)
bulk_region = (z_lower_bulk, z_upper_bulk)

df["density_smooth_g_cm3"] = rho_smooth

print("=" * 70)
print("Automatically detected density regions")
print("=" * 70)
print(f"Bulk density estimate:       {bulk_density:.4f} g/cm^3")
print(f"Plateau threshold:           {plateau_threshold:.4f} g/cm^3")
print(f"Vapor/water threshold:       {vapor_threshold:.4f} g/cm^3")
print()
print("Bulk-like water region:")
print(f"  z = {bulk_region[0]:.3f} to {bulk_region[1]:.3f} Å")
print()
print("Lower interfacial water region:")
print(f"  z = {lower_interface[0]:.3f} to {lower_interface[1]:.3f} Å")
print()
print("Upper interfacial water region:")
print(f"  z = {upper_interface[0]:.3f} to {upper_interface[1]:.3f} Å")
print("=" * 70)

# ============================================================
# PLOT
# ============================================================

fig, ax = plt.subplots(figsize=(6, 8))

ax.plot(df["density_g_cm3"], df["z"], linewidth=2, label="Density")
ax.plot(df["density_smooth_g_cm3"], df["z"], linewidth=2, linestyle="--", label="Smoothed density")

ax.fill_betweenx(
    df["z"],
    df["density_g_cm3"] - df["sd_g_cm3"],
    df["density_g_cm3"] + df["sd_g_cm3"],
    alpha=0.3,
    label="±1 SD"
)

ax.axhline(0, color="black", linestyle="--", alpha=0.5)

ax.axhline(z_lower_bulk, color="red", linestyle="--", linewidth=2, label="Bulk/interfacial boundary")
ax.axhline(z_upper_bulk, color="red", linestyle="--", linewidth=2)

ax.axhspan(z_lower_bulk, z_upper_bulk, alpha=0.12, color="red", label="Bulk-like region")
ax.axhspan(lower_interface[0], lower_interface[1], alpha=0.12, color="orange", label="Interfacial region")
ax.axhspan(upper_interface[0], upper_interface[1], alpha=0.12, color="orange")

ax.set_xlabel(r"Oxygen Atom Density (g/cm$^3$)")
ax.set_ylabel("Z Coordinate (Å)")
ax.set_ylim(-100, 100)
ax.set_title("TIP3P Water Density Profile Along Z @ 78K")
ax.grid(alpha=0.3)
ax.legend()

plt.tight_layout()
plt.savefig("water_density_profile_g_cm3.png", dpi=300)
plt.show()
