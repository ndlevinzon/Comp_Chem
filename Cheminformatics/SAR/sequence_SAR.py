# Updated SAR script with interaction terms between flanking windows

import pandas as pd
import numpy as np
from protlearn.features import aaindex1
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import f_regression
from sklearn.metrics import r2_score
from scipy.stats import pearsonr
import matplotlib.pyplot as plt
import seaborn as sns

# Load data
df = pd.read_csv("C:/Users/ndlev/OneDrive/Documents/Research/diehl/pancake_test/biosensor_ml/data/experimental_data.csv")
sequences = df['sequence'].tolist()
y = df['fluor_%Dif'].values

# Define linker positions and window size
linker_positions = [70, 71, 313, 314]
window_size = 7

# Helper to extract window with padding
def extract_window(seq, center, window):
    start = center - window
    end = center + window + 1
    padded = ("X" * max(0, -start)) + seq[max(0, start):min(len(seq), end)] + ("X" * max(0, end - len(seq)))
    return padded[:2 * window + 1]

# Extract and encode windows for each position
all_windows = {}
for pos in linker_positions:
    windows = [extract_window(seq, pos, window_size) for seq in sequences]
    encoded, desc = aaindex1(windows)
    all_windows[pos] = (encoded, list(desc))

# Flatten features across positions
X_combined = np.concatenate([all_windows[pos][0] for pos in linker_positions], axis=1)
features_combined = []
for pos in linker_positions:
    features_combined += [f"{pos}_{d}" for d in all_windows[pos][1]]

# --- Interaction terms between positions ---
interactions = [(70, 71), (313, 314)]
for a, b in interactions:
    X_a = all_windows[a][0]
    X_b = all_windows[b][0]
    desc_a = all_windows[a][1]
    desc_b = all_windows[b][1]

    X_inter = np.einsum('ij,ik->ijk', X_a, X_b).reshape(X_a.shape[0], -1)
    X_combined = np.concatenate([X_combined, X_inter], axis=1)
    features_combined += [f"{a}_{da}*{b}_{db}" for da in desc_a for db in desc_b]

# Train-test split
X_train, X_test, y_train, y_test = train_test_split(X_combined, y, test_size=0.2, random_state=42)

# Standardize
scaler = StandardScaler()
X_train_scaled = scaler.fit_transform(X_train)
X_test_scaled = scaler.transform(X_test)

# Train model
model = RandomForestRegressor(random_state=42, n_estimators=100, max_features='sqrt', n_jobs=-1)
model.fit(X_train_scaled, y_train)
y_pred = model.predict(X_test_scaled)
r2 = r2_score(y_test, y_pred)
print(f"Overall R² score: {r2:.3f}")

# Feature importance + Pearson r
importances = model.feature_importances_
r_vals = [pearsonr(X_combined[:, i], y)[0] for i in range(X_combined.shape[1])]
r_squared = np.array(r_vals) ** 2

results = pd.DataFrame({
    "Feature": features_combined,
    "Importance": importances,
    "abs_r": np.abs(r_vals),
    "R_squared": r_squared,
    "Combined": np.abs(r_vals) * importances
}).sort_values("Combined", ascending=False)

results.to_csv("C:/Users/ndlev/OneDrive/Documents/Research/diehl/pancake_test/biosensor_ml/results/SAR_3/linker_flank_interactions_results.csv", index=False)

# Histogram
plt.figure(figsize=(10, 6))
bin_data = pd.cut(results["abs_r"], bins=np.linspace(0, 1, 21))
sns.barplot(x=bin_data, y=results["Importance"], estimator=np.mean, ci=None)
plt.xticks(rotation=90)
plt.title("Importance vs |Pearson r| with Windowed Interactions")
plt.tight_layout()
plt.savefig("C:/Users/ndlev/OneDrive/Documents/Research/diehl/pancake_test/biosensor_ml/results/SAR_3/linker_flank_interactions_hist.png", dpi=300)
plt.close()

# Top 100 SAR grid
top_idx = np.argsort(results["Combined"])[-100:][::-1]
fig, axes = plt.subplots(10, 10, figsize=(25, 25))
fig.suptitle("Top 100 SARs (Flanking Linkers + Interactions)", fontsize=20)
for idx, ax in zip(top_idx, axes.flat):
    sns.regplot(x=X_combined[:, idx], y=y, ax=ax, scatter_kws={'s': 10, 'alpha': 0.7}, line_kws={'color': 'red'})
    ax.set_title(f"{results['Feature'].iloc[idx]}\nR²={results['R_squared'].iloc[idx]:.2f}", fontsize=7)
    ax.set_xticks([])
    ax.set_yticks([])
plt.tight_layout(rect=[0, 0, 1, 0.97])
plt.savefig("C:/Users/ndlev/OneDrive/Documents/Research/diehl/pancake_test/biosensor_ml/results/SAR_3/linker_flank_interactions_top100_grid.png", dpi=300)
plt.close()
