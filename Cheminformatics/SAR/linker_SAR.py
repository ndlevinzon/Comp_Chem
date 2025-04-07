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
linkers = df['linker'].tolist()
y = df['fluor_%Dif'].values

# Sanity check: all linkers should be 4 AAs
assert all(len(seq) == 4 for seq in linkers), "Non-4mer linker found!"

# Position labels
position_labels = ["AA70", "AA71", "AA313", "AA314"]
position_encodings = []

# Position-wise encoding and analysis
for pos in range(4):
    label = position_labels[pos]
    print(f"\n--- Analyzing {label} ---")

    # Extract and encode
    position_aa = [seq[pos] for seq in linkers]
    X_encoded, descriptors = aaindex1(position_aa)
    descriptors = list(descriptors)
    position_encodings.append(X_encoded)

    # Train-test split
    X_train, X_test, y_train, y_test = train_test_split(X_encoded, y, test_size=0.2, random_state=42)

    # Standardize
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)

    # Random Forest
    model = RandomForestRegressor(random_state=42)
    model.fit(X_train_scaled, y_train)
    y_pred = model.predict(X_test_scaled)
    r2 = r2_score(y_test, y_pred)
    print(f"R² on test set ({label}): {r2:.2f}")

    # Feature importance + F-test
    importances = model.feature_importances_
    F, p_vals = f_regression(X_train_scaled, y_train)
    results = pd.DataFrame({
        "Feature": descriptors,
        "Importance": importances,
        "F_score": F,
        "p_value": p_vals
    }).sort_values("Importance", ascending=False)

    print(results.head(10))

    # Plot top 10 importance
    plt.figure(figsize=(10, 6))
    sns.barplot(x="Importance", y="Feature", data=results.head(10))
    plt.title(f"Top 10 Features at {label} (R² = {r2:.2f})")
    plt.xlabel("Feature Importance")
    plt.ylabel("AAIndex1 Descriptor")
    plt.tight_layout()
    plt.show()

    # --- Vectorized Linear R² ---
    X_centered = X_encoded - X_encoded.mean(axis=0)
    y_centered = y - y.mean()
    Xy = np.dot(X_centered.T, y_centered)
    X2 = np.sum(X_centered ** 2, axis=0)
    y2 = np.sum(y_centered ** 2)
    r_squared = (Xy ** 2) / (X2 * y2)
    r_squared = np.nan_to_num(r_squared)

    # --- Plot Top 100 SARs in Grid ---
    top_indices = np.argsort(np.abs(r_squared))[-100:][::-1]
    fig, axes = plt.subplots(10, 10, figsize=(25, 25))
    fig.suptitle(f"Top 100 SAR Features by |R²| ({label})", fontsize=20)

    for idx, ax in zip(top_indices, axes.flat):
        sns.regplot(x=X_encoded[:, idx], y=y, ax=ax, scatter_kws={'s': 10, 'alpha': 0.7}, line_kws={'color': 'red'})
        ax.set_title(f"{descriptors[idx]}\nR²={r_squared[idx]:.2f}", fontsize=8)
        ax.set_xticks([])
        ax.set_yticks([])


    plt.tight_layout(rect=[0, 0, 1, 0.97])
    plt.savefig(f"C:/Users/ndlev/OneDrive/Documents/Research/diehl/pancake_test/biosensor_ml/results/SAR_2/{label}_top100_SAR_grid.png", dpi=300)
    plt.close()

    # --- Reconcile with Pearson correlation ---
    r_vals = []
    p_vals = []

    X_target = X_encoded if 'X_encoded' in locals() and X_encoded.shape[0] == len(y) else X_combined

    for i in range(X_target.shape[1]):
        r, p = pearsonr(X_target[:, i], y)
        r_vals.append(r)
        p_vals.append(p)

    r_vals = np.nan_to_num(r_vals)
    r_squared_pearson = np.array(r_vals) ** 2

    # Add to results dataframe
    results["Pearson_r"] = r_vals
    results["R_squared"] = r_squared_pearson
    results["Pearson_pval"] = p_vals
    results["abs_r"] = np.abs(r_vals)
    results["combined_score"] = results["Importance"] * results["abs_r"]

    plt.figure(figsize=(6, 6))
    sns.scatterplot(data=results, x="abs_r", y="Importance")
    plt.title(f"Feature Importance vs |Pearson r| ({label})")
    plt.xlabel("|Pearson Correlation|")
    plt.ylabel("Random Forest Importance")
    plt.tight_layout()
    plt.show()


# --- Interaction Analysis ---
interaction_pairs = [(0, 1), (2, 3)]
for i, j in interaction_pairs:
    label_i = position_labels[i]
    label_j = position_labels[j]
    print(f"\n--- Analyzing Interaction: {label_i} × {label_j} ---")

    Xi = position_encodings[i]
    Xj = position_encodings[j]
    n_samples, n_features = Xi.shape

    # Efficient einsum-based cross feature generation
    X_combined = np.einsum('ij,ik->ijk', Xi, Xj).reshape(n_samples, -1)
    combo_labels = [f"{a}*{b}" for a in descriptors for b in descriptors]

    # Train-test split
    X_train, X_test, y_train, y_test = train_test_split(X_combined, y, test_size=0.2, random_state=42)

    # Standardize
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)

    # Random Forest
    model = RandomForestRegressor(random_state=42, n_estimators=100, max_features="sqrt", n_jobs=-1)
    model.fit(X_train_scaled, y_train)
    y_pred = model.predict(X_test_scaled)
    r2 = r2_score(y_test, y_pred)
    print(f"R² on test set ({label_i}×{label_j}): {r2:.2f}")

    # Feature importance + F-test
    importances = model.feature_importances_
    F, p_vals = f_regression(X_train_scaled, y_train)
    results = pd.DataFrame({
        "Feature": combo_labels,
        "Importance": importances,
        "F_score": F,
        "p_value": p_vals
    }).sort_values("Importance", ascending=False)

    print(results.head(10))

    # Plot top 10 combo features
    plt.figure(figsize=(10, 6))
    sns.barplot(x="Importance", y="Feature", data=results.head(10))
    plt.title(f"Top 10 Interaction Features: {label_i} × {label_j} (R² = {r2:.2f})")
    plt.xlabel("Interaction Importance")
    plt.ylabel("Feature Pair")
    plt.tight_layout()
    plt.show()

    # --- Vectorized Linear R² for interactions ---
    X_centered = X_combined - X_combined.mean(axis=0)
    y_centered = y - y.mean()
    Xy = np.dot(X_centered.T, y_centered)
    X2 = np.sum(X_centered ** 2, axis=0)
    y2 = np.sum(y_centered ** 2)
    r_squared = (Xy ** 2) / (X2 * y2)
    r_squared = np.nan_to_num(r_squared)

    # --- Plot Top 100 SARs in Grid ---
    top_indices = np.argsort(np.abs(r_squared))[-100:][::-1]
    fig, axes = plt.subplots(10, 10, figsize=(25, 25))
    fig.suptitle(f"Top 100 Interaction SARs by |R²| ({label_i}×{label_j})", fontsize=20)

    for idx, ax in zip(top_indices, axes.flat):
        sns.regplot(x=X_combined[:, idx], y=y, ax=ax, scatter_kws={'s': 10, 'alpha': 0.7}, line_kws={'color': 'red'})
        ax.set_title(f"{combo_labels[idx]}\nR²={r_squared[idx]:.2f}", fontsize=6)
        ax.set_xticks([])
        ax.set_yticks([])

    plt.tight_layout(rect=[0, 0, 1, 0.97])
    plt.savefig(f"C:/Users/ndlev/OneDrive/Documents/Research/diehl/pancake_test/biosensor_ml/results/SAR_2/{label_i}_{label_j}_top100_SAR_grid.png", dpi=300)
    plt.close()

    # --- Reconcile with Pearson correlation ---
    r_vals = []
    p_vals = []

    # Determine correct matrix and label set
    if 'X_combined' in locals() and X_combined.shape[1] == len(results):
        X_target = X_combined
        labels = combo_labels
    else:
        X_target = X_encoded
        labels = descriptors

    # Ensure feature count matches
    assert X_target.shape[1] == len(results), f"Mismatch: {X_target.shape[1]} features vs {len(results)} result rows"

    # Compute Pearson stats
    for i in range(X_target.shape[1]):
        r, p = pearsonr(X_target[:, i], y)
        r_vals.append(r)
        p_vals.append(p)

    r_vals = np.nan_to_num(r_vals)
    r_squared_pearson = np.array(r_vals) ** 2

    # Add to results dataframe
    results["Pearson_r"] = r_vals
    results["R_squared"] = r_squared_pearson
    results["Pearson_pval"] = p_vals
    results["abs_r"] = np.abs(r_vals)
    results["combined_score"] = results["Importance"] * results["abs_r"]

    plt.figure(figsize=(6, 6))
    sns.scatterplot(data=results, x="abs_r", y="Importance")
    plt.title(f"Feature Importance vs |Pearson r| ({label})")
    plt.xlabel("|Pearson Correlation|")
    plt.ylabel("Random Forest Importance")
    plt.tight_layout()
    plt.show()
