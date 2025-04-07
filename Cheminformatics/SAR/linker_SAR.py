import pandas as pd
import numpy as np
from protlearn.features import aaindex1
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import f_regression
from sklearn.metrics import r2_score
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

    # Extract single AA per position
    position_aa = [seq[pos] for seq in linkers]

    # Encode AA at this position
    X_encoded, descriptors = aaindex1(position_aa)
    descriptors = list(descriptors)
    position_encodings.append(X_encoded)  # Store for pairwise later

    # Train-test split
    X_train, X_test, y_train, y_test = train_test_split(X_encoded, y, test_size=0.2, random_state=42)

    # Standardize
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)

    # Fit model
    model = RandomForestRegressor(random_state=42)
    model.fit(X_train_scaled, y_train)
    y_pred = model.predict(X_test_scaled)
    r2 = r2_score(y_test, y_pred)
    print(f"R² on test set ({label}): {r2:.2f}")

    # Feature importance and statistics
    importances = model.feature_importances_
    F, p_vals = f_regression(X_train_scaled, y_train)

    results = pd.DataFrame({
        "Feature": descriptors,
        "Importance": importances,
        "F_score": F,
        "p_value": p_vals
    }).sort_values("Importance", ascending=False)

    print(results.head(10))

    # --- Plot: Feature Importance Barplot ---
    plt.figure(figsize=(10, 6))
    sns.barplot(x="Importance", y="Feature", data=results.head(10))
    plt.title(f"Top 10 Features at {label} (R² = {r2:.2f})")
    plt.xlabel("Feature Importance")
    plt.ylabel("AAIndex1 Descriptor")
    plt.tight_layout()
    plt.show()

    # --- Plot: SAR Line Plot for Top Significant Features (p < 0.05) ---
    sig_results = results[results['p_value'] < 0.05].head(3)
    for feat in sig_results['Feature']:
        feat_index = descriptors.index(feat)
        plt.figure(figsize=(6, 4))
        sns.regplot(x=X_encoded[:, feat_index], y=y)
        plt.title(f"SAR: {feat} vs. Fluorescence ({label})")
        plt.xlabel(feat)
        plt.ylabel("Fluorescence")
        plt.tight_layout()
        plt.show()

# --- Combined Interaction Analysis: AA70*AA71 and AA313*AA314 ---
interaction_pairs = [(0, 1), (2, 3)]
for i, j in interaction_pairs:
    label_i = position_labels[i]
    label_j = position_labels[j]
    print(f"\n--- Analyzing Interaction: {label_i} × {label_j} ---")

    Xi = position_encodings[i]
    Xj = position_encodings[j]

    # Element-wise multiplication of descriptor vectors
    X_combined = Xi * Xj
    combo_labels = [f"{a}*{b}" for a, b in zip(descriptors, descriptors)]

    # Train-test split
    X_train, X_test, y_train, y_test = train_test_split(X_combined, y, test_size=0.2, random_state=42)

    # Standardize
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)

    # Fit model
    model = RandomForestRegressor(random_state=42)
    model.fit(X_train_scaled, y_train)
    y_pred = model.predict(X_test_scaled)
    r2 = r2_score(y_test, y_pred)
    print(f"R² on test set ({label_i}×{label_j}): {r2:.2f}")

    # Feature importance
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

    # SAR line plots for top 3 significant combos
    sig_results = results[results['p_value'] < 0.05].head(3)
    for feat in sig_results['Feature']:
        feat_index = combo_labels.index(feat)
        plt.figure(figsize=(6, 4))
        sns.regplot(x=X_combined[:, feat_index], y=y)
        plt.title(f"SAR: {feat} vs. Fluorescence ({label_i}×{label_j})")
        plt.xlabel(feat)
        plt.ylabel("Fluorescence")
        plt.tight_layout()
        plt.show()
