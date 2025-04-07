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

# Position-wise encoding and analysis
for pos in range(4):
    print(f"\n--- Analyzing Position {pos + 1} ---")

    # Extract single AA per position
    position_aa = [seq[pos] for seq in linkers]

    # Encode AA at this position
    X_encoded, descriptors = aaindex1(position_aa)
    descriptors = list(descriptors)

    feature_names = descriptors  # each feature is already named

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
    print(f"R² on test set (position {pos+1}): {r2:.2f}")

    # Feature importance and statistics
    importances = model.feature_importances_
    F, p_vals = f_regression(X_train_scaled, y_train)

    results = pd.DataFrame({
        "Feature": feature_names,
        "Importance": importances,
        "F_score": F,
        "p_value": p_vals
    }).sort_values("Importance", ascending=False)

    print(results.head(10))

    # --- Plot: Feature Importance Barplot ---
    plt.figure(figsize=(10, 6))
    sns.barplot(x="Importance", y="Feature", data=results.head(10))
    plt.title(f"Top 10 Features at Position {pos+1} (R² = {r2:.2f})")
    plt.xlabel("Feature Importance")
    plt.ylabel("AAIndex1 Descriptor")
    plt.tight_layout()
    plt.show()

    # --- Plot: SAR Line Plot for Top 3 Features ---
    top_feats = results.head(3)['Feature'].values
    for feat in top_feats:
        feat_index = descriptors.index(feat)
        plt.figure(figsize=(6, 4))
        sns.regplot(x=X_encoded[:, feat_index], y=y)
        plt.title(f"SAR: {feat} vs. Fluorescence (Position {pos+1})")
        plt.xlabel(feat)
        plt.ylabel("Fluorescence")
        plt.tight_layout()
        plt.show()
