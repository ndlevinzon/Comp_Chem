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

# Load your data
df = pd.read_csv("C:/Users/ndlev/OneDrive/Documents/Research/diehl/pancake_test/biosensor_ml/data/experimental_data.csv")
linkers = df['linker'].tolist()
y = df['fluor_%Dif'].values

# Encode using aaindex1 (physicochemical properties)
X_encoded, descriptors = aaindex1(linkers)
print("Number of sequences:", len(linkers))
print("Encoded shape:", X_encoded.shape)

# Skip reshaping — just use X_encoded as flat feature vector
X = X_encoded
feature_names = [f"feat_{i}" for i in range(X.shape[1])]

# Train-test split
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# Standardize
scaler = StandardScaler()
X_train_scaled = scaler.fit_transform(X_train)
X_test_scaled = scaler.transform(X_test)

# Fit regression model
model = RandomForestRegressor(random_state=42)
model.fit(X_train_scaled, y_train)
y_pred = model.predict(X_test_scaled)
r2 = r2_score(y_test, y_pred)
print(f"R² on test set: {r2:.2f}")

# Feature importance
importances = model.feature_importances_
F, p_vals = f_regression(X_train_scaled, y_train)

# Combine into dataframe
results = pd.DataFrame({
    "Feature": feature_names,
    "Importance": importances,
    "F_score": F,
    "p_value": p_vals
}).sort_values("Importance", ascending=False)

print(results.head(10))

# Plot top 10 features
plt.figure(figsize=(10, 6))
sns.barplot(x="Importance", y="Feature", data=results.head(10))
plt.title(f"Top 10 Physicochemical Features (R² = {r2:.2f})")
plt.xlabel("Feature Importance")
plt.ylabel("AAIndex1 Feature Index")
plt.tight_layout()
plt.show()
