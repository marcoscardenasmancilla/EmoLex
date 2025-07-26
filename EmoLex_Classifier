import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split, cross_val_score
from sklearn.preprocessing import LabelEncoder
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report, confusion_matrix
import joblib
import matplotlib.pyplot as plt
import seaborn as sns
import os
import datetime

# Load data
df = pd.read_csv('long_format_sub-clustering.csv')

# Drop rows with missing values
df = df.dropna(subset=['subgroup', 'Valence_Mean', 'Arousal_Mean', 'Concreteness_Mean', 'Emotionality', 'Zipf_EsPal', 'Balanced_Integration_Score'])

# Define features and target
X = df[['Valence_Mean', 'Arousal_Mean', 'Concreteness_Mean', 'Emotionality', 'Zipf_EsPal', 'Balanced_Integration_Score']]
y = df['subgroup']

# Encode target labels if they are not numeric
if y.dtype == 'O':
    le = LabelEncoder()
    y = le.fit_transform(y)
    label_mapping = dict(zip(le.classes_, le.transform(le.classes_)))
    print("Label mapping:", label_mapping)

# Split data into training and testing sets (80/20)
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# Train Random Forest Classifier
rf = RandomForestClassifier(n_estimators=100, random_state=42)
rf.fit(X_train, y_train)

# Predict on test data
y_pred = rf.predict(X_test)

# Evaluate model
print("Classification Report:\n", classification_report(y_test, y_pred))
print("Confusion Matrix:\n", confusion_matrix(y_test, y_pred))

# Save model
timestamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
model_filename = f"random_forest_model_{timestamp}.joblib"
joblib.dump(rf, model_filename)
print(f"Model saved to {model_filename}")

# Feature importance plot
importances = rf.feature_importances_
features = X.columns
indices = np.argsort(importances)[::-1]

plt.figure(figsize=(10, 6))
sns.barplot(x=importances[indices], y=features[indices])
plt.title("Feature Importances - Random Forest")
plt.tight_layout()
plt.show()
