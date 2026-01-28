import pandas as pd
import matplotlib
#matplotlib.use("Agg")   # Use non-GUI backend

import matplotlib.pyplot as plt
from sklearn.metrics import mean_squared_error
from scipy import stats

df = pd.read_csv("IDO1_exp_inhibition.csv")

col_x = "inhibition_rate_efo21_1"
col_y = "inhibition_rate_sk0v3_1"

data = df[[col_x, col_y]].dropna()
x = data[col_x]
y = data[col_y]

# -------------------------------
# 1. Scatter plot + correlation line
# -------------------------------
slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
line = slope * x + intercept

plt.figure(figsize=(6, 6))
plt.scatter(x, y, alpha=0.7, label="Data points")
plt.plot(x, line, color="red", label=f"Fit line (R²={r_value**2:.2f})")
plt.xlabel(col_x)
plt.ylabel(col_y)
plt.legend()
plt.title("Correlation between experiments")
plt.grid(True)
plt.show()
plt.savefig("scatter_plot.png", dpi=300)

# -------------------------------
# 2. Calculate metrics
# -------------------------------
mse = mean_squared_error(y, line)
print(f"Correlation coefficient (r): {r_value:.3f}")
print(f"P-value: {p_value:.3e}")
print(f"Mean Squared Error (MSE): {mse:.4f}")
