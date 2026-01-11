import numpy as np
import matplotlib.pyplot as plt

# Antoine constants  for Benzene–Toluene system
A1, B1, C1 = 6.90565, 1211.033, 220.790   # Benzene
A2, B2, C2 = 6.95464, 1344.800, 219.480   # Toluene

# User Input
P_total = float(input("Enter total pressure in mmHg : "))
T_min = float(input("Enter minimum temperature in °C : "))
T_max = float(input("Enter maximum temperature in °C : "))
n_points = int(input("Enter number of calculation points : "))

# Temperature range
T = np.linspace(T_min, T_max, n_points)

# Saturation pressures using Antoine equation
P1_sat = 10**(A1 - B1 / (C1 + T))  # Benzene
P2_sat = 10**(A2 - B2 / (C2 + T))  # Toluene

# Liquid phase composition
x1 = np.linspace(0, 1, n_points)
x2 = 1 - x1

# Average saturation pressures 
P1_avg = np.mean(P1_sat)
P2_avg = np.mean(P2_sat)

# Vapor phase composition using Raoult’s law
y1 = (x1 * P1_avg) / P_total


y1 = np.clip(y1, 0, 1)

# Plot VLE x–y diagram
plt.figure(figsize=(7, 7))
plt.plot(x1, y1, label="VLE Curve (Raoult’s Law)", linewidth=2)
plt.plot(x1, x1, '--', label="x = y Line")
plt.xlabel("Liquid mole fraction of Benzene (x₁)")
plt.ylabel("Vapor mole fraction of Benzene (y₁)")
plt.title("VLE Diagram for Benzene–Toluene System")
plt.legend()
plt.grid(True)
plt.show()
