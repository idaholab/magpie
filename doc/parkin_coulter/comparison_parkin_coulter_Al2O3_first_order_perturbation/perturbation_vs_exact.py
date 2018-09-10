import numpy as np
import matplotlib.pyplot as plt

base_file = "parkin_coulter_Al2O3_base.csv"
mod_file = "parkin_coulter_Al2O3_mod.csv"

derivative_file = "parkin_coulter_Al2O3_derivative.csv"

base = np.loadtxt(base_file, delimiter = ",", skiprows = 1)
mod = np.loadtxt(mod_file, delimiter = ",", skiprows = 1)
derivative = np.loadtxt(derivative_file, delimiter = ",", skiprows = 1)

plt.xlabel('Energy (eV)')
plt.ylabel('Net displacements (gij)')
plt.semilogx(base[:, 0], base[:, 1], 'b-', label = "Al - Al")
plt.semilogx(base[:, 0], base[:, 2], 'g-', label = "Al - O")
plt.semilogx(base[:, 0], base[:, 3], 'r-', label = "O - Al")
plt.semilogx(base[:, 0], base[:, 4], 'k-', label = "O - O")

plt.semilogx(mod[:, 0], mod[:, 1], 'bx', label = "Al - Al, mod")
plt.semilogx(mod[:, 0], mod[:, 2], 'gx', label = "Al - O, mod")
plt.semilogx(mod[:, 0], mod[:, 3], 'rx', label = "O - Al, mod")
plt.semilogx(mod[:, 0], mod[:, 4], 'kx', label = "O - O, mod")

# compute changes using first order perturbations
delta = [-0.1, 0.1]
fop = mod
for i in range(len(fop)):
    fop[i, 1:] = 0.0

for i in range(len(fop)):
    for j in range(1, len(fop[i, :])):
        fop[i, j] = base[i, j]
        for l in [1, 2]:
            col = (j - 1) * 2 + l
            df = delta[l - 1]
            fop[i, j] = fop[i, j] + df * derivative[i, col]

plt.semilogx(fop[:, 0], fop[:, 1], 'b-.', label = "Al - Al, fop")
plt.semilogx(fop[:, 0], fop[:, 2], 'g-.', label = "Al - O, fop")
plt.semilogx(fop[:, 0], fop[:, 3], 'r-.', label = "O - Al, fop")
plt.semilogx(fop[:, 0], fop[:, 4], 'k-.', label = "O - O, fop")
plt.legend(loc=0)
plt.savefig("NRT_first_order_perturbation_theory.pdf")
plt.show()
