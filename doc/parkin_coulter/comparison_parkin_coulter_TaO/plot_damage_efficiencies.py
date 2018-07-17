import matplotlib.pyplot as plt
import numpy as np
import sys

def nuFit(energy, params):
    return energy / (1 + params[0] * energy**0.15 + params[1] * energy**0.75 + params[2] * energy)

gij = np.loadtxt("pc_TaO_net.csv", skiprows = 1, delimiter = ",")
nui = np.loadtxt("pc_TaO_energy.csv",\
                 skiprows = 1, delimiter = ",")

Edisp = [60, 60]
paramsAl = [1.038e-1, 5.917e-6, 3.887e-6]
paramsO = [1.108e-1, 1.129e-4, 9.258e-6]

eff_11 = []
eff_12 = []
eff_21 = []
eff_22 = []
for j in range(len(gij)):
    energy = gij[j, 0]

    # interpolate damage energy values
    #nu1 = nuFit(energy, paramsAl)
    #nu2 = nuFit(energy, paramsO) #

    nu1 = np.interp(energy, nui[:, 0], nui[:, 1])
    nu2 = np.interp(energy, nui[:, 0], nui[:, 2])

    #print energy, nu1 / energy, nu2 / energy

    eff_11.append((gij[j, 1] - 1.0) * Edisp[0] / (0.5 * nu1))
    eff_12.append(gij[j, 2] * Edisp[1] / (0.5 * nu1))
    eff_21.append(gij[j, 3] * Edisp[0] / (0.5 * nu2))
    eff_22.append((gij[j, 4] - 1.0) * Edisp[1] / (0.5 * nu2))

Ta_Ta = np.loadtxt("Ta-Ta.csv", skiprows = 1, delimiter = ",")
Ta_O = np.loadtxt("Ta-O.csv", skiprows = 1, delimiter = ",")
O_Ta = np.loadtxt("O-Ta.csv", skiprows = 1, delimiter = ",")
O_O = np.loadtxt("O-O.csv", skiprows = 1, delimiter = ",")

plt.xlabel('Energy (eV)')
plt.ylabel('Displacement efficiency (kij)')
plt.semilogx(gij[:, 0], eff_11, 'b-', label = "Ta - Ta")
plt.semilogx(gij[:, 0], eff_12, 'g-', label = "Ta - O")
plt.semilogx(gij[:, 0], eff_21, 'r-', label = "O - Ta")
plt.semilogx(gij[:, 0], eff_22, 'k-', label = "O - O")
plt.semilogx(Ta_Ta[:, 0], Ta_Ta[:, 1], 'b--')
plt.semilogx(Ta_O[:, 0], Ta_O[:, 1], 'g--')
plt.semilogx(O_Ta[:, 0], O_Ta[:, 1], 'r--')
plt.semilogx(O_O[:, 0], O_O[:, 1], 'k--')


plt.xlim((50, 1e7))
plt.ylim((0, 0.8))
plt.legend(loc=0)

plt.savefig("comparison_PC_1981_Magpie_TaO.pdf")
plt.show()
