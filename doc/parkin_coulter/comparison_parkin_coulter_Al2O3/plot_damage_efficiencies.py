import matplotlib.pyplot as plt
import numpy as np
import sys

def nuFit(energy, params):
    return energy / (1 + params[0] * energy**0.15 + params[1] * energy**0.75 + params[2] * energy)

gij = np.loadtxt("pc_Al2O3_net.csv", skiprows = 1, delimiter = ",")
nui = np.loadtxt("pc_Al2O3_energy.csv",\
                 skiprows = 1, delimiter = ",")

Edisp = [18.0, 72.0]
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

    eff_11.append((gij[j, 1] - 1.0) * Edisp[0] / (0.4 * nu1))
    eff_12.append(gij[j, 2] * Edisp[1] / (0.6 * nu1))
    eff_21.append(gij[j, 3] * Edisp[0] / (0.4 * nu2))
    eff_22.append((gij[j, 4] - 1.0) * Edisp[1] / (0.6 * nu2))

pk = np.loadtxt("PC_paper_Al2O3.dat", skiprows = 1)

plt.xlabel('Energy (eV)')
plt.ylabel('Displacement efficiency (kij)')
plt.semilogx(gij[:, 0], eff_11, 'b-', label = "Al - Al")
plt.semilogx(gij[:, 0], eff_12, 'g-', label = "Al - O")
plt.semilogx(gij[:, 0], eff_21, 'r-', label = "O - Al")
plt.semilogx(gij[:, 0], eff_22, 'k-', label = "O - O")
plt.semilogx(pk[:, 0], pk[:, 1], 'bs')
plt.semilogx(pk[:, 0], pk[:, 2], 'gs')
plt.semilogx(pk[:, 0], pk[:, 3], 'rs')
plt.semilogx(pk[:, 0], pk[:, 4], 'ks')


plt.xlim((50, 1e7))
plt.ylim((0, 0.5))
plt.legend(loc=0)

plt.savefig("comparison_PC_1981_Magpie_Al2O3.pdf")
plt.show()
