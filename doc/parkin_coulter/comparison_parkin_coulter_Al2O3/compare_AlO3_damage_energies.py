import numpy as np

def nuFit(energy, params):
    return energy / (1 + params[0] * energy**0.15 + params[1] * energy**0.75 + params[2] * energy)

data = np.loadtxt("pc_Al2O3_energy.csv", skiprows = 1, delimiter = ",")

energy_points = [10, 100, 1000, 1e4, 1e5, 1e6, 1e7]
Al_eff = [0.878, 0.83, 0.766, 0.664, 0.441, 0.132, 0.0196]
O_eff = [0.868, 0.817, 0.745, 0.613, 0.318, 0.0682, 0.0088]

paramsAl = [1.03736e-1, 5.91664e-6, 3.887e-6]
paramsO = [1.10761e-1, 1.12891e-4, 9.25824e-6]

if True:
    Al_eff = []
    O_eff = []
    for e in energy_points:
        Al_eff.append(nuFit(e, paramsAl) / e)
        O_eff.append(nuFit(e, paramsO) / e)


for j in range(len(energy_points)):
    al_diff = abs(1.0 - np.interp(energy_points[j], data[:, 0], data[:, 1]) / energy_points[j] / Al_eff[j]) * 100
    o_diff = abs(1.0 - np.interp(energy_points[j], data[:, 0], data[:, 2]) / energy_points[j] / O_eff[j]) * 100
    string = "%8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f \\\\" % (energy_points[j],\
             Al_eff[j], np.interp(energy_points[j], data[:, 0], data[:, 1]) / energy_points[j], al_diff,\
             O_eff[j], np.interp(energy_points[j], data[:, 0], data[:, 2]) / energy_points[j], o_diff)
    print string
