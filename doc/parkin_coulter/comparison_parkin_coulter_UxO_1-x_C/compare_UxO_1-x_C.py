import numpy as np

case1 = np.loadtxt("parkin_coulter_UxZr1-xC_case1.csv", skiprows = 1, delimiter = ",")
case2 = np.loadtxt("parkin_coulter_UxZr1-xC_case2.csv", skiprows = 1, delimiter = ",")
case3 = np.loadtxt("parkin_coulter_UxZr1-xC_case3.csv", skiprows = 1, delimiter = ",")

energy_points = [1.12e1, 1.07e2, 1.02e3, 1.15e4, 1.1e5, 1.05e6, 1e7]

print "----- Magpie Results -----"
for j in range(len(energy_points)):
    Uc1 = np.interp(energy_points[j], case1[:, 0], case1[:, 1]) / energy_points[j]
    Zrc1 = np.interp(energy_points[j], case1[:, 0], case1[:, 2]) / energy_points[j]
    Cc1 = np.interp(energy_points[j], case1[:, 0], case1[:, 3]) / energy_points[j]
    Uc2 = np.interp(energy_points[j], case2[:, 0], case2[:, 1]) / energy_points[j]
    Zrc2 = np.interp(energy_points[j], case2[:, 0], case2[:, 2]) / energy_points[j]
    Cc2 = np.interp(energy_points[j], case2[:, 0], case2[:, 3]) / energy_points[j]
    Uc3 = np.interp(energy_points[j], case3[:, 0], case3[:, 1]) / energy_points[j]
    Zrc3 = np.interp(energy_points[j], case3[:, 0], case3[:, 2]) / energy_points[j]
    Cc3 = np.interp(energy_points[j], case3[:, 0], case3[:, 3]) / energy_points[j]

    string = "%8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f \\\\" %\
             (energy_points[j], Uc1, Uc2, Uc3, Zrc1, Zrc2, Zrc3, Cc1, Cc2, Cc3)
    print string

print
print "----- Relative Difference to PC 1980 -----"

U_case_1 = [0.894, 0.855, 0.803, 0.732, 0.639, 0.498, 0.28]
U_case_2 = [0.889, 0.848, 0.796, 0.722, 0.631, 0.501, 0.297]
U_case_3 = [0.886, 0.846, 0.793, 0.720, 0.631, 0.510, 0.314]

Zr_case_1 = [0.878, 0.835, 0.777, 0.696, 0.582, 0.387, 0.128]
Zr_case_2 = [0.870, 0.824, 0.764, 0.680, 0.570, 0.388, 0.137]
Zr_case_3 = [0.863, 0.815, 0.754, 0.670, 0.564, 0.394, 0.146]

C_case_1 = [0.820, 0.758, 0.672, 0.507, 0.235, 0.050, 0.007]
C_case_2 = [0.797, 0.730, 0.638, 0.472, 0.222, 0.051, 0.007]
C_case_3 = [0.770, 0.697, 0.600, 0.436, 0.211, 0.052, 0.008]

for j in range(len(energy_points)):
    Uc1 = np.interp(energy_points[j], case1[:, 0], case1[:, 1]) / energy_points[j]
    Uc1 = abs(1 - Uc1 / U_case_1[j]) * 100
    Zrc1 = np.interp(energy_points[j], case1[:, 0], case1[:, 2]) / energy_points[j]
    Zrc1 = abs(1 - Zrc1 / Zr_case_1[j]) * 100
    Cc1 = np.interp(energy_points[j], case1[:, 0], case1[:, 3]) / energy_points[j]
    Cc1 = abs(1 - Cc1 / C_case_1[j]) * 100
    Uc2 = np.interp(energy_points[j], case2[:, 0], case2[:, 1]) / energy_points[j]
    Uc2 = abs(1 - Uc2 / U_case_2[j]) * 100
    Zrc2 = np.interp(energy_points[j], case2[:, 0], case2[:, 2]) / energy_points[j]
    Zrc2 = abs(1 - Zrc2 / Zr_case_2[j]) * 100
    Cc2 = np.interp(energy_points[j], case2[:, 0], case2[:, 3]) / energy_points[j]
    Cc2 = abs(1 - Cc2 / C_case_2[j]) * 100
    Uc3 = np.interp(energy_points[j], case3[:, 0], case3[:, 1]) / energy_points[j]
    Uc3 = abs(1 - Uc3 / U_case_3[j]) * 100
    Zrc3 = np.interp(energy_points[j], case3[:, 0], case3[:, 2]) / energy_points[j]
    Zrc3 = abs(1 - Zrc3 / Zr_case_3[j]) * 100
    Cc3 = np.interp(energy_points[j], case3[:, 0], case3[:, 3]) / energy_points[j]
    Cc3 = abs(1 - Cc3 / C_case_3[j]) * 100

    string = "%8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f \\\\" %\
             (energy_points[j], Uc1, Uc2, Uc3, Zrc1, Zrc2, Zrc3, Cc1, Cc2, Cc3)
    print string
