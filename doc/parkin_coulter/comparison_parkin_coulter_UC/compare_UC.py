import numpy as np

case1 = np.loadtxt("pc_UC_net_case1.csv", skiprows = 1, delimiter = ",")
case2 = np.loadtxt("pc_UC_net_case2.csv", skiprows = 1, delimiter = ",")
case3 = np.loadtxt("pc_UC_net_case3.csv", skiprows = 1, delimiter = ",")
case4 = np.loadtxt("pc_UC_net_case4.csv", skiprows = 1, delimiter = ",")
nu = np.loadtxt("pc_UC_energy.csv", skiprows = 1, delimiter = ",")
energy = 1e7

Ed1 = [1.0, 1.0]
Ed2 = [60.0, 60.0]
Ed3 = [10.0, 60.0]
Ed4 = [60.0, 10.0]

case1_eff = [0 for j in range(4)]
case2_eff = [0 for j in range(4)]
case3_eff = [0 for j in range(4)]
case4_eff = [0 for j in range(4)]

case1_err = [0 for j in range(4)]
case2_err = [0 for j in range(4)]
case3_err = [0 for j in range(4)]
case4_err = [0 for j in range(4)]

case1_ref = [0.454, 0.266, 0.461, 0.258]
case2_ref = [0.545, 0.237, 0.548, 0.205]
case3_ref = [0.498, 0.230, 0.508, 0.199]
case4_ref = [0.532, 0.261, 0.535, 0.236]

for i in range(2):
    for j in range(2):
        n = 2 * i + j + 1
        case1_eff[n - 1] = (case1[len(case1) - 1, n] - 1.0) * Ed1[j] / (0.5 * nu[len(nu) - 1, i + 1])
        case2_eff[n - 1] = (case2[len(case2) - 1, n]) * Ed2[j] / (0.5 * nu[len(nu) - 1, i + 1])
        case3_eff[n - 1] = (case3[len(case3) - 1, n]) * Ed3[j] / (0.5 * nu[len(nu) - 1, i + 1])
        case4_eff[n - 1] = (case4[len(case4) - 1, n] - 1.0) * Ed4[j] / (0.5 * nu[len(nu) - 1, i + 1])

        case1_err[n - 1] = abs(1.0 - case1_eff[n - 1] / case1_ref[n - 1]) * 100.0
        case2_err[n - 1] = abs(1.0 - case2_eff[n - 1] / case2_ref[n - 1]) * 100.0
        case3_err[n - 1] = abs(1.0 - case3_eff[n - 1] / case3_ref[n - 1]) * 100.0
        case4_err[n - 1] = abs(1.0 - case4_eff[n - 1] / case4_ref[n - 1]) * 100.0


print "----- Magpie results -----"
print "%8.3f %8.3f %8.3f %8.3f" % (case1_eff[0], case1_eff[1], case1_eff[2], case1_eff[3])
print "%8.3f %8.3f %8.3f %8.3f" % (case2_eff[0], case2_eff[1], case2_eff[2], case2_eff[3])
print "%8.3f %8.3f %8.3f %8.3f" % (case3_eff[0], case3_eff[1], case3_eff[2], case3_eff[3])
print "%8.3f %8.3f %8.3f %8.3f" % (case4_eff[0], case4_eff[1], case4_eff[2], case4_eff[3])

print "----- Relative difference w/ PC 1980 -----"
print "%8.2f %8.2f %8.2f %8.2f" % (case1_err[0], case1_err[1], case1_err[2], case1_err[3])
print "%8.2f %8.2f %8.2f %8.2f" % (case2_err[0], case2_err[1], case2_err[2], case2_err[3])
print "%8.2f %8.2f %8.2f %8.2f" % (case3_err[0], case3_err[1], case3_err[2], case3_err[3])
print "%8.2f %8.2f %8.2f %8.2f" % (case4_err[0], case4_err[1], case4_err[2], case4_err[3])
