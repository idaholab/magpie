#
# Case from Parkin Coulter "Damage Energy Functions in Polyatomic Materials",
# JNM 88, (1980).
# U_x Zr_{1-x} C; x = 0.02 & 0.05
#
[Mesh]
 type = GeneratedMesh
 dim = 1
 xmin = 0
 xmax = 1
 nx = 2
[]

[Problem]
  type = FEProblem
  kernel_coverage_check = false
[]

[Variables]
  [./test_var]
  [../]
[]

[UserObjects]
  [./parkin_coulter]
    type = PolyatomicRecoil
    Z = '92      40     6'
    A = '238.028 91.224 12.011'
    # case 1: x = 0.02
    # number_fraction = '0.01 0.49 0.5'
    # case 2: x = 0.5
    # number_fraction = '0.25 0.25 0.5'
    # case 3: x = 0.98
    number_fraction = '0.49 0.01 0.5'

    displacement_thresholds = '0 0 0'
    logarithmic_energy_spacing = 1.05
    Emax = 1e7
    uniform_energy_spacing = 0.1
    uniform_energy_spacing_threshold = 10
    displacement_file_base = parkin_coulter_UxZr1-xC_case3
    damage_type = ENERGY
  [../]
[]

[Executioner]
  type = Steady
[]
