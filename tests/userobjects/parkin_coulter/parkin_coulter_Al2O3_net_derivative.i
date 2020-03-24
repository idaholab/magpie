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
    Z = '13  8'
    A = '27  16'
    number_fraction = '0.4 0.6'
    displacement_thresholds = '18 72'
    Ecap = '18 45; 45 72'
    logarithmic_energy_spacing = 1.15
    Emax = 500
    uniform_energy_spacing = 0.2
    uniform_energy_spacing_threshold = 10
    displacement_file_base = parkin_coulter_Al2O3
    damage_type = NET
  [../]
[]

[Executioner]
  type = Steady
[]
