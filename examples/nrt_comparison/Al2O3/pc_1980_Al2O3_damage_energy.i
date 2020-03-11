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
    displacement_thresholds = '0 0'
    logarithmic_energy_spacing = 1.05
    Emax = 1e7
    uniform_energy_spacing = 0.1
    uniform_energy_spacing_threshold = 10
    displacement_file_base = pc_Al2O3_energy
    damage_type = ENERGY
  [../]
[]

[Executioner]
  type = Steady
[]
