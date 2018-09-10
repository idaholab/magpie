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
    Z = '73     8'
    A = '180.95 16'
    number_fraction = '0.5 0.5'
    displacement_thresholds = '0 0'
    logarithmic_energy_spacing = 1.05
    Emax = 1e7
    uniform_energy_spacing = 0.1
    uniform_energy_spacing_threshold = 10
    execute_on = timestep_end
    displacement_file_base = pc_TaO_energy
    damage_type = ENERGY
  [../]
[]

[Executioner]
  type = Steady
[]
