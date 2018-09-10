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
    Z = '92     6'
    A = '238.03 12.01'
    number_fraction = '0.5 0.5'

    # case 1: (1, 1, 1, 1)
    # displacement_thresholds = '1 1'
    # Ecap = '1 1; 1 1'

    # case 2: (60, 60, 60, 60)
    # displacement_thresholds = '60 60'
    # Ecap = '60 60; 60 60'

    # case 3: (10, 60, 60, 60)
    #displacement_thresholds = '10 60'
    #Ecap = '10 60; 60 60'

    # case 4: (60, 60, 60, 10)
    displacement_thresholds = '60 10'
    Ecap = '60 60; 60 10'

    logarithmic_energy_spacing = 1.05
    Emax = 1e7
    uniform_energy_spacing = 0.1
    uniform_energy_spacing_threshold = 10
    execute_on = timestep_end
    displacement_file_base = pc_UC_net_case4
    damage_type = NET
  [../]
[]

[Executioner]
  type = Steady
[]
