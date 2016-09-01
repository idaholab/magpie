[Mesh]
  file = 3x3-Standard.e
  parallel_type = REPLICATED  # must be serial with MatID
  uniform_refine = 0
[]

[Variables]
  [./u]
  [../]
[]

[Kernels]
  [./diff]
    type = Diffusion
    variable = u
  [../]
[]

[BCs]
  [./left]
    type = DirichletBC
    variable = u
    boundary = 1
    value = 0
  [../]

  [./right]
    type = DirichletBC
    variable = u
    boundary = 2
    value = 1
  [../]
[]


[AuxVariables]
  [./N92235]
    order = CONSTANT
    family = MONOMIAL
    [./InitialCondition]
      type = ConstantIC
      value = 8.68490 # note: multiplied by 1e6
    [../]
  [../]

  [./N92238]
    order = CONSTANT
    family = MONOMIAL
    [./InitialCondition]
      type = ConstantIC
      value = 6.29670E-01 # note: multiplied by 1e6
    [../]
  [../]

  ## scalar fluxes
  [./sflux_g0]
    family = LAGRANGE
    order = FIRST
  [../]
  [./sflux_g1]
    family = LAGRANGE
    order = FIRST
  [../]
  [./sflux_g2]
    family = LAGRANGE
    order = FIRST
  [../]
  [./sflux_g3]
    family = LAGRANGE
    order = FIRST
  [../]
  [./sflux_g4]
    family = LAGRANGE
    order = FIRST
  [../]
  [./sflux_g5]
    family = LAGRANGE
    order = FIRST
  [../]
  [./sflux_g6]
    family = LAGRANGE
    order = FIRST
  [../]
  [./sflux_g7]
    family = LAGRANGE
    order = FIRST
  [../]
  [./sflux_g8]
    family = LAGRANGE
    order = FIRST
  [../]
  [./sflux_g9]
    family = LAGRANGE
    order = FIRST
  [../]
  [./sflux_g10]
    family = LAGRANGE
    order = FIRST
  [../]
[]

[UserObjects]
  [./fission_damage_sampler]
    type = NeutronicsSpectrumSamplerFission
    points = '0 0 143.552375'
    #points = '0 0 143.552375
    #          0 0 163.79425
    #          0 0 183.79425'
    energy_group_boundaries = '0.0 2.001e-8 4.7302e-08 7.6497e-8 2.09610e-7 6.25e-7
                               8.10003e-6 1.327e-4 3.4811e-3 1.1562e-1 3.3287 20'
    target_isotope_names = 'U235 U238'
    number_densities = 'N92235 N92238'
    scalar_fluxes = 'sflux_g10 sflux_g9 sflux_g8 sflux_g7 sflux_g6 sflux_g5 sflux_g4 sflux_g3 sflux_g2 sflux_g1 sflux_g0'
    # each row are the recoil cross sections for a single isotope && all group (g=0,1)
    # I also set a lot of the tiny thermal fission XS for 92238 to zero (slowest five groups E < 6.25e-7) because we have no product data for them
    fission_cross_sections = '
                              9.4511e2 5.12074e2 3.52511e2 2.37127e2 1.4114e2 3.14337e1  4.48499e1  1.18495e1  2.55139    1.24034   1.17273
                              0        0         0         0         0        2.33593e-5 7.24493e-5 5.65586e-4 6.63056e-5 1.1122e-1 6.01374e-1
                            '
    #fission_cross_sections = '
    #                          9.4511e2 5.12074e2 3.52511e2 2.37127e2 1.4114e2 3.14337e1 4.48499e1 1.18495e1 2.55139 1.24034 1.17273
    #                          2.63065e-5 1.49392e-5 1.09581e-5 8.08171e-6 4.72438e-6 2.33593e-5 7.24493e-5 5.65586e-4 6.63056e-5 1.1122e-1 6.01374e-1
    #                          9.4511e2 5.12074e2 3.52511e2 2.37127e2 1.4114e2 3.14337e1 4.48499e1 1.18495e1 2.55139 1.24034 1.17273
    #                          2.63065e-5 1.49392e-5 1.09581e-5 8.08171e-6 4.72438e-6 2.33593e-5 7.24493e-5 5.65586e-4 6.63056e-5 1.1122e-1 6.01374e-1
    #                          9.4511e2 5.12074e2 3.52511e2 2.37127e2 1.4114e2 3.14337e1 4.48499e1 1.18495e1 2.55139 1.24034 1.17273
    #                          2.63065e-5 1.49392e-5 1.09581e-5 8.08171e-6 4.72438e-6 2.33593e-5 7.24493e-5 5.65586e-4 6.63056e-5 1.1122e-1 6.01374e-1
    #                        '
    execute_on = 'initial nonlinear'
  [../]
[]

[MultiApps]
  [./eigenvalue_solve]
    type = FullSolveMultiApp
    input_files = 'TREAT_eigenvalue_assembly_11g.i '
    execute_on = 'initial'
  [../]

  [./radiation_damage_app]
    type = FullSolveMultiApp
    input_files = 'TREAT_2D_radiation_damage_sub.i'
    #input_files = 'TREAT_3D_radiation_damage_sub.i'
    positions = '0 0 143.552375'
    #positions = '0 0 143.552375
    #             0 0 163.79425
    #             0 0 183.79425'
    execute_on = 'nonlinear'
  []
[]

[Transfers]
  [./transfer_g0]
    type = MultiAppProjectionTransfer
    multi_app = eigenvalue_solve
    direction = FROM_MULTIAPP
    variable = sflux_g0
    source_variable = sflux_g0
  [../]
  [./transfer_g1]
    type = MultiAppProjectionTransfer
    multi_app = eigenvalue_solve
    direction = FROM_MULTIAPP
    variable = sflux_g1
    source_variable = sflux_g1
  [../]
  [./transfer_g2]
    type = MultiAppProjectionTransfer
    multi_app = eigenvalue_solve
    direction = FROM_MULTIAPP
    variable = sflux_g2
    source_variable = sflux_g2
  [../]
  [./transfer_g3]
    type = MultiAppProjectionTransfer
    multi_app = eigenvalue_solve
    direction = FROM_MULTIAPP
    variable = sflux_g3
    source_variable = sflux_g3
  [../]
  [./transfer_g4]
    type = MultiAppProjectionTransfer
    multi_app = eigenvalue_solve
    direction = FROM_MULTIAPP
    variable = sflux_g4
    source_variable = sflux_g4
  [../]
  [./transfer_g5]
    type = MultiAppProjectionTransfer
    multi_app = eigenvalue_solve
    direction = FROM_MULTIAPP
    variable = sflux_g5
    source_variable = sflux_g5
  [../]
  [./transfer_g6]
    type = MultiAppProjectionTransfer
    multi_app = eigenvalue_solve
    direction = FROM_MULTIAPP
    variable = sflux_g6
    source_variable = sflux_g6
  [../]
  [./transfer_g7]
    type = MultiAppProjectionTransfer
    multi_app = eigenvalue_solve
    direction = FROM_MULTIAPP
    variable = sflux_g7
    source_variable = sflux_g7
  [../]
  [./transfer_g8]
    type = MultiAppProjectionTransfer
    multi_app = eigenvalue_solve
    direction = FROM_MULTIAPP
    variable = sflux_g8
    source_variable = sflux_g8
  [../]
  [./transfer_g9]
    type = MultiAppProjectionTransfer
    multi_app = eigenvalue_solve
    direction = FROM_MULTIAPP
    variable = sflux_g9
    source_variable = sflux_g9
  [../]
  [./transfer_g10]
    type = MultiAppProjectionTransfer
    multi_app = eigenvalue_solve
    direction = FROM_MULTIAPP
    variable = sflux_g10
    source_variable = sflux_g10
  [../]
  [./radiation_damage_transfer]
    type = MultiAppNeutronicsSpectrumTransfer
    multi_app = radiation_damage_app
    pka_neutronics = neutronics_fission_generator
    radiation_damage_sampler = fission_damage_sampler
    direction = to_multiapp
  [../]
[]

[Executioner]
  type = Transient
  num_steps = 1
[]

[Outputs]
  file_base = TREAT_damage_master
  exodus = true
  print_linear_residuals = true
  print_perf_log = true
  csv = true
[]
