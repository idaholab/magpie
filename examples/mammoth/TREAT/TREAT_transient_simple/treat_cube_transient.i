[Mesh]
 type = GeneratedBIDMesh
 dim = 3
 xmin = 0
 xmax = 400.0
 ymin = 0
 ymax = 400.0
 zmin = 0
 zmax = 400.0
 nx = 8
 ny = 8
 nz = 8
 elem_type = HEX8
 subdomain = '
               1 1 1 1 1 1 1 1
               1 1 1 1 1 1 1 1
               1 1 1 1 1 1 1 1
               1 1 1 1 1 1 1 1
               1 1 1 1 1 1 1 1
               1 1 1 1 1 1 1 1
               1 1 1 1 1 1 1 1
               1 1 1 1 1 1 1 1
               1 1 1 1 1 1 1 1
               1 1 1 1 1 1 1 1
               1 1 1 1 1 1 1 1
               1 1 1 1 1 1 1 1
               1 1 1 1 1 1 1 1
               1 1 1 1 1 1 1 1
               1 1 1 1 1 1 1 1
               1 1 1 1 1 1 1 1
               1 1 1 1 1 1 1 1
               1 1 1 1 1 1 1 1
               1 1 2 2 2 2 1 1
               1 1 2 2 2 2 1 1
               1 1 2 2 2 2 1 1
               1 1 2 2 2 2 1 1
               1 1 1 1 1 1 1 1
               1 1 1 1 1 1 1 1
               1 1 1 1 1 1 1 1
               1 1 1 1 1 1 1 1
               1 1 2 2 2 2 1 1
               1 1 2 2 2 2 1 1
               1 1 2 2 2 2 1 1
               1 1 2 2 2 2 1 1
               1 1 1 1 1 1 1 1
               1 1 1 1 1 1 1 1
               1 1 1 1 1 1 1 1
               1 1 1 1 1 1 1 1
               1 1 2 2 2 2 1 1
               1 1 2 2 2 2 1 1
               1 1 2 2 2 2 1 1
               1 1 2 2 2 2 1 1
               1 1 1 1 1 1 1 1
               1 1 1 1 1 1 1 1
               1 1 1 1 1 1 1 1
               1 1 1 1 1 1 1 1
               1 1 2 2 2 2 1 1
               1 1 2 2 2 2 1 1
               1 1 2 2 2 2 1 1
               1 1 2 2 2 2 1 1
               1 1 1 1 1 1 1 1
               1 1 1 1 1 1 1 1
               1 1 1 1 1 1 1 1
               1 1 1 1 1 1 1 1
               1 1 1 1 1 1 1 1
               1 1 1 1 1 1 1 1
               1 1 1 1 1 1 1 1
               1 1 1 1 1 1 1 1
               1 1 1 1 1 1 1 1
               1 1 1 1 1 1 1 1
               1 1 1 1 1 1 1 1
               1 1 1 1 1 1 1 1
               1 1 1 1 1 1 1 1
               1 1 1 1 1 1 1 1
               1 1 1 1 1 1 1 1
               1 1 1 1 1 1 1 1
               1 1 1 1 1 1 1 1
               1 1 1 1 1 1 1 1
             '
 uniform_refine = 1
[]

[GlobalParams]
  isotopes = 'pseudo'
  densities = 1.0
  plus = true
  isMeter = false
  dbgmat = false
  MGLibObject = BlockLib
  grid_names = 'Tfuel CR_Boron'
  grid_variables = 'Tfuel CR_Boron'
  scalar_flux = 'sflux_g0 sflux_g1 sflux_g2 sflux_g3 sflux_g4 sflux_g5 sflux_g6 sflux_g7 sflux_g8 sflux_g9 sflux_g10'
  evaluate_on = element
[]

[TransportSystems]
  particle = neutron
  equation_type = transient
  G = 11
  VacuumBoundary = '0 1 2 3 4 5'
  [./ts]
    scheme = CFEM-Diffusion
    family = LAGRANGE
    order = FIRST
    n_delay_groups = 6
    fission_source_as_material = true
    collapse_scattering = false
    vacuum_extrapolation_factor = 0.7104

    assemble_fission_jacobian = true
    assemble_scattering_jacobian = true
    assemble_delay_jacobian = true
  [../]
[]

[Variables]
  [./Tfuel]               # Units [=] Kelvin
    order = FIRST
    family = LAGRANGE
    scaling = 1e-4
    [./InitialCondition]
      type = ConstantIC
      value = 300.0
    [../]
  [../]
[]

[Kernels]
  [./HeatConduction]          # Units [=] J/cm^3-sec = Watt/cm^3
    type = HeatConduction
    variable = Tfuel
  [../]
  [./HeatStorage]                          # Heat Storage, units [=] Watt/cm^3
    type = HeatConductionTimeDerivative
    variable =Tfuel
  [../]
  [./HeatSource]            # Heat Source, units [=] Watt/cm^3
    type = CoupledForce
    v = PowerDensity
    block = '2'   # Heat generation only in the core
    variable = Tfuel
  [../]
[]

[BCs]
  [./TempBC]
    type = DirichletBC
    variable = Tfuel
    boundary = '0 1 2 3 4 5'
    value = 300.0
  [../]
[]

[AuxVariables]
  [./CR_Boron]
    order = CONSTANT
    family = MONOMIAL
    [./InitialCondition]
      type = ConstantIC
      value = 0
    [../]
  [../]
  [./PowerDensity]
    order = CONSTANT
    family = MONOMIAL
  [../]

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
[]

[AuxKernels]
  [./PowerDensityCalc]
    type = VectorReactionRate
    variable = PowerDensity
    cross_section = kappa_sigma_fission
    scale_factor = PowerScaling
    execute_on = 'initial linear'
    block ='2'
  [../]
[]

[YAKXSLibraries]
  [./BlockLib]
    type = BaseLibObject
    library_file = FC159.xml
    library_name = FC159_FhomNch
    library_type = MultigroupLibrary
  [../]
[]

[Materials]
  [./graphite]
    type = ConstantDensityFeedbackNeutronicsMaterial
    block = '1'
    material_id = 104
  [../]

  [./fuel]
    type = ConstantDensityFeedbackNeutronicsMaterial
    block ='2'
    material_id = 106
  [../]

  [./ThermalProperties]
    type = HeatConductionMaterial
    temp = Tfuel
    thermal_conductivity_temperature_function = Set_k    # Units [=] W/cm-K
    specific_heat_temperature_function = SetCp     # Units [=] J/g-Kelvin
  [../]

  [./density]
    type = Density
    density = 1.52639336
  [../]
[]

[Functions]
  [./Set_k]
    type = ParsedFunction
    value = '8.178024E-08*t*t - 1.942915E-04*t + 2.816049E-1'   # Units [=] W/cm-K
  [../]
  [./SetCp]
    type = ParsedFunction
    value = '-5.8219E-10*t*t*t - 4.3694E-7*t*t + 2.8369E-3*t -1.009E-2'
  [../]
[]

[Preconditioning]
  [./SMP_jfnk]
    type = SMP
    full = true
  [../]
[]

[Executioner]
  type = Transient
  scheme = crank-nicolson
  start_time = 0.0
  end_time = 2.5 # Obtain period and peak

  l_tol = 1e-3
  l_max_its = 100
  nl_max_its = 200
  nl_abs_tol = 1e-8
  nl_rel_tol = 1e-8
  solve_type = 'PJFNK'
  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart -pc_hypre_boomeramg_max_iter -pc_hypre_boomeramg_tol'
  petsc_options_value = 'hypre boomeramg 100 20 1.0e-6'

  [./TimeStepper]
    type = ConstantDT
    dt = 1.0E-02
    growth_factor = 1.5
  [../]
[]

[UserObjects]
  ## unfortunately we haev to use two UOs because we need to transfer the power density at each point
  ## and PP Transfer only accepts a single PP so we have to split the Multiapps
  [./fission_damage_sampler_P1]
    type = NeutronicsSpectrumSamplerFission
    points = '200 200 200'
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
    execute_on = 'timestep_end'
  [../]

  [./fission_damage_sampler_P2]
    type = NeutronicsSpectrumSamplerFission
    points = '299 200 200'
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
    execute_on = 'timestep_end'
  [../]
[]

[MultiApps]
  [./initial_solve]
    type = FullSolveMultiApp
    app_type = MammothApp
    execute_on = initial
    positions = '0 0 0'
    input_files = treat_cube_steady.i
  [../]

  [./radiation_damage_app_P1]
    type = TransientMultiApp
    input_files = 'TREAT_3D_radiation_damage_sub_P1.i'
    positions = '200 200 200'
    execute_on = 'timestep_end'
  [../]

  [./radiation_damage_app_P2]
    type = TransientMultiApp
    input_files = 'TREAT_3D_radiation_damage_sub_P2.i'
    positions = '299 200 200'
    execute_on = 'timestep_end'
  [../]
[]

[Transfers]
  ## transfers from initial solve
  [./copy_solution]
    type = TransportSystemVariableTransfer
    direction = from_multiapp
    multi_app = initial_solve
    from_transport_system = ts
    to_transport_system = ts
    scale_with_keff = false
  [../]

  [./Copy_pp]
    type = MultiAppPostprocessorTransfer
    direction = from_multiapp
    reduction_type = minimum
    from_postprocessor = PowerScaling
    to_postprocessor = PowerScaling
    multi_app = initial_solve
  [../]

  ## transfers to the marmot sub app 1
  [./radiation_damage_transfer_P1]
    type = MultiAppNeutronicsSpectrumTransfer
    multi_app = radiation_damage_app_P1
    pka_neutronics = neutronics_fission_generator
    radiation_damage_sampler = fission_damage_sampler_P1
    direction = to_multiapp
  [../]

  [./local_power_density_P1]
    type = MultiAppPostprocessorTransfer
    direction = to_multiapp
    reduction_type = minimum
    from_postprocessor = power_density_P1
    to_postprocessor = local_power_density
    multi_app = radiation_damage_app_P1
  [../]

  ## transfers to the marmot sub app 2
  [./radiation_damage_transfer_P2]
    type = MultiAppNeutronicsSpectrumTransfer
    multi_app = radiation_damage_app_P2
    pka_neutronics = neutronics_fission_generator
    radiation_damage_sampler = fission_damage_sampler_P2
    direction = to_multiapp
  [../]

  [./local_power_density_P2]
    type = MultiAppPostprocessorTransfer
    direction = to_multiapp
    reduction_type = minimum
    from_postprocessor = power_density_P2
    to_postprocessor = local_power_density
    multi_app = radiation_damage_app_P2
  [../]
[]

[Postprocessors]
  [./PowerScaling]
    type = Receiver
    outputs = none
  [../]

  [./Power]
    type = ElementIntegralVariablePostprocessor
    variable = PowerDensity
    execute_on = 'initial linear'
    outputs = all
  [../]

  [./avg_Tfuel]
    # Units [=] Kelvin
    type = ElementAverageValue
    execute_on = 'initial timestep_end'
    variable = Tfuel
    block = '2'
    outputs = all
  [../]

  [./max_Tfuel]
    # (Not spatially Dependent: Constant Value)
    type = ElementExtremeValue
    execute_on = 'initial timestep_end'
    variable = Tfuel
    value_type = max
    block = '2'
    outputs = all
  [../]

  [./power_density_P1]
    type = PointValue
    variable = PowerDensity
    point = '200 200 200'
    execute_on = 'initial timestep_end'
  [../]

  [./power_density_P2]
    type = PointValue
    variable = PowerDensity
    point = '299 200 200'
    execute_on = 'initial timestep_end'
  [../]
[]

[Outputs]
  csv = true
  print_linear_residuals = true
  print_perf_log = true
  exodus = true
[]
