[Mesh]
  file = 3x3-Standard.e
  distribution = SERIAL  # must be serial with MatID
  uniform_refine = 0
[]

[Debug]
  show_actions = 0
  show_material_props = 0
  show_neutronics_material_coverage = 1
  print_block_volume = 0
[]

[GlobalParams]
  isotopes = 'pseudo'
  forDiffusion = true
  densities = 1.0
  ngroup = 11
  plus = true
  isMeter = false
  dbgmat = false
  grid = '1'
  grid_names = 'Tfuel'
  material_id_variable = MatID
  library_file = SPH.xml
  scalar_flux = 'sflux_g0 sflux_g1 sflux_g2 sflux_g3 sflux_g4 sflux_g5 sflux_g6 sflux_g7 sflux_g8 sflux_g9 sflux_g10'
  material_ids = '1 2 3 4 5 6 7 8 9 10 11 12 13'
[]

[TransportSystems]
 particle = neutron
 equation_type = eigenvalue

 G = 11
 ReflectingBoundary ='1 2 3 4'
 VacuumBoundary = '5 6'

 [./diff]
   scheme = CFEM-Diffusion
   n_delay_groups = 0
   family = LAGRANGE
   order = FIRST
   #fission_source_as_material = true
   #diffusion_coefficient_type = tensor
   diffusion_coefficient_scheme = local
   balance_table = true
[../]
[]

[AuxVariables]
  [./PowerDensity]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./Tfuel]
    order = CONSTANT
    family = MONOMIAL
      [./InitialCondition]
        type = ConstantIC
        value = 293.6  # Tabulated value
      [../]
  [../]

  [./MatID]
   order = CONSTANT
   family = MONOMIAL
   initial_from_file_var = material_id
   initial_from_file_timestep = 1
  [../]
[]

[AuxKernels]
  [./Power]
    #type = FissionSource
    type = VectorReactionRate
    variable = PowerDensity
    #nusigf = kappa_sigma_fission
    cross_section = kappa_sigma_fission
    block = '11000'
  # A's, B's, 6100=6
  [../]
[]

[UserObjects]

[./Balance]
 type = FluxCartesianCoreMap
 transport_system = diff
 block = '11000'
 print= 'block material region'
 power_map_from = kappa_sigma_f
 print_groupflux = 1
 print_assemblywise_fluxes = true
 # grid_coord are for end pts for each assembly
 output_in = 'Balance.out'
 execute_on = 'timestep_end'
[../]

[./fission_damage_sampler]
  type = RadiationDamageFission
  points = '0 0 143.552375
            0 0 163.79425
            0 0 183.79425'
  energy_group_boundaries = '1.0e-11 2.001e-8 4.7302e-08 7.6497e-8 2.09610e-7 6.25e-7
                             8.10003e-6 1.327e-4 3.4811e-3 1.1562e-1 3.3287 20'
  target_isotope_names = 'U235 U238'
  number_densities = 'N92235 N92238'
  scalar_fluxes = 'sflux_g10 sflux_g9 sflux_g8 sflux_g7 sflux_g6 sflux_g5 sflux_g4 sflux_g3 sflux_g2 sflux_g1 sflux_g0'
  # each row are the recoil cross sections for a single isotope && all group (g=0,1)
  fission_cross_sections = '
                            9.4511e2 5.12074e2 3.52511e2 2.37127e2 1.4114e2 3.14337e1 4.48499e1 1.18495e1 2.55139 1.24034 1.17273
                            2.63065e-5 1.49392e-5 1.09581e-5 8.08171e-6 4.72438e-6 2.33593e-5 7.24493e-5 5.65586e-4 6.63056e-5 1.1122e-1 6.01374e-1
                            9.4511e2 5.12074e2 3.52511e2 2.37127e2 1.4114e2 3.14337e1 4.48499e1 1.18495e1 2.55139 1.24034 1.17273
                            2.63065e-5 1.49392e-5 1.09581e-5 8.08171e-6 4.72438e-6 2.33593e-5 7.24493e-5 5.65586e-4 6.63056e-5 1.1122e-1 6.01374e-1
                            9.4511e2 5.12074e2 3.52511e2 2.37127e2 1.4114e2 3.14337e1 4.48499e1 1.18495e1 2.55139 1.24034 1.17273
                            2.63065e-5 1.49392e-5 1.09581e-5 8.08171e-6 4.72438e-6 2.33593e-5 7.24493e-5 5.65586e-4 6.63056e-5 1.1122e-1 6.01374e-1
                            9.4511e2 5.12074e2 3.52511e2 2.37127e2 1.4114e2 3.14337e1 4.48499e1 1.18495e1 2.55139 1.24034 1.17273
                            2.63065e-5 1.49392e-5 1.09581e-5 8.08171e-6 4.72438e-6 2.33593e-5 7.24493e-5 5.65586e-4 6.63056e-5 1.1122e-1 6.01374e-1
                          '
  execute_on = initial
[../]


[./Power]
  type = VariableCartesianCoreMap
  variables='PowerDensity'
  regular_grid = true
  grid_coord_x = '-5.08 5.08'
  grid_coord_y = '-5.08 5.08'
  grid_coord_z = '62.82675 72.82675 82.82675 92.9476875 103.068625 113.1895625 123.3105 133.4314375 143.552375 153.6733125 163.79425 173.79425 183.79425'
  output_in = 'Power.out'
  execute_on = 'timestep_end'
[../]

[]


[YAKXSLibraries]
  [./A0]
     type = BaseLibObject
     library_name = A0_A11000
     library_type = MultigroupLibrary
  [../]
[]

[MultiApps]
  [./radiation_damage_app]
    type = FullSolveMultiApp
    input_files = TREAT_sub.i
    points = '0 0 143.552375
              0 0 163.79425
              0 0 183.79425'
    execute_on = nonlinear
  []
[]

[Transfers]
  [./radiation_damage_transfer]
    type = MultiAppRadiationDamageTransfer
    multi_app = radiation_damage_app
    pka_neutronics = neutronics_fission_generator
    radiation_damage_sampler = fission_damage_sampler
    direction = to_multiapp
  [../]
[]

[Materials]
 [./Mat_1]
    type = MixedMatIDNeutronicsMaterial
    block = '11000'
    multigroup_library_object = A0
 [../]
[]

[Executioner]
  type = NonlinearEigen
  solve_type = 'PJFNK'
#petsc_options_iname =  '-pc_type -pc_hypre_type -pc_hypre_boomeramg_tol -pc_hypre_boomeramg_max_iter -ksp_gmres_restart'
#petsc_options_value =  'hypre boomeramg 1e-4 10 50'
petsc_options_iname =  '-pc_type -pc_hypre_type -ksp_gmres_restart'
petsc_options_value =  'hypre boomeramg 50'

  l_tol = 1e-4
  free_power_iterations = 4.
  source_abs_tol = 1e-8
  normal_factor = 100.0
  normalization = TotalPower
[]

[VectorPostprocessors]
  [./Fission Source]
    type = PointValueSampler
    variable = 'fission_source'
    points = '0 0 143.552375 0 0 163.79425 0 0 183.79425'
    sort_by = x
    execute_on = timestep_end
  [../]
  [./Scalar Flux (G0)]
    type = PointValueSampler
    variable = 'sflux_g0'
    points = '0 0 143.552375 0 0 163.79425 0 0 183.79425'
    sort_by = x
    execute_on = timestep_end
  [../]
  [./Scalar Flux (G1)]
    type = PointValueSampler
    variable = 'sflux_g1'
    points = '0 0 143.552375 0 0 163.79425 0 0 183.79425'
    sort_by = x
    execute_on = timestep_end
  [../]
  [./Scalar Flux (G2)]
    type = PointValueSampler
    variable = 'sflux_g2'
    points = '0 0 143.552375 0 0 163.79425 0 0 183.79425'
    sort_by = x
    execute_on = timestep_end
  [../]
  [./Scalar Flux (G3)]
    type = PointValueSampler
    variable = 'sflux_g3'
    points = '0 0 143.552375 0 0 163.79425 0 0 183.79425'
    sort_by = x
    execute_on = timestep_end
  [../]
  [./Scalar Flux (G4)]
    type = PointValueSampler
    variable = 'sflux_g4'
    points = '0 0 143.552375 0 0 163.79425 0 0 183.79425'
    sort_by = x
    execute_on = timestep_end
  [../]
  [./Scalar Flux (G5)]
    type = PointValueSampler
    variable = 'sflux_g5'
    points = '0 0 143.552375 0 0 163.79425 0 0 183.79425'
    sort_by = x
    execute_on = timestep_end
  [../]
  [./Scalar Flux (G6)]
    type = PointValueSampler
    variable = 'sflux_g6'
    points = '0 0 143.552375 0 0 163.79425 0 0 183.79425'
    sort_by = x
    execute_on = timestep_end
  [../]
  [./Scalar Flux (G7)]
    type = PointValueSampler
    variable = 'sflux_g7'
    points = '0 0 143.552375 0 0 163.79425 0 0 183.79425'
    sort_by = x
    execute_on = timestep_end
  [../]
  [./Scalar Flux (G8)]
    type = PointValueSampler
    variable = 'sflux_g8'
    points = '0 0 143.552375 0 0 163.79425 0 0 183.79425'
    sort_by = x
    execute_on = timestep_end
  [../]
  [./Scalar Flux (G9)]
    type = PointValueSampler
    variable = 'sflux_g9'
    points = '0 0 143.552375 0 0 163.79425 0 0 183.79425'
    sort_by = x
    execute_on = timestep_end
  [../]
  [./Scalar Flux (G10)]
    type = PointValueSampler
    variable = 'sflux_g10'
    points = '0 0 143.552375 0 0 163.79425 0 0 183.79425'
    sort_by = x
    execute_on = timestep_end
  [../]
[]

[Postprocessors]
  [./TotalPower]
    type = ElementIntegralVariablePostprocessor
    block = '11000'
    variable = PowerDensity
    execute_on = linear
  [../]
[]


[Outputs]
  file_base = A0
  exodus = true
  print_linear_residuals = true
  print_perf_log = true
  csv = true
[]
