[Mesh]
  file = 3x3-Standard.e
  parallel_type = REPLICATED  # must be serial with MatID
  uniform_refine = 0
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
  [./Power]
    type = VectorReactionRate
    variable = PowerDensity
    #nusigf = kappa_sigma_fission
    cross_section = kappa_sigma_fission
    block = '11000'
  [../]
[]

[YAKXSLibraries]
  [./A0]
     type = BaseLibObject
     library_name = A0_A11000
     library_type = MultigroupLibrary
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
  petsc_options_iname =  '-pc_type -pc_hypre_type -pc_hypre_boomeramg_tol -pc_hypre_boomeramg_max_iter -ksp_gmres_restart'
  petsc_options_value =  'hypre boomeramg 1e-4 10 50'
  line_search = none

  l_tol = 1e-4
  free_power_iterations = 4
  source_abs_tol = 1e-8
  normal_factor = 100.0
  normalization = TotalPower
[]

[VectorPostprocessors]
  active = ''
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
