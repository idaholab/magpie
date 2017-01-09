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
             1 1 2 2 1 1 1 1
             1 1 2 2 1 1 1 1
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
  evaluate_on = element
[]

[TransportSystems]
  particle = neutron
  equation_type = eigenvalue
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
  [../]
[]

[AuxVariables]
  [./CR_Boron]
    order = CONSTANT
    family = MONOMIAL
    [./InitialCondition]
      type = ConstantIC
      value = 0.0
    [../]
  [../]

  [./Tfuel]               # Units [=] Kelvin
    order = CONSTANT
    family = MONOMIAL
    [./InitialCondition]
      type = ConstantIC
      value = 300.0
    [../]
  [../]

  [./PowerDensity]
    order = CONSTANT
    family = MONOMIAL
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
    type = CoupledFeedbackNeutronicsMaterial
    block = '1'
    material_id = 104
  [../]

  [./fuel]
    type = CoupledFeedbackNeutronicsMaterial
    block ='2'
    material_id = 106
  [../]
[]

[Postprocessors]
  [./UnscaledTotalPower]
    type = FluxRxnIntegral
    block = '2'
    coupled_flux_groups = 'sflux_g0 sflux_g1 sflux_g2 sflux_g3 sflux_g4 sflux_g5 sflux_g6 sflux_g7 sflux_g8 sflux_g9 sflux_g10'
    cross_section  = kappa_sigma_fission
    execute_on = 'initial linear'
  [../]

  [./PowerScaling]
    type = PowerModulateFactor
    power_pp = UnscaledTotalPower
    rated_power = 1000000
    power_modulating = 1
    execute_on = 'initial linear'
    outputs = all
  [../]
[]

[Executioner]
  type = NonlinearEigen
  solve_type = 'PJFNK'
  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart -pc_hypre_boomeramg_max_iter -pc_hypre_boomeramg_tol'
  petsc_options_value = 'hypre boomeramg 100 20 1.0e-6'
  free_power_iterations = 2
  l_max_its = 100
  pfactor = 1.0e-3
  source_abs_tol = 1e-8
  line_search = 'none'
[]

[Outputs]
  print_linear_residuals = true
  print_perf_log = true
  exodus = true
[]
