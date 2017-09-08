[Mesh]
  type = CartesianMesh
  dim = 2
  dx = '0.15 0.95 0.15'
  dy = '0.15 0.95 0.15'

  subdomain_id = '2 2 2
                 2 1 2
                 2 2 2'
[]

[TransportSystems]
  particle = neutron
  equation_type = eigenvalue
  G = 2
  ReflectingBoundary = 'right top left bottom'

  [./sn]
    scheme = SAAF-CFEM-SN
    family = LAGRANGE
    order = FIRST
    AQorder = 4
    AQtype = Level-Symmetric
    hide_angular_flux = true
  [../]
[]

[Materials]
  [./fuel]
    type = ConstantNeutronicsMaterial
    block = 1
    ngroup = 2
    sigma_t = '0.2222   0.8333'
    L = 0
    sigma_s = '
             0.1921 0.00001
             0.020  0.7483'
    fissile = true
    nu_sigma_f = '0.000001   0.1350'
    chi = '0.999999 0.000001'
    fromFile = false
  [../]

  [./reflector]
    type = ConstantNeutronicsMaterial
    block = 2
    ngroup = 2
    sigma_t = '0.1666 1.1111'
    L = 0
    sigma_s = '
             0.1265 0.00001
             0.04  1.1010'
    fromFile = false
  [../]

[]

[AuxVariables]
  [./N92235]
    order = CONSTANT
    family = MONOMIAL
    initial_condition = 0.02
  [../]

  [./N92238]
    order = CONSTANT
    family = MONOMIAL
    initial_condition = 0.32
  [../]
[]

[UserObjects]
  [./pka]
    type = RadiationDamageFission
    points = '0.625 0.625 0.0'
    energy_group_boundaries = '20e6 1 0'
    target_isotope_names = 'U235 U238'
    number_densities = 'N92235 N92238'
    scalar_fluxes = 'flux_moment_g0_L0_M0 flux_moment_g1_L0_M0'
    # each row are the recoil cross sections for a single isotope && all group (g=0,1)
    fission_cross_sections = '
                              0.0 1.0
                              0.01 0.0
                            '
  [../]
[]

#[MultiApps]
#  [./magpie]
#    type = TransientMultiApp
#    app_type = MagpieApp
#    execute_on = timestep_end
#  [../]
#[]

[Executioner]
  type = NonlinearEigen

  free_power_iterations = 1
  bx_norm = 'fission_source_integral'
  source_abs_tol = 1e-8

  #Preconditioned JFNK (default)
  solve_type = 'PJFNK'
  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart '
  petsc_options_value = 'hypre boomeramg 50'
[]

[Outputs]
  file_base = pka_distribution_fission
  exodus = true
[]
