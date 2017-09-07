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

  [./N08016]
    order = CONSTANT
    family = MONOMIAL
    initial_condition = 0.66
  [../]
[]

[UserObjects]
  [./pka]
    type = RadiationDamageSN
    L = 1
    points = '0.625 0.625 0.0'
    energy_group_boundaries = '20e6 1 0'
    target_isotope_names = 'U235 U238 O16'
    recoil_isotope_names = 'U235 U238 O16'
    number_densities = 'N92235 N92238 N08016'
    angular_variables = 'aflux_g0_d0 aflux_g0_d1 aflux_g0_d2 aflux_g0_d3 aflux_g0_d4 aflux_g0_d5
                         aflux_g0_d6 aflux_g0_d7 aflux_g0_d8 aflux_g0_d9 aflux_g0_d10 aflux_g0_d11
                         aflux_g1_d0 aflux_g1_d1 aflux_g1_d2 aflux_g1_d3 aflux_g1_d4 aflux_g1_d5
                         aflux_g1_d6 aflux_g1_d7 aflux_g1_d8 aflux_g1_d9 aflux_g1_d10 aflux_g1_d11'
    # each row are the recoil cross sections for a single isotope && l (l=0,1)
    # order: l=0 1->1 2->1 1->2 2->2
    #        l=1 1->1 2->1 1->2 2->2
    # this order is consistent with sigma_s in ConstantNeutronicsMaterial
    recoil_cross_sections = '
                              0.1 0.0 0.5 0.0
                              0.0 0.0 0.0 0.0
                              0.09 0.0 0.45 0.0
                              0.0 0.0 0.0 0.0
                              0.4 0.0 1.2 0.2
                              0.01 0.0 0.1 0.05
                            '
    aqdata = sn_aqdata
  [../]
[]

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
  file_base = pka_distribution_sn
  exodus = true
[]
