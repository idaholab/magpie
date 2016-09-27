[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 10
  ny = 10
  xmin = 0
  xmax = 1
  ymin = 0
  ymax = 1
[]

[Variables]
  [./u]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[Kernels]
  active = 'diff'

  [./diff]
    type = Diffusion
    variable = u
  [../]
[]

[BCs]

  [./left]
    type = DirichletBC
    variable = u
    boundary = '1 2 3'
    value = 0
  [../]

  [./right]
    type = DirichletBC
    variable = u
    boundary = 0
    value = 1
  [../]

[]

[AuxVariables]
  [./scalar_flux_g0]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = FunctionIC
      function = set_phi_g0
    [../]
  [../]

  [./scalar_flux_g1]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = FunctionIC
      function = set_phi_g1
    [../]
  [../]

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

[Functions]
  [./set_phi_g0]
    type = ParsedFunction
    value = 'sin(pi*x)*sin(pi*y)'
  [../]
  [./set_phi_g1]
    type = ParsedFunction
    value = '0.5+0.05*(x*x+y*y)'
  [../]
[]

[UserObjects]
  [./fission_damage_sampler]
    type = NeutronicsSpectrumSamplerFission
    points = '0.5 0.5 0.0'
    energy_group_boundaries = '0.0 0.5e-6 0.75'
    target_isotope_names = 'U235 U238'
    number_densities = 'N92235 N92238'
    scalar_fluxes = 'scalar_flux_g1 scalar_flux_g0'
    # each row are the recoil cross sections for a single isotope && all group (g=0,1)
    fission_cross_sections = '1 0
                              0 0.01'
    execute_on = initial
  [../]
[]

[MultiApps]
  [./radiation_damage_app]
    type = FullSolveMultiApp
    input_files = damage_sub.i
    positions = '0.5 0.5 0.0'
    execute_on = nonlinear
  []
[]

[Transfers]
  [./radiation_damage_transfer]
    type = MultiAppNeutronicsSpectrumTransfer
    multi_app = radiation_damage_app
    pka_neutronics = neutronics_fission_generator
    radiation_damage_sampler = fission_damage_sampler
    direction = to_multiapp
  [../]
[]

[Executioner]
  type = Steady
  solve_type = 'NEWTON'
[]

[Outputs]
  file_base = coupled_fission_mockup
  exodus = true
[]
