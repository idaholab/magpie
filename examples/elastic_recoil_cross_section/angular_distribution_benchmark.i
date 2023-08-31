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
  [test_var]
  []
[]

[UserObjects]
  [XSGenerator]
    type = ElasticRecoil

    # Inputs
    cross_section_output_filename = angular_distribution_out.csv
    mu_L_output_filename = angular_distribution_mu_L_out.csv
    atomic_mass = 2
    legendre_order = 20
    quadrature_order = 800
    neutron_energy_limits = '1.2 1 0.8'
    recoil_energy_limits = '1 0.95 0.9 0.85 0.8 0.75 0.7 0.65 0.6 0.55
                             0.5 0.45 0.4 0.35 0.3 0.25 0.2 0.15 0.1 0.05 0'

    # Functions
    neutron_spectrum = neutron_spectrum
    scattering_law = scattering_law
    scattering_xs = scattering_xs

    execute_on = timestep_end
  []
[]

# t is equal to Ei
[Functions]
  [neutron_spectrum]
    type = ParsedFunction
    expression = '1'
  []

  # t is equal to Ei
  [scattering_xs]
    type = ParsedFunction
    expression = '1'
  []

  # t is equal to mu_c
  [scattering_law]
    type = ParsedFunction
    expression = '0.5'
  []
[]

[Postprocessors]
  [test]
    type = AverageNodalVariableValue
    variable = test_var
  []
[]

[Executioner]
  type = Steady
[]
