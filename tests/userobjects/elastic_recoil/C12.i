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
    cross_section_output_filename = C12_out.csv
    mu_L_output_filename = C12_mu_L_out.csv
    atomic_mass = 12
    legendre_order = 4
    neutron_energy_limits = '203240 145000'
    recoil_energy_limits = '1e4 1000 0'
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
    expression = '1/t'
  []

  # t is equal to Ei
  [scattering_xs]
    type = PiecewiseLinear
    data_file = C12_xs.csv
    xy_in_file_only = false
    format = columns
  []

  # t is equal to mu_c
  [scattering_law]
    type = ParsedFunction
    expression = '0.5+0.1*t'
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
