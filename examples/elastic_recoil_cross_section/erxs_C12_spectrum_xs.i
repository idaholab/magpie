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
    cross_section_output_filename = erxs_C_spectrum_xs_out.csv
    mu_L_output_filename = erxs_C_spectrum_xs_mu_L_out.csv
    atomic_mass = 12
    legendre_order = 7
    # 14 group CASMO lib
    neutron_energy_limits = '10000000  2231000  821000  5530  48.052  4  1.15  0.972  0.625  0.35  0.28  0.14  0.058  0.03  0.00001'
    # 25 group CASMO lib
    recoil_energy_limits = '10000000  6065500  3679000  2231000  1353000  821000  500000  111000  9118  5530  148.728  15.968  9.877
                             4  1.855  1.15  1.097  1.02  0.972  0.625  0.35  0.28  0.14  0.058  0.03  0.00001'
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
    data_file = C12_input.csv
    xy_in_file_only = false
    format = columns
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
