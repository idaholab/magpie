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
    cross_section_output_filename = erxs_H_simple_out.csv
    mu_L_output_filename = erxs_H_simple_mu_L_out.csv
    atomic_mass = 1
    legendre_order = 7
    neutron_energy_limits = '1e7 1e6 1e5 1e4 1e3 1e2 1e1 1e0'
    recoil_energy_limits = '2097152 1048576 524288 262144 131072 65536 32768 16384 8192 4096 2048 1024 512 256 128 64 32 16 8 4 2 1 0'

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
