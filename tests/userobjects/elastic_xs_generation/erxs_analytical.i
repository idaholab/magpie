
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
  [./test_var]
  [../]
[]

[UserObjects]
  [./XSGenerator]
     type = ElasticRecoilCrossSectionUserObject
     erxs_output_file_name = erxs_analytical_out.csv
     mu_L_output_file_name = erxs_analytical_mu_L_out.csv
     atomic_mass = 12
     legendre_order = 2
     neutron_energy_limits = '2 1'
     recoil_energy_limits = '1.0 0.9 0.8 0.7 0.6 0.5 0.4 0.3'
     neutron_spectrum = neutron_spectrum
     scattering_law = scattering_law
     elastic_xs = elastic_xs
     execute_on = timestep_end
  [../]
[]

# t is equal to Ei
[Functions]
  [./neutron_spectrum]
     type = ParsedFunction
     value = '1'
  [../]

  # t is equal to Ei
  [./elastic_xs]
    type = ParsedFunction
    value = '1'
  [../]

  # t is equal to mu_c
  [./scattering_law]
      type = ParsedFunction
     value = '0.5'
  [../]
[]

[Postprocessors]
  [./test]
    type = AverageNodalVariableValue
    variable = test_var
  [../]
[]

[Executioner]
  type = Steady
[]
