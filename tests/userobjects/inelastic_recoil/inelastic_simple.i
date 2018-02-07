
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
     type = InelasticRecoil
     cross_section_output_filename = inelastic_simple_out.csv
     mu_L_output_filename = inelastic_simple_mu_L_out.csv
     atomic_mass = 1
     legendre_order = 2
     Q = '0.4'
     neutron_energy_limits = '2 1'
     recoil_energy_limits = '1.2 1.1 1.0 0.9 0.8 0.7'
     neutron_spectrum = neutron_spectrum
     scattering_xs = scattering_xs
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
  [./scattering_xs]
    type = ParsedFunction
    value = '1'
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
