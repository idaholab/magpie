#
# This test follows the Huang & Ghoniem "Neutron displacement damage cross sections for SiC"
# J. of Nucl. Mat. 199, (1993), 221-230
#
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

[Functions]
  [./HG_SiSi]
    type = PiecewiseLinear
    data_file = 'HG_SiSi.csv'
    format = columns
  [../]

  [./HG_SiC]
    type = PiecewiseLinear
    data_file = 'HG_SiC.csv'
    format = columns
  [../]

  [./HG_CSi]
    type = PiecewiseLinear
    data_file = 'HG_CSi.csv'
    format = columns
  [../]

  [./HG_CC]
    type = PiecewiseLinear
    data_file = 'HG_CC.csv'
    format = columns
  [../]
[]

[UserObjects]
  [./parkin_coulter]
     type = ParkinCoulterDPAUserObject
     damage_reaction_types = 'elastic inelastic'
     irradiation_time = 1

     Z = '6 14'
     A = '12 28'
     number_densities = '0.5 0.5'
     scalar_flux = '0.048138819 0.16809268 0.142519665 0.122144812 0.104682829
                    0.08971709 0.04904444 0.04213553 0.233524136'
     energy_group_boundaries = '3.455e6 1e6 1e5 1e4 1e3 1e2 10.0 2.4 0.625 1e-5'
     cross_section = '2.136 3.7 4.622 4.737 4.748 4.75 4.75 4.757 13.198;
                      2.938 4.205 1.980 1.918 1.952 1.956 1.957 1.959 4.034;
                      0     0     0     0     0     0     0     0     0;
                      0.10848 0   0     0     0     0     0     0     0'
     Q = '0 0; 0 -1.779e+6'

     displacement_thresholds = '16.3 92.6'
     logarithmic_energy_spacing = 1.25
  [../]

  [./huang_ghoniem_CSV]
     type = FunctionDPAUserObject
     damage_reaction_types = 'elastic inelastic'
     damage_functions = 'HG_CC HG_CSi; HG_SiC HG_SiSi'
     irradiation_time = 1

     Z = '6 14'
     A = '12 28'
     number_densities = '0.5 0.5'
     scalar_flux = '0.048138819 0.16809268 0.142519665 0.122144812 0.104682829
                    0.08971709 0.04904444 0.04213553 0.233524136'
     energy_group_boundaries = '3.455e6 1e6 1e5 1e4 1e3 1e2 10.0 2.4 0.625 1e-5'
     cross_section = '2.136 3.7 4.622 4.737 4.748 4.75 4.75 4.757 13.198;
                      2.938 4.205 1.980 1.918 1.952 1.956 1.957 1.959 4.034;
                      0     0     0     0     0     0     0     0     0;
                      0.10848 0   0     0     0     0     0     0     0'
     Q = '0 0; 0 -1.779e+6'
  [../]
[]

[Postprocessors]
  [./dpa_from_function]
    type = DPAPostprocessor
    dpa_object = huang_ghoniem_CSV
  [../]

  [./dpa_from_parkin_coulter]
    type = DPAPostprocessor
    dpa_object = parkin_coulter
  [../]
[]

[Executioner]
  type = Steady
[]

[Outputs]
  csv = true
[]
