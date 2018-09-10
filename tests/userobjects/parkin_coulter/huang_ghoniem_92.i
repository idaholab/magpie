#
# This test follows the Huang & Ghoniem "Neutron displacement damage cross sections for SiC"
# J. of Nucl. Mat. 199, (1993), 221-230
# Results (with extended energy range) from this test case are included in the writeup
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

[UserObjects]
  [./parkin_coulter]
     type = PolyatomicRecoil

     Z = '6 14'
     A = '12 28'
     number_fraction = '0.5 0.5'
     displacement_thresholds = '16.3 92.6'

     logarithmic_energy_spacing = 1.1

     Emax = 1.0e3
     execute_on = timestep_end

     displacement_file_base = 'huang_ghoniem_n_ij'
  [../]
[]

[Executioner]
  type = Steady
[]
