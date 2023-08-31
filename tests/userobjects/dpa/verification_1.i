#
# Constant spectrum (one everywhere between 0 and 1e-6), constant cross section (one)
# K-P model, C-12 atoms, Ed = 30 eV
#
# correct answer 1.183 x 10^9
# computed: 1.183432e+09
#
Ed = 30

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

[Functions]
  [KP]
    type = ParsedFunction
    expression = 't / 2 / ${Ed}'
  []
[]

[UserObjects]
  [KPtest]
    type = FunctionDPAUserObject
    damage_reaction_types = 'elastic'
    damage_functions = 'KP'
    irradiation_time = 1
    Z = '6'
    A = '12'
    number_densities = '1'

    scalar_flux = '1e6'
    energy_group_boundaries = '1e6 0'
    cross_section = '1'
  []
[]

[Postprocessors]
  [dpa]
    type = DPAPostprocessor
    dpa_object = KPtest
  []
[]

[Executioner]
  type = Steady
[]

[Outputs]
  csv = true
[]
