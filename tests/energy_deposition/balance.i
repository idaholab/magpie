#
# This input logs the total energy of all PKAs each timestep and compares
# it to the deposited electronic and nuclear energy (due to stopping) and
# the potential energy of the vacancies that are created.
# Epka - Estopping - Evac != 0 to ensure conservation of energy!
# The equation above only holds under periodic boundary conditions or if
# no recoil leaves the simulation box!
#

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 50
  ny = 50
  xmin = 0
  xmax = 100
  ymin = 0
  ymax = 100
[]

[Functions]
  [./ebalance]
    type = ParsedFunction
    vars = 'Edep Evac Epka'
    vals = 'Edep Evac Epka'
    value = 'Epka - Edep - Evac'
  [../]
[]

[Variables]
  [./c]
    initial_condition = 1.0
  [../]
[]

[BCs]
  # the composition field used by the rasterizer is marked as periodic
  [./Periodic]
    [./all]
      variable = c
      auto_direction = 'x y'
    [../]
  [../]
[]

[AuxVariables]
  [./edep]
    # deposited energy density (electronic and nuclear stopping)
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./vac]
    # number density of vacancies created
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./evac]
    # potential energy density of the vacancies created
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[AuxKernels]
  [./edep]
    type = MyTRIMElementEnergyAux
    variable = edep
    runner = runner
    execute_on = timestep_end
  [../]
  [./int]
    type = MyTRIMElementResultAux
    variable = vac
    runner = runner
    defect = VAC
    ivar = 0
    execute_on = timestep_end
  [../]
  [./evac]
    type = ParsedAux
    variable = evac
    args = vac
    # binding energy is 3eV, so that is what we take as vacancy formation energy
    function = 3.0*vac
    execute_on = timestep_end
  [../]
[]

[UserObjects]
  [./constant]
    # fixed energy (10keV) model ions (approx Calcium)
    type = PKAConstant
    E = 10000
    m = 40
    Z = 20
    pka_rate = 1e-3
  [../]
  [./rasterizer]
    type = MyTRIMRasterizer
    var = c
    M = 40
    Z = 20
    site_volume = 0.0404 # nm^3 per atom
    pka_generator = constant
    trim_module = ENERGY_DEPOSITION
  [../]
  [./runner]
    type = MyTRIMElementRun
    rasterizer = rasterizer
  [../]
[]

[Postprocessors]
  [./Edep]
    type = ElementIntegralVariablePostprocessor
    variable = edep
    execute_on = timestep_end
  [../]
  [./Evac]
    type = ElementIntegralVariablePostprocessor
    variable = evac
    execute_on = timestep_end
  [../]
  [./Epka]
    type = MyTRIMPKAInfo
    rasterizer = rasterizer
    value_type = TOTAL_ENERGY
  [../]

  [./Ebalance]
    type = FunctionValuePostprocessor
    function = ebalance
  [../]
[]

[Problem]
  kernel_coverage_check = false
  solve = true
[]

[Executioner]
  type = Transient
  num_steps = 3
  nl_abs_tol = 1e-10
[]

[Outputs]
  csv = true
[]
