[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 25
  ny = 25
  nz = 25
  xmin = -250
  xmax = 250
  ymin = -250
  ymax = 250
  zmin = -250
  zmax = 250
[]

[Problem]
  kernel_coverage_check = false
[]

[Variables]
  [dummy]
  []
[]

[AuxVariables]
  [f]
    order = CONSTANT
    family = MONOMIAL
    [InitialCondition]
      type = FunctionIC
      function = setf
    []
  []
[]

[Functions]
  [setf]
    type = ParsedFunction
    expression = 'r := sqrt(x*x + y*y + z*z); exp(-r / 10)'
  []
[]

[Postprocessors]
  [rms]
    type = RMSDistance
    variable = f
    point = '0 0 0'
  []

  [intf]
    type = ElementIntegralVariablePostprocessor
    variable = f
  []
[]

[Executioner]
  type = Steady
[]

[Outputs]
  csv = true
[]
