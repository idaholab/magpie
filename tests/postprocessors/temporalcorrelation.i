[Mesh]
  type = MyTRIMMesh
  dim = 2
  nx = 10
  ny = 10
  xmax = 6
  ymax = 9
[]

[Variables]
  [./c1]
  [../]
  [./c2]
    [./InitialCondition]
      type = FunctionIC
      function = x
    [../]
  [../]
[]

[Kernels]
  [./dt1]
    type = TimeDerivative
    variable = c1
  [../]
  [./dt2]
    type = TimeDerivative
    variable = c2
  [../]
  [./s1]
    type = BodyForce
    variable = c1
    value = 1.23
  [../]
  [./s2]
    type = Diffusion
    variable = c2
  [../]
[]

[Executioner]
  type = Transient
  num_steps = 5
[]


[Postprocessors]
  [./C1]
    type = TemporalCorrelation
    variable = c1
  [../]
  [./C2]
    type = TemporalCorrelation
    variable = c2
  [../]
[]

[Outputs]
  csv = true
[]
