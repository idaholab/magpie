[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 60
  ny = 60
  xmin = -100
  xmax = 100
  ymin = 0
  ymax = 200
  elem_type = TRI3
[]

[AuxVariables]
  [./c]
    initial_condition = 1.0
  [../]
[]

[Variables]
  [./int]
  [../]
  [./vac]
  [../]
[]

[Kernels]
  [./dt_int]
    type = TimeDerivative
    variable = int
  [../]
  [./diff_int]
    type = Diffusion
    variable = int
  [../]
  [./dt_vac]
    type = TimeDerivative
    variable = vac
  [../]
  [./diff_vac]
    type = Diffusion
    variable = vac
  [../]
[]

[DiracKernels]
  [./int]
    variable = int
    type = MyTRIMDiracSource
    runner = runner
    ivar = 0
    defect = INT
  [../]
  [./vac]
    variable = vac
    type = MyTRIMDiracSource
    runner = runner
    ivar = 0
    defect = VAC
  [../]
[]

#[BCs]
#  [./Periodic]
#    [./all]
#      auto_direction = 'x y'
#    [../]
#  [../]
#[]

[UserObjects]
  [./thermal_fission]
    type = PKAFissionFragmentEmpirical
    relative_density = 1
    fission_rate = 0.001
  [../]
  [./rasterizer]
    type = MyTRIMRasterizer
    var = c
    periodic_var = int
    M = 40
    Z = 20
    site_volume = 0.0404 # nm^3 per UO2 unit
    pka_generator = thermal_fission
  [../]
  [./runner]
    type = MyTRIMDiracRun
    rasterizer = rasterizer
  [../]
[]

[Postprocessors]
  [./int]
    type = ElementIntegralVariablePostprocessor
    variable = int
    execute_on = timestep_end
  [../]
  [./vac]
    type = ElementIntegralVariablePostprocessor
    variable = vac
    execute_on = timestep_end
  [../]
[]

[Executioner]
  type = Transient
  num_steps = 1
  nl_abs_tol = 1e-10
[]

[Outputs]
  exodus = true
  csv = true
[]
