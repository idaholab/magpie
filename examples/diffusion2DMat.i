
[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 20
  ny = 20
  xmax = 100
  ymax = 100
[]

[Variables]
  [./diffused]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = SmoothCircleIC
      x1 = 50
      y1 = 50
      radius = 10
      int_width = 1
      invalue = 1
      outvalue = 0.0001
    [../]
  [../]
[]

[Kernels]
  [./diff]
    type = MatDiffusion
    variable = diffused
    diffusivity = diffusivity
  [../]
  [./dt]
    type = CoefTimeDerivative
    variable = diffused
    Coefficient = 1
  [../]
[]

[BCs]
  [./Periodic]
    #Note: Enable either "auto" or both "manual" conditions for this example
    # active = 'manual_x manual_y'

    # Can use auto_direction with Generated Meshes
    [./auto]
      variable = diffused
      auto_direction = 'x y'
    [../]
  [../]
[]

[Materials]
  [./k]
    type = GenericConstantMaterial
    prop_names = 'diffusivity'
    prop_values = '400'
  [../]
[]

[Executioner]
  type = Transient
  dt = 0.005
  num_steps = 21
  nl_rel_tol = 1e-10
  nl_abs_tol = 1.0e-10
[]

[VectorPostprocessors]
  [./linevalue]
    type = LineValueSampler
    variable = 'diffused'
    start_point = '0 50 0'
    end_point = '100 50 0'
    num_points = 101
    sort_by = x
  [../]
[]

[Outputs]
  execute_on = 'timestep_end'
  exodus = true
  perf_graph = true
  [./something]
    type = CSV
    show = 'linevalue'
  [../]
[]
