#
# Example problem showing how to use the DerivativeParsedMaterial with SplitCHParsed.
# The free energy is identical to that from SplitCHMath, f_bulk = 1/4*(1-c)^2*(1+c)^2.
#

[Mesh]
  type = GeneratedMesh
  dim = ${dim}
  nx = 20
  ny = 5
  nz = 5
  xmin = -25
  ymin = -10
  zmin = -10
  xmax = 25
  ymax = 10
  zmax = 10
[]

[Variables]
  [./c]
  [../]
[]

[ICs]
  [./xstep]
    type = FunctionIC
    variable = c
    function = 'if(x>0,
      if(x<1, 1+2*x^3-3*x^2, 0),
      if(x<-24, 3*(x+25)^2-2*(x+25)^3, 1)
    )'
  [../]
[]

[Kernels]
  [./c_dot]
    type = TimeDerivative
    variable = c
  [../]
  [./green]
    type = RadialGreensSource
    variable = c
    gamma = 1
    convolution = green
  [../]
[]

[UserObjects]
  [./green]
    type = RadialGreensConvolution
    v = c
    r_cut = 10
    normalize = true
    function = 'exp(-x/5)'
  [../]
[]

[BCs]
  [./Periodic]
    [./all]
      auto_direction = ${pbc}
    [../]
  [../]
[]

[VectorPostprocessors]
  [./c]
    type = LineValueSampler
    execute_on = 'INITIAL TIMESTEP_END'
    variable = c
    start_point = '-25 0 0'
    end_point = '25 0 0'
    num_points = 40
    sort_by = x
  [../]
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  scheme = bdf2

  petsc_options_iname = '-pc_type -sub_pc_type'
  petsc_options_value = 'asm      ilu         '

  l_max_its = 60
  l_tol = 1e-4
  nl_max_its = 30
  nl_rel_tol = 1e-9
  nl_abs_tol = 1e-10

  dt = 1
  num_steps = 1
[]

[Outputs]
  execute_on = FINAL
  file_base = dimension${dim}_out
  csv = true
  print_linear_residuals = false
[]
