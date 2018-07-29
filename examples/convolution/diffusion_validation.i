#
# Validate the Green's Function convolution user object by comparing the
# time integration of the Diffusion equation using the TimeDerivative and Diffusion
# Kernels (c1) to the repeated application of the Green's Function of the Diffusion
# equation (c2).
#
# Compare the results in gnuplot using
#   gnuplot> set datafile separator ','; set key autotitle columnhead
#   gnuplot> pl 'diffusion_validation_out_c_0020.csv' u 1:5 w l, 'diffusion_validation_out_c_0020.csv' u 1:6 w l
#

[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 200
  xmin = -10
  xmax = 10
[]

[Functions]
  [./ic]
    type = ParsedFunction
    # box profile initial condition (Note: this causes Gibbs oscillations that get amplified in the c2 solution)
    value = 'if(abs(x)<=2,1,0)'
  [../]
  [./G]
    type = ParsedFunction
    # adjust the simulation timestep here!
    vars = dt
    vals = 0.1
    # the diffusion Green's function for D=1 is
    # exp(-r^2/(4*dt))
    # Note1: We need to multiply with r^2 to compensate for the r^-2 geometrical
    # attenuation term that is multiplied on by the userobject (x <=> r)
    # Note2: we use normalize = true below to autonormalize the Green's function
    # so we can leave off the prefactor
    value = exp(-x^2/(4*dt))*x^2
  [../]
[]

[Variables]
  [./c1]
    [./InitialCondition]
      type = FunctionIC
      function = ic
    [../]
  [../]
  [./c2]
    [./InitialCondition]
      type = FunctionIC
      function = ic
    [../]
  [../]
[]

[BCs]
  [./Periodic]
    [./all]
      auto_direction = x
    [../]
  [../]
[]

[Kernels]
  [./dt1]
    type = TimeDerivative
    variable = c1
  [../]
  [./diff1]
    type = Diffusion
    variable = c1
  [../]

  [./dt2]
    type = TimeDerivative
    variable = c2
  [../]
  [./source2]
    type = RadialGreensSource
    variable = c2
    # set gamma = 1/dt
    # otherwise the convolution result is not fully applied. We are not using the
    # convolution result as a rate, but want to change the variable field into the
    # convolution result
    gamma = 10
    convolution = convolution
  [../]
[]

[UserObjects]
  [./convolution]
    type = RadialGreensConvolution
    normalize = true
    v = c2
    function = G
    r_cut = 2.0
  [../]
[]

[VectorPostprocessors]
  [./c]
    type = LineValueSampler
    execute_on = 'INITIAL TIMESTEP_END'
    variable = 'c1 c2'
    start_point = '-10 0 0'
    end_point = '10 0 0'
    # sample at element edges (nodes) and centers |---*---|---*---|
    num_points = 399
    sort_by = x
  [../]
[]

[Executioner]
  type = Transient
  # use a small timestep to reduce time integration error
  dt = 0.1
  num_steps = 20
[]

[Outputs]
  csv = true
[]
