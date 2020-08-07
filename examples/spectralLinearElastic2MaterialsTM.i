[GlobalParams]
  displacements = 'u_x u_y u_z'
[]

[Mesh]
  [generated_mesh]
    type = GeneratedMeshGenerator
    dim = 3
    xmin = -5
    xmax = 5
    ymin = -5
    ymax = 5
    zmin = -5
    zmax = 5
    nx = 32
    ny = 32
    nz = 32
  []

  [cnode]
    type = ExtraNodesetGenerator
    coord = '0 0 0'
    new_boundary = 100
    input = generated_mesh
  []
[]

[AuxVariables]
  [./disp_x]
  [../]
  [./disp_y]
  [../]
  [./disp_z]
  [../]

[./stress_xy]
  order = CONSTANT
  family = MONOMIAL
[../]
[./elastic_strain_xy]
  order = CONSTANT
  family = MONOMIAL
[../]

[]

[Variables]
  [./u_x]
  [../]
  [./u_y]
  [../]
  [./u_z]
  [../]
  [./global_strain]
    order = SIXTH
    family = SCALAR
  [../]
[]

[Kernels]
  [./TensorMechanics]
  [../]
[]

[Functions]
  [./young_modulus_function]
    type = ParsedFunction
    vars = 'radius ym_min ym_max'
    vals = '2.0      1.0  4.0'
    value = 'if(sqrt(x*x+y*y+z*z)<radius,ym_min,ym_max)'
  [../]
[]

[ScalarKernels]
  [./global_strain]
    type = GlobalStrain
    variable = global_strain
    global_strain_uo = global_strain_uo
  [../]
[]

[Materials]
  [./strain]
    type = ComputeSmallStrain
    global_strain = global_strain
  [../]
  [./elasticity_tensor]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 1.0e4
    poissons_ratio = 0.3
    elasticity_tensor_prefactor = young_modulus_function
    # Is the line above during the righ thing by scaling the young modulus?
  [../]
  [./stress]
    type = ComputeLinearElasticStress
  [../]
  [./global_strain]
    type = ComputeGlobalStrain
    scalar_global_strain = global_strain
    global_strain_uo = global_strain_uo
  [../]
[]

[UserObjects]
  [./global_strain_uo]
    type = GlobalStrainUserObject
    applied_stress_tensor = '0 0 0 0 0 2.959e2'
    execute_on = 'Initial Linear Nonlinear'
  [../]
[]

[AuxKernels]
# global strain
[./disp_x]
    type = GlobalDisplacementAux
    variable = disp_x
    scalar_global_strain = global_strain
    global_strain_uo = global_strain_uo
    component = 1
  [../]
  [./disp_y]
    type = GlobalDisplacementAux
    variable = disp_y
    scalar_global_strain = global_strain
    global_strain_uo = global_strain_uo
    component = 1
  [../]
  [./disp_z]
    type = GlobalDisplacementAux
    variable = disp_z
    scalar_global_strain = global_strain
    global_strain_uo = global_strain_uo
    component = 2
  [../]

    [./stress_xy]
      type = RankTwoAux
      rank_two_tensor = stress
      variable = stress_xy
      index_i = 0
      index_j = 1
    [../]
    [./elastic_strain_xy]
        type = RankTwoAux
        rank_two_tensor = elastic_strain
        variable = elastic_strain_xy
        index_i = 0
        index_j = 1
    [../]
[]

[BCs]
  [./Periodic]
    [./all]
      variable = ' u_x u_y u_z'
      auto_direction = 'x y z'
    [../]
  [../]

  # fix center point location
  [./centerfix_x]
    type = DirichletBC
    boundary = 100
    variable = u_x
    value = 0
  [../]
  [./centerfix_y]
    type = DirichletBC
    boundary = 100
    variable = u_y
    value = 0
  [../]
  [./centerfix_z]
    type = DirichletBC
    boundary = 100
    variable = u_z
    value = 0
  [../]
[]

[Preconditioning]
  [./SMP]
    type = SMP
    full = true
  [../]
[]

[Executioner]
  type = Transient
  solve_type = 'PJFNK'

  line_search = basic

  petsc_options_iname = '-pc_type -ksp_gmres_restart -sub_ksp_type -sub_pc_type -pc_asm_overlap'
  petsc_options_value = 'asm         31   preonly   lu      1'

  l_max_its = 30
  nl_max_its = 12

  l_tol = 1.0e-4

  nl_rel_tol = 1.0e-09
  nl_abs_tol = 1.0e-09

  start_time = 0.0
  num_steps = 2
[]

[VectorPostprocessors]
  [./linevaluex]
    type = LineValueSampler
    variable = 'elastic_strain_xy stress_xy'
    start_point = '-5 0 0'
    end_point = '5 0 0'
    num_points = 101
    sort_by = x
    execute_on = final
  [../]
  [./linevaluey]
    type = LineValueSampler
    variable = 'elastic_strain_xy stress_xy'
    start_point = '0 -5 0'
    end_point = '0 5 0'
    num_points = 101
    sort_by = y
    execute_on = final
  [../]
  [./linevaluez]
    type = LineValueSampler
    variable = 'elastic_strain_xy stress_xy'
    start_point = '0 0 -5'
    end_point = '0 0 5'
    num_points = 101
    sort_by = z
    execute_on = final
  [../]
[]

[Outputs]
  exodus = true
  execute_on = 'INITIAL FINAL'
  perf_graph = true
  [./comp]
    type = CSV
  [../]
[]
