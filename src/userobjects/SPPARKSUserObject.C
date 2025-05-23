/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#include "SPPARKSUserObject.h"

#include "MooseRandom.h"
#include "libmesh/mesh_tools.h"
#include <numeric>

registerMooseObject("MagpieApp", SPPARKSUserObject);

InputParameters
SPPARKSUserObject::validParams()
{
  InputParameters params = GeneralUserObject::validParams();
  params.addParam<std::string>("file", "", "SPPARKS input file");
  params.addParam<bool>("spparks_only", false, "Whether to run SPPARKS independently of MOOSE");

  params.addParam<std::vector<unsigned int>>(
      "from_ivar", {}, "Index into SPPARKS iarray.  This data will be extracted from SPPARKS.");
  params.addParam<std::vector<unsigned int>>(
      "from_dvar", {}, "Index into SPPARKS darray.  This data will be extracted from SPPARKS.");
  params.addParam<std::vector<unsigned int>>(
      "to_ivar", {}, "Index into SPPARKS iarray.  This data will be copied to SPPARKS.");
  params.addParam<std::vector<unsigned int>>(
      "to_dvar", {}, "Index into SPPARKS darray.  This data will be copied to SPPARKS.");

  params.addParam<std::vector<std::string>>("int_vars", {}, "Integer Vars to send to SPPARKS.");
  params.addParam<std::vector<std::string>>("double_vars", {}, "Double Vars to send to SPPARKS.");
  params.addParam<std::vector<std::string>>(
      "sol_vars", {}, "Double solving-for Vars obtained from SPPARKS.");

  params.addParam<Real>("xmin", 0.0, "Lower X Coordinate of the generated mesh");
  params.addParam<Real>("ymin", 0.0, "Lower Y Coordinate of the generated mesh");
  params.addParam<Real>("zmin", 0.0, "Lower Z Coordinate of the generated mesh");
  params.addParam<Real>("xmax", 1.0, "Upper X Coordinate of the generated mesh");
  params.addParam<Real>("ymax", 1.0, "Upper Y Coordinate of the generated mesh");
  params.addParam<Real>("zmax", 1.0, "Upper Z Coordinate of the generated mesh");

  params.addParam<bool>("init_spparks", false, "Set values of SPPARKS variables at time zero.");
  params.addParam<bool>(
      "one_time_run", false, "Run SPPARKS just once at time zero."); // added by YF
  params.addParam<Real>("time_spparks_time_ratio", 1, "Ratio of time to SPPARKS time.");

  // Hide from input file dump
  // params.addPrivateParam<std::string>("built_by_action", "");
  return params;
}

SPPARKSUserObject::SPPARKSUserObject(const InputParameters & params)
  : GeneralUserObject(params),
    _spparks(nullptr),
    _file(getParam<std::string>("file")),
    _spparks_only(getParam<bool>("spparks_only")),
    _from_ivar(getParam<std::vector<unsigned int>>("from_ivar")),
    _from_dvar(getParam<std::vector<unsigned int>>("from_dvar")),
    _to_ivar(getParam<std::vector<unsigned int>>("to_ivar")),
    _to_dvar(getParam<std::vector<unsigned int>>("to_dvar")),
    _xmin(getParam<Real>("xmin")),
    _ymin(getParam<Real>("ymin")),
    _zmin(getParam<Real>("zmin")),
    _xmax(getParam<Real>("xmax")),
    _ymax(getParam<Real>("ymax")),
    _zmax(getParam<Real>("zmax")),
    _init_spparks(getParam<bool>("init_spparks")),
    _time_ratio(getParam<Real>("time_spparks_time_ratio")),
    _initialized(false),
    _one_time_run(getParam<bool>("one_time_run")),
    _n_spparks_run(0),
    _int_vars(),
    _double_vars(),
    _sol_vars(),
    _last_time(std::numeric_limits<Real>::min())
{
  for (const auto & name : getParam<std::vector<std::string>>("int_vars"))
    _int_vars.push_back(&_fe_problem.getStandardVariable(0, name));
  if (_int_vars.size() != _to_ivar.size())
    mooseError("Mismatch with int_vars and to_ivar");

  for (const auto & name : getParam<std::vector<std::string>>("double_vars"))
    _double_vars.push_back(&_fe_problem.getStandardVariable(0, name));
  if (_double_vars.size() != _to_dvar.size())
    mooseError("Mismatch with double_vars and to_dvar");

  for (const auto & name : getParam<std::vector<std::string>>("sol_vars"))
    _sol_vars.push_back(&_fe_problem.getStandardVariable(0, name));

  _console << "\n>>>> STARTING SPPARKS <<<<\n";
  spparks_open(0, nullptr, _communicator.get(), &_spparks);
  if (!_spparks)
    mooseError("Error initializing SPPARKS");

  _console << "\n>>>> RUNNING SPPARKS FILE " << _file << " <<<<\n";
  spparks_file(_spparks, _file.c_str());

  int * iptr;
  double * dptr;

  // Extract and print information about the SPPARKS internals
  getSPPARKSPointer(iptr, "dimension");
  _dim = *iptr;
  _console << "\n>>>> SPPARKS DIMENSION: " << _dim << " <<<<\n";

  getSPPARKSPointer(dptr, "boxxlo");
  _console << "\n>>>> SPPARKS BOXXLO: " << *dptr << " <<<<\n";

  getSPPARKSPointer(iptr, "nlocal");
  _console << "\n>>>> SPPARKS NLOCAL: " << *iptr << " <<<<\n";
}

SPPARKSUserObject::~SPPARKSUserObject() { spparks_close(_spparks); }

void
SPPARKSUserObject::initialize()
{
  if (_spparks_only)
    return;

  // getSPPARKSData();
  // setSPPARKSData();
}

int
SPPARKSUserObject::getIntValue(unsigned int fem_node_id, unsigned int index) const
{
  return _spparks_only ? 0 : getValue(_int_data_for_fem, fem_node_id, index);
}

Real
SPPARKSUserObject::getDoubleValue(unsigned int fem_node_id, unsigned int index) const
{
  return _spparks_only ? 0.0 : getValue(_double_data_for_fem, fem_node_id, index);
}

char *
SPPARKSUserObject::runSPPARKSCommand(const std::string & cmd)
{
  return spparks_command(_spparks, cmd.c_str());
}

void
SPPARKSUserObject::getSPPARKSData()
{
  // Update the integer data
  for (unsigned int i = 0; i < _from_ivar.size(); ++i)
    getSPPARKSData(_int_data_for_fem[_from_ivar[i]], "iarray", _from_ivar[i]);

  // Update the double data
  for (unsigned int i = 0; i < _from_dvar.size(); ++i)
    getSPPARKSData(_double_data_for_fem[_from_dvar[i]], "darray", _from_dvar[i]);
}

void
SPPARKSUserObject::setSPPARKSData()
{
  // Update the integer data
  int * ip = nullptr;
  for (unsigned int i = 0; i < _to_ivar.size(); ++i)
    setSPPARKSData(ip, "iarray", _to_ivar[i], *_int_vars[i]);

  // Update the double data
  double * dp = nullptr;
  for (unsigned int i = 0; i < _to_dvar.size(); ++i)
    setSPPARKSData(dp, "darray", _to_dvar[i], *_double_vars[i]);
}

void
SPPARKSUserObject::setFEMData()
{
  // Update the double solving variables using the data from SPPARKS, added by YF
  for (unsigned int i = 0; i < _sol_vars.size(); ++i)
    setFEMData<double>("darray", _from_dvar[i], *_sol_vars[i]);
}

void
SPPARKSUserObject::execute()
{
  if ((_one_time_run && _n_spparks_run > 0) || _spparks_only || _t == _last_time)
    return;

  initSPPARKS();

  _last_time = _t;

  // set variables in SPPARKS defined by to_ivar and to_dvar
  setSPPARKSData();

  // Run SPPARKS over a certain time defined by sp_time
  const Real sp_time = getSPPARKSTime(_dt);
  std::stringstream cmd;
  cmd << "run ";
  cmd << sp_time;
  cmd << " pre no" << std::endl;
  runSPPARKSCommand(cmd.str());

  // times that SPPARKS has been called
  _n_spparks_run++;

  // obtain data from SPPARKS
  // getSPPARKSData update the auxvariables defined by from_ivar and from_dvar
  // setFEMData update the MOOSEvariables defined by sol__vars
  getSPPARKSData();
  setFEMData();
}

void
SPPARKSUserObject::initSPPARKS()
{
  // Set SPPARKS values based on 3.3.2 in Veena's paper
  // This is needed only to couple POTTS model with MARMOT

  // Loop over local FEM nodes
  // For each node, pick a random spin
  // Half alpha phase, half beta phase
  // Composition for alpha phase is 0.25; for beta phase, 0.75.

  if (!_init_spparks)
    return;

  if (_from_ivar.size() != 2)
    mooseError("Must have two integer variables from SPPARKS");

  if (_double_vars.size() != 1)
    mooseError("Must have one double variable to send to SPPARKS");

  SystemBase & sys = _double_vars[0]->sys();
  NumericVector<Number> & solution = sys.solution();

  std::map<libMesh::dof_id_type, int> ints[2];
  std::map<libMesh::dof_id_type, Real> comp;
  ConstNodeRange & node_range = *_fe_problem.mesh().getLocalNodeRange();
  for (auto i = node_range.begin(); i < node_range.end(); ++i)
  {
    int value = std::floor(100 * MooseRandom::rand()) + 1;
    ints[0][(*i)->id()] = value;
    ints[1][(*i)->id()] = value < 51;

    Real comp = 0.75;
    if (value < 51)
      comp = 0.25;

    // Set data
    solution.set((*i)->dof_number(sys.number(), _double_vars[0]->number(), 0), comp);
  }
  solution.close();

  int * pint = nullptr;
  char iarray[] = "iarray";
  for (unsigned int i = 0; i < 2; ++i)
  {
    getSPPARKSDataPointer(pint, iarray, _from_ivar[i]);

    // Index into data is SPPARKS node id.
    for (auto it = _fem_to_spparks.begin(); it != _fem_to_spparks.end(); ++it)
      pint[it->second.id] = ints[i][it->first.id];

    // Copy data across processors
    if (n_processors() > 1)
      sendRecvFEMData(ints[i], pint);
  }

  // do not re-init SPPARKS during this simulation
  _init_spparks = false;
}

void
SPPARKSUserObject::initialSetup()
{
  if (_spparks_only || _initialized)
    return;

  _initialized = true;

  // Initialize communication maps
  // _console << "\ninitialSetup: begin\n";

  // 1. Get on-processor map from SPPARKS ID to FEM ID
  int * iptr;
  getSPPARKSPointer(iptr, "nlocal");
  int nlcl = *iptr;

  double ** xyz_array;
  getSPPARKSPointer(xyz_array, "xyz");

  std::set<SPPARKSID> spparks_id; // Local SPPARKS nodes
  for (int i = 0; i < nlcl; ++i)
  {
    Point p(xyz_array[i][0], _dim > 1 ? xyz_array[i][1] : 0, _dim > 2 ? xyz_array[i][2] : 0);
    spparks_id.insert(SPPARKSID(i, p));
  }
  _num_local_spparks_nodes = spparks_id.size();

  // MOOSE FEM nodes are set up with periodic bcs (nodes at each end of periodicity),
  // but SPPARKS nodes are not.  We may have two or more MOOSE FEM nodes at SPPARKS
  // node locations.
  std::multiset<FEMID> fem_id; // Local MOOSE FEM nodes
  for (auto & node : *_fe_problem.mesh().getLocalNodeRange())
  {
    Point coor(*node);

    // TODO: this is horrible
    if (coor(0) == _xmax)
      coor(0) = _xmin;
    if (coor(1) == _ymax)
      coor(1) = _ymin;
    if (coor(2) == _zmax)
      coor(2) = _zmin;

    fem_id.insert(FEMID(node->id(), coor));
  }
  _num_local_fem_nodes = fem_id.size();

  // SPPARKS nodes not found on this processor
  std::set<SPPARKSID> unmatched_spparks;

  for (auto & id : spparks_id)
  {
    auto fem_iter = fem_id.find(id);
    if (fem_iter != fem_id.end())
      _spparks_to_fem.insert(std::pair<SPPARKSID, FEMID>(id, *fem_iter)); // SPPARKSID to FEMID
    else
      unmatched_spparks.insert(id);
  }

  for (auto & id : fem_id)
  {
    auto spparks_iter = spparks_id.find(id);
    if (spparks_iter != spparks_id.end())
      _fem_to_spparks.insert(std::pair<FEMID, SPPARKSID>(id, *spparks_iter)); // FEMID to SPPARKSID
    // else
    //   _spparks_to_proc.insert(std::pair<SPPARKSID, unsigned int>(*i, -1));
  }

  const unsigned int num_procs = n_processors();

  if (num_procs == 1)
  {
    if (spparks_id.size() != _spparks_to_fem.size())
    {
      // Error skipped to allow unmatched meshes
      // mooseError("Did not find MOOSE FEM node for each SPPARKS node, ", spparks_id.size()
      //           , ", ", fem_id.size(), ", ", _spparks_to_fem.size());

      _console << "\n  spparks size  " << spparks_id.size() << "not equal to fem size  "
               << fem_id.size() << '\n';
    }

    return;
  }

  // 2. Get send map (spparks id -> proc id)

  // A. Get local SPPARKS bounding box
  // TODO: Expand bounding boxes by half cell width?
  //       May not be necessary.

  unsigned int proc_id = processor_id();

  std::vector<Real> spparks_bounds(num_procs * 6, 0.0);
  unsigned int offset = proc_id * 6;
  Real * s_bounds = &spparks_bounds[0] + offset;
  s_bounds[0] = s_bounds[1] = s_bounds[2] = std::numeric_limits<Real>::max();
  s_bounds[3] = s_bounds[4] = s_bounds[5] = std::numeric_limits<Real>::min();
  for (auto & id : unmatched_spparks)
  {
    s_bounds[0] = std::min(s_bounds[0], id.coor(0));
    s_bounds[1] = std::min(s_bounds[1], id.coor(1));
    s_bounds[2] = std::min(s_bounds[2], id.coor(2));
    s_bounds[3] = std::max(s_bounds[3], id.coor(0));
    s_bounds[4] = std::max(s_bounds[4], id.coor(1));
    s_bounds[5] = std::max(s_bounds[5], id.coor(2));
  }
  libMesh::BoundingBox spparks_bb(Point(s_bounds[0], s_bounds[1], s_bounds[2]),
                                  Point(s_bounds[3], s_bounds[4], s_bounds[5]));
  _communicator.sum(spparks_bounds);

  //
  // B: Get MOOSE bounding boxes
  //

  //
  // TODO: These bboxes could/should be based on non-matched fem nodes.
  //
  std::vector<Real> fem_bounds(num_procs * 6, 0.0);
  Real * e_bounds = &fem_bounds[0] + offset;
  e_bounds[0] = e_bounds[1] = e_bounds[2] = std::numeric_limits<Real>::max();
  e_bounds[3] = e_bounds[4] = e_bounds[5] = std::numeric_limits<Real>::min();
  for (auto & id : fem_id)
  {
    e_bounds[0] = std::min(e_bounds[0], id.coor(0));
    e_bounds[1] = std::min(e_bounds[1], id.coor(1));
    e_bounds[2] = std::min(e_bounds[2], id.coor(2));
    e_bounds[3] = std::max(e_bounds[3], id.coor(0));
    e_bounds[4] = std::max(e_bounds[4], id.coor(1));
    e_bounds[5] = std::max(e_bounds[5], id.coor(2));
  }
  libMesh::BoundingBox fem_bb(Point(e_bounds[0], e_bounds[1], e_bounds[2]),
                              Point(e_bounds[3], e_bounds[4], e_bounds[5]));
  _communicator.sum(fem_bounds);

  //
  // C: Get number of processors that overlap my SPPARKS and MOOSE domains
  //

  std::vector<unsigned int> procs_overlapping_spparks_domain;
  std::vector<unsigned int> procs_overlapping_fem_domain;
  for (unsigned int i = 0; i < num_procs; ++i)
  {
    if (i == proc_id)
      continue;

    offset = i * 6;
    e_bounds = &fem_bounds[0] + offset;
    libMesh::BoundingBox e_box(Point(e_bounds[0], e_bounds[1], e_bounds[2]),
                               Point(e_bounds[3], e_bounds[4], e_bounds[5]));
    if (spparks_bb.intersects(e_box))
      procs_overlapping_spparks_domain.push_back(i);

    s_bounds = &spparks_bounds[0] + offset;
    libMesh::BoundingBox s_box(Point(s_bounds[0], s_bounds[1], s_bounds[2]),
                               Point(s_bounds[3], s_bounds[4], s_bounds[5]));
    if (fem_bb.intersects(s_box))
      procs_overlapping_fem_domain.push_back(i);
  }

  //
  // D: Communicate number of MOOSE FEM nodes, number of SPPARKS nodes
  //

  // Number remote fem nodes for each proc overlapping my spparks
  std::vector<unsigned int> num_fem_nodes(procs_overlapping_spparks_domain.size(), 0);
  std::vector<unsigned int> num_spparks_nodes(procs_overlapping_fem_domain.size(), 0);

  std::vector<MPI_Request> recv_request1(
      std::max(procs_overlapping_spparks_domain.size(), procs_overlapping_fem_domain.size()));
  std::vector<MPI_Request> recv_request2(
      std::max(procs_overlapping_spparks_domain.size(), procs_overlapping_fem_domain.size()));
  int comm_tag = 100;

  for (unsigned int i = 0; i < procs_overlapping_spparks_domain.size(); ++i)
    // if (num_fem_nodes.size() && procs_overlapping_spparks_domain.size())
    MPI_Irecv(&num_fem_nodes[i],
              1,
              MPI_UNSIGNED,
              procs_overlapping_spparks_domain[i],
              comm_tag,
              _communicator.get(),
              &recv_request1[i]);

  for (unsigned int i = 0; i < procs_overlapping_fem_domain.size(); ++i)
    // if (num_spparks_nodes.size() && procs_overlapping_fem_domain.size())
    MPI_Irecv(&num_spparks_nodes[i],
              1,
              MPI_UNSIGNED,
              procs_overlapping_fem_domain[i],
              comm_tag + 11,
              _communicator.get(),
              &recv_request2[i]);

  for (unsigned int i = 0; i < procs_overlapping_fem_domain.size(); ++i)
    // if (procs_overlapping_fem_domain.size())
    MPI_Send(&_num_local_fem_nodes,
             1,
             MPI_UNSIGNED,
             procs_overlapping_fem_domain[i],
             comm_tag,
             _communicator.get());

  for (unsigned int i = 0; i < procs_overlapping_spparks_domain.size(); ++i)
    // if (procs_overlapping_spparks_domain.size())
    MPI_Send(&_num_local_spparks_nodes,
             1,
             MPI_UNSIGNED,
             procs_overlapping_spparks_domain[i],
             comm_tag + 11,
             _communicator.get());

  std::vector<MPI_Status> recv_status1(
      std::max(procs_overlapping_spparks_domain.size(), procs_overlapping_fem_domain.size()));
  std::vector<MPI_Status> recv_status2(
      std::max(procs_overlapping_spparks_domain.size(), procs_overlapping_fem_domain.size()));
  MPI_Waitall(procs_overlapping_spparks_domain.size(), &recv_request1[0], &recv_status1[0]);
  MPI_Waitall(procs_overlapping_fem_domain.size(), &recv_request2[0], &recv_status2[0]);

  //
  // E: Communicate MOOSE FEM nodes, SPPARKS nodes
  //

  comm_tag = 200;
  int comm_tag_double = comm_tag + 1;
  const unsigned int num_fem_nodes_total =
      std::accumulate(&num_fem_nodes[0], &num_fem_nodes[0] + num_fem_nodes.size(), 0);
  const unsigned int num_spparks_nodes_total =
      std::accumulate(&num_spparks_nodes[0], &num_spparks_nodes[0] + num_spparks_nodes.size(), 0);
  std::vector<unsigned int> remote_fem_nodes(num_fem_nodes_total);
  std::vector<Real> remote_fem_coords(num_fem_nodes_total * 3);
  std::vector<unsigned int> remote_spparks_nodes(num_spparks_nodes_total);
  std::vector<Real> remote_spparks_coords(num_spparks_nodes_total * 3);
  offset = 0;

  // Receive MOOSE FEM nodes
  std::vector<MPI_Request> recv_request_coor1(procs_overlapping_spparks_domain.size());
  for (unsigned int i = 0; i < procs_overlapping_spparks_domain.size(); ++i)
  {
    MPI_Irecv(&remote_fem_nodes[offset],
              num_fem_nodes[i],
              MPI_UNSIGNED,
              procs_overlapping_spparks_domain[i],
              comm_tag,
              _communicator.get(),
              &recv_request1[i]);
    MPI_Irecv(&remote_fem_coords[offset * 3],
              num_fem_nodes[i] * 3,
              MPI_DOUBLE,
              procs_overlapping_spparks_domain[i],
              comm_tag_double,
              _communicator.get(),
              &recv_request_coor1[i]);
    offset += num_fem_nodes[i];
  }

  // Receive SPPARKS nodes
  std::vector<MPI_Request> recv_request_coor2(procs_overlapping_fem_domain.size());
  offset = 0;
  for (unsigned int i = 0; i < procs_overlapping_fem_domain.size(); ++i)
  {
    MPI_Irecv(&remote_spparks_nodes[offset],
              num_spparks_nodes[i],
              MPI_UNSIGNED,
              procs_overlapping_fem_domain[i],
              comm_tag + 11,
              _communicator.get(),
              &recv_request2[i]);
    MPI_Irecv(&remote_spparks_coords[offset * 3],
              num_spparks_nodes[i] * 3,
              MPI_DOUBLE,
              procs_overlapping_fem_domain[i],
              comm_tag_double + 11,
              _communicator.get(),
              &recv_request_coor2[i]);
    offset += num_spparks_nodes[i];
  }

  // Prepare vectors of MOOSE FEM ids, coordinates
  std::vector<unsigned int> fem_ids(_num_local_fem_nodes);
  std::vector<Real> fem_coords(3 * _num_local_fem_nodes);
  offset = 0;
  for (auto & id : fem_id)
  {
    fem_ids[offset] = id.id;
    fem_coords[offset * 3 + 0] = id.coor(0);
    fem_coords[offset * 3 + 1] = id.coor(1);
    fem_coords[offset * 3 + 2] = id.coor(2);
    ++offset;
  }

  // Send MOOSE FEM ids, coordinates
  offset = 0;
  for (unsigned int i = 0; i < procs_overlapping_fem_domain.size(); ++i)
  {
    MPI_Send(&fem_ids[0],
             _num_local_fem_nodes,
             MPI_UNSIGNED,
             procs_overlapping_fem_domain[i],
             comm_tag,
             _communicator.get());
    MPI_Send(&fem_coords[0],
             _num_local_fem_nodes * 3,
             MPI_DOUBLE,
             procs_overlapping_fem_domain[i],
             comm_tag_double,
             _communicator.get());
    ++offset;
  }

  // Prepare SPPARKS ids, coordinates
  std::vector<unsigned int> spparks_ids(_num_local_spparks_nodes);
  std::vector<Real> spparks_coords(3 * _num_local_spparks_nodes);
  offset = 0;
  for (auto & id : spparks_id)
  {
    spparks_ids[offset] = id.id;
    spparks_coords[offset * 3 + 0] = id.coor(0);
    spparks_coords[offset * 3 + 1] = id.coor(1);
    spparks_coords[offset * 3 + 2] = id.coor(2);
    ++offset;
  }

  // Send SPPARKS ids, coordinates
  offset = 0;
  for (unsigned int i = 0; i < procs_overlapping_spparks_domain.size(); ++i)
  {
    MPI_Send(&spparks_ids[0],
             _num_local_spparks_nodes,
             MPI_UNSIGNED,
             procs_overlapping_spparks_domain[i],
             comm_tag + 11,
             _communicator.get());
    MPI_Send(&spparks_coords[0],
             _num_local_spparks_nodes * 3,
             MPI_DOUBLE,
             procs_overlapping_spparks_domain[i],
             comm_tag_double + 11,
             _communicator.get());
    ++offset;
  }

  MPI_Waitall(procs_overlapping_spparks_domain.size(), &recv_request1[0], &recv_status1[0]);
  MPI_Waitall(procs_overlapping_spparks_domain.size(), &recv_request_coor1[0], &recv_status1[0]);
  MPI_Waitall(procs_overlapping_fem_domain.size(), &recv_request2[0], &recv_status2[0]);
  MPI_Waitall(procs_overlapping_fem_domain.size(), &recv_request_coor2[0], &recv_status2[0]);

  //
  // F: Count matching nodes for each proc that sent MOOSE FEM nodes, SPPARKS nodes
  //

  comm_tag = 300;
  // Receive number of MOOSE FEM nodes on this processor that match SPPARKS nodes on another
  std::vector<unsigned int> num_remote_fem_matches(procs_overlapping_fem_domain.size());
  for (unsigned int i = 0; i < procs_overlapping_fem_domain.size(); ++i)
    MPI_Irecv(&num_remote_fem_matches[i],
              1,
              MPI_UNSIGNED,
              procs_overlapping_fem_domain[i],
              comm_tag,
              _communicator.get(),
              &recv_request1[i]);

  // Receive number of SPPARKS nodes on this processor that match MOOSE FEM nodes on another
  std::vector<unsigned int> num_remote_spparks_matches(procs_overlapping_spparks_domain.size());
  for (unsigned int i = 0; i < procs_overlapping_spparks_domain.size(); ++i)
    MPI_Irecv(&num_remote_spparks_matches[i],
              1,
              MPI_UNSIGNED,
              procs_overlapping_spparks_domain[i],
              comm_tag + 11,
              _communicator.get(),
              &recv_request2[i]);

  // Count number of remote MOOSE FEM nodes that match the processor's SPPARKS nodes
  // Store the matches
  std::vector<std::vector<unsigned int>> spparks_matches(procs_overlapping_spparks_domain.size());
  offset = 0;
  for (unsigned int i = 0; i < procs_overlapping_spparks_domain.size(); ++i)
    for (unsigned int j = 0; j < num_fem_nodes[i]; ++j)
    {
      SPPARKSID tmp(remote_fem_nodes[offset],
                    Point(remote_fem_coords[offset * 3 + 0],
                          remote_fem_coords[offset * 3 + 1],
                          remote_fem_coords[offset * 3 + 2]));
      auto iter = spparks_id.find(tmp);
      if (iter != spparks_id.end())
      {
        spparks_matches[i].push_back(tmp.id);
        _spparks_to_proc[*iter].push_back(procs_overlapping_spparks_domain[i]);
      }

      ++offset;
    }

  // Count number of remote SPPARKS nodes that match this processor's MOOSE FEM nodes
  // Store the matches
  std::vector<std::vector<unsigned int>> fem_matches(procs_overlapping_fem_domain.size());
  offset = 0;
  for (unsigned int i = 0; i < procs_overlapping_fem_domain.size(); ++i)
    for (unsigned int j = 0; j < num_spparks_nodes[i]; ++j)
    {
      FEMID tmp(remote_spparks_nodes[offset],
                Point(remote_spparks_coords[offset * 3 + 0],
                      remote_spparks_coords[offset * 3 + 1],
                      remote_spparks_coords[offset * 3 + 2]));
      auto iter = fem_id.find(tmp);
      if (iter != fem_id.end())
      {
        fem_matches[i].push_back(tmp.id);
        _fem_to_proc[*iter].push_back(procs_overlapping_fem_domain[i]);
      }

      ++offset;
    }

  // Send number of MOOSE FEM nodes that match this processor's SPPARKS nodes
  std::vector<unsigned int> spparks_sizes(procs_overlapping_spparks_domain.size());
  for (unsigned int i = 0; i < procs_overlapping_spparks_domain.size(); ++i)
  {
    spparks_sizes[i] = spparks_matches[i].size();
    MPI_Send(&spparks_sizes[i],
             1,
             MPI_UNSIGNED,
             procs_overlapping_spparks_domain[i],
             comm_tag,
             _communicator.get());
  }

  // Send number of SPPARKS nodes that match this processor's MOOSE FEM nodes
  std::vector<unsigned int> fem_sizes(procs_overlapping_fem_domain.size());
  for (unsigned int i = 0; i < procs_overlapping_fem_domain.size(); ++i)
  {
    fem_sizes[i] = fem_matches[i].size();
    MPI_Send(&fem_sizes[i],
             1,
             MPI_UNSIGNED,
             procs_overlapping_fem_domain[i],
             comm_tag + 11,
             _communicator.get());
  }

  MPI_Waitall(procs_overlapping_fem_domain.size(), &recv_request1[0], &recv_status1[0]);
  MPI_Waitall(procs_overlapping_spparks_domain.size(), &recv_request2[0], &recv_status2[0]);

  //
  // G: Communicate matching nodes
  //

  comm_tag = 400;
  // Receive MOOSE DEM ids that match SPPARKS nodes on another processor
  std::vector<std::vector<unsigned int>> matched_fem_ids(procs_overlapping_fem_domain.size());
  for (unsigned int i = 0; i < procs_overlapping_fem_domain.size(); ++i)
  {
    matched_fem_ids[i].resize(num_remote_fem_matches[i]);
    if (matched_fem_ids[i].size())
      MPI_Irecv(&matched_fem_ids[i][0],
                num_remote_fem_matches[i],
                MPI_UNSIGNED,
                procs_overlapping_fem_domain[i],
                comm_tag,
                _communicator.get(),
                &recv_request1[i]);
  }

  // Receive SPPARKS ids that match MOOSE FEM nodes on another processor
  std::vector<std::vector<unsigned int>> matched_spparks_ids(
      procs_overlapping_spparks_domain.size());
  for (unsigned int i = 0; i < procs_overlapping_spparks_domain.size(); ++i)
  {
    matched_spparks_ids[i].resize(num_remote_spparks_matches[i]);
    if (matched_spparks_ids[i].size())
      MPI_Irecv(&matched_spparks_ids[i][0],
                num_remote_spparks_matches[i],
                MPI_UNSIGNED,
                procs_overlapping_spparks_domain[i],
                comm_tag + 11,
                _communicator.get(),
                &recv_request2[i]);
  }

  // Send remote MOOSE FEM ids that match SPPARKS nodes on this processor
  for (unsigned int i = 0; i < procs_overlapping_spparks_domain.size(); ++i)
  {
    spparks_matches[i].resize(spparks_sizes[i]);
    if (spparks_matches[i].size())
      MPI_Send(&spparks_matches[i][0],
               spparks_sizes[i],
               MPI_UNSIGNED,
               procs_overlapping_spparks_domain[i],
               comm_tag,
               _communicator.get());
  }
  // Send remote SPPARKS ids that match MOOSE FEM nodes on this processor
  for (unsigned int i = 0; i < procs_overlapping_fem_domain.size(); ++i)
  {
    fem_matches[i].resize(fem_sizes[i]);
    if (fem_matches[i].size())
      MPI_Send(&fem_matches[i][0],
               fem_sizes[i],
               MPI_UNSIGNED,
               procs_overlapping_fem_domain[i],
               comm_tag + 11,
               _communicator.get());
  }

  MPI_Waitall(procs_overlapping_fem_domain.size(), &recv_request1[0], &recv_status1[0]);
  MPI_Waitall(procs_overlapping_spparks_domain.size(), &recv_request2[0], &recv_status2[0]);

  //
  // H: Generate final recv communication maps
  //

  for (unsigned int i = 0; i < procs_overlapping_fem_domain.size(); ++i)
    for (unsigned int j = 0; j < num_remote_fem_matches[i]; ++j)
      _sending_proc_to_fem_id[procs_overlapping_fem_domain[i]].push_back(matched_fem_ids[i][j]);

  for (unsigned int i = 0; i < procs_overlapping_spparks_domain.size(); ++i)
    for (unsigned int j = 0; j < num_remote_spparks_matches[i]; ++j)
      _sending_proc_to_spparks_id[procs_overlapping_spparks_domain[i]].push_back(
          matched_spparks_ids[i][j]);
}
