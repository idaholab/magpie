#ifndef SPPARKSUSEROBJECT_H
#define SPPARKSUSEROBJECT_H

#include "GeneralUserObject.h"
#include "spparks/src/library.h"

// forward declarations
class SPPARKSUserObject;

template<>
InputParameters validParams<SPPARKSUserObject>();

class SPPARKSUserObject : public GeneralUserObject
{
public:
  SPPARKSUserObject(const std::string & name, InputParameters parameters);
  virtual ~SPPARKSUserObject();

  virtual void initialSetup();
  virtual void residualSetup() {}
  virtual void timestepSetup() {}

  virtual void initialize();
  virtual void execute();
  virtual void finalize() {}

  int  getIntValue(unsigned elk_node_id, unsigned index) const;
  Real getDoubleValue(unsigned elk_node_id, unsigned index) const;

protected:
  struct FEMID;
  typedef FEMID SPPARKSID;

  void initSPPARKS();

  char * runSPPARKSCommand(const std::string & cmd);

  Real getSPPARKSTime(Real dt) { return dt / _time_ratio; }

  template <typename T>
  Real getValue(const std::map<unsigned, std::map<unsigned, T> > & data, unsigned elk_node_id, unsigned index) const;

  template <typename T>
  void getSPPARKSDataPointer(T *& ptr, const std::string & name, unsigned index);

  template <typename T>
  void getSPPARKSPointer(T *& ptr, const std::string & name) const;

  void getSPPARKSData();
  void setSPPARKSData();

  template <typename T>
  void getSPPARKSData(std::map<unsigned, T> & storage, const std::string & name, unsigned index);

  template <typename T>
  void sendRecvSPPARKSData(const T * const data, std::map<unsigned, T> & storage);

  template <typename T>
  void setSPPARKSData(T * data, const std::string & name, unsigned index, MooseVariable & var);

  // Initiate ELK based on data from SPPARKS, added by YF.
  void setELKData();

  template <typename T>
  void setELKData(const std::string & name, unsigned index, MooseVariable & var);

  template <typename T>
  void sendRecvELKData(const std::map<libMesh::dof_id_type, T> & storage, T * const data);

  /// SPPARKS instance pointer
  void * _spparks;

  const std::string & _file;
  const bool _spparks_only;
  const std::vector<unsigned> _from_ivar;
  const std::vector<unsigned> _from_dvar;
  const std::vector<unsigned> _to_ivar;
  const std::vector<unsigned> _to_dvar;
  const Real _xmin;
  const Real _ymin;
  const Real _zmin;
  const Real _xmax;
  const Real _ymax;
  const Real _zmax;
  bool _init_spparks;
  const Real _time_ratio;
  bool _initialized;

  // added by YF, do one time SPPARKS run
  bool _one_time_run;
  int _times_of_run;

  int _dim;

  std::vector<MooseVariable*> _int_vars;
  std::vector<MooseVariable*> _double_vars;
  std::vector<MooseVariable*> _sol_vars; //added by YF

  // Communication maps
  std::map<SPPARKSID, std::vector<unsigned> > _spparks_to_proc; // SPPARKSID to vector of procs that need the value
  std::map<unsigned, std::vector<libMesh::dof_id_type> > _sending_proc_to_elk_id; // Processor to list of ELK ids
  std::map<FEMID, std::vector<unsigned> > _elk_to_proc; // FEMID to vector of procs that need the value
  std::map<unsigned, std::vector<unsigned> > _sending_proc_to_spparks_id; // Processor to list of SPPARKS ids

  unsigned _num_local_elk_nodes;
  unsigned _num_local_spparks_nodes;

  std::map<SPPARKSID, FEMID> _spparks_to_elk;      // Local SPPARKSID to local FEMID
  std::multimap<FEMID, SPPARKSID> _elk_to_spparks; // Local FEMID to local SPPARKSID

  // Maps from variable index to ELK node id to value
  std::map<unsigned, std::map<unsigned, int> > _int_data_for_elk;
  std::map<unsigned, std::map<unsigned, double> > _double_data_for_elk;

  // Maps from variable index to value (vector index is array index)
  std::map<unsigned, std::vector<int> > _int_data_for_spparks;
  std::map<unsigned, std::vector<double> > _double_data_for_spparks;

  Real _last_time;
};

struct SPPARKSUserObject::FEMID
{
  FEMID(libMesh::dof_id_type ident, const Point & p) :
    id(ident), coor(p)
  {}

  libMesh::dof_id_type id;
  Point coor;
  bool operator<(const FEMID & rhs) const
  {
    // return coor < rhs.coor;
    const Real tol = 1e-12;
    if (!equalCoor(coor(0), rhs.coor(0), tol))
      return (coor(0) + tol) < rhs.coor(0);
    if (!equalCoor(coor(1), rhs.coor(1), tol))
      return (coor(1) + tol) < rhs.coor(1);
    return (coor(2) + tol) < rhs.coor(2);
  }

private:
  bool equalCoor(Real f, Real s, Real tol) const
  {
    return f < (s+tol) && (f > s-tol);
  }
};

template <typename T>
Real
SPPARKSUserObject::getValue(const std::map<unsigned, std::map<unsigned, T> > & data, unsigned elk_node_id, unsigned index) const
{
  // Extract the data
  const typename std::map<unsigned, std::map<unsigned, T> >::const_iterator it = data.find(index);
  if (it == data.end())
    mooseError("SPPARKSUserObject error: unknown index " << index);

  const typename std::map<unsigned, T>::const_iterator it2 = it->second.find(elk_node_id);
  if (it2 == it->second.end())
  {
    mooseWarning("SPPARKSUserObject error: unknown elk node id " << elk_node_id);
    return 0;
  }

  return it2->second;
}

template <typename T>
void
SPPARKSUserObject::getSPPARKSDataPointer(T *& ptr, const std::string & name, unsigned index)
{
  std::stringstream fullname;
  fullname << name;
  fullname << index;
  getSPPARKSPointer(ptr, fullname.str().c_str());
}

template <typename T>
void
SPPARKSUserObject::getSPPARKSPointer(T *& ptr, const std::string & name) const
{
  void * p = spparks_extract(_spparks, name.c_str());
  if (!p)
    mooseError("SPPARKS returned NULL pointer for " << name);
  ptr = static_cast<T*>(p);
}

template <typename T>
void
SPPARKSUserObject::getSPPARKSData(std::map<unsigned, T> & storage, const std::string & name, unsigned index)
{
  T * data;
  getSPPARKSDataPointer(data, name.c_str(), index);

  // Copy data from local SPPARKS node to local ELK node
  // Index into storage is ELK node id.
  for (std::multimap<FEMID, SPPARKSID>::const_iterator i = _elk_to_spparks.begin(); i != _elk_to_spparks.end(); ++i)
    storage[i->first.id] = data[i->second.id];

  // Copy data across processors
  if (n_processors() > 1)
    sendRecvSPPARKSData(data, storage);
}

template <typename T>
void
SPPARKSUserObject::sendRecvSPPARKSData(const T * const data, std::map<unsigned, T> & storage)
{
  Parallel::MessageTag comm_tag(101);

  const unsigned num_recvs = _sending_proc_to_elk_id.size();
  std::vector<Parallel::Request> recv_request(num_recvs);

  std::map<unsigned, std::vector<T> > data_to_me; // sending proc, vector of SPPARKS values (one value per SPPARKS node)
  unsigned offset = 0;
  for (std::map<unsigned, std::vector<libMesh::dof_id_type> >::const_iterator i = _sending_proc_to_elk_id.begin();
       i != _sending_proc_to_elk_id.end(); ++i)
  {
    data_to_me[i->first].resize(i->second.size());
    _communicator.receive(i->first, data_to_me[i->first], recv_request[offset], comm_tag);
    ++offset;
  }

  std::map<unsigned, std::vector<T> > data_from_me; // Processor, list of SPPARKS values
  for (std::map<SPPARKSID, std::vector<unsigned> >::const_iterator i = _spparks_to_proc.begin();
       i != _spparks_to_proc.end(); ++i)
  {
    for (unsigned j = 0; j < i->second.size(); ++j)
      data_from_me[i->second[j]].push_back(data[i->first.id]);
  }

  for (typename std::map<unsigned, std::vector<T> >::const_iterator i = data_from_me.begin(); i != data_from_me.end(); ++i)
    _communicator.send(i->first, data_from_me[i->first], comm_tag);

  Parallel::wait(recv_request);

  // Move data into storage
  for (std::map<unsigned, std::vector<libMesh::dof_id_type> >::const_iterator i = _sending_proc_to_elk_id.begin();
       i != _sending_proc_to_elk_id.end(); ++i)
  {
    const std::vector<libMesh::dof_id_type> & id = i->second;
    const std::vector<T> & v = data_to_me[i->first];

    // storage is ELK node id, value
    for (unsigned j = 0; j < v.size(); ++j)
      storage[id[j]] = v[j];
  }
}

template <typename T>
void
SPPARKSUserObject::setSPPARKSData(T * data, const std::string & name, unsigned index, MooseVariable & var)
{
  getSPPARKSDataPointer(data, name, index);

  SystemBase & sys = var.sys();
  NumericVector<Number> & solution = sys.solution();

  // Extract ELK data
  std::map<libMesh::dof_id_type, T> elk_data;
  ConstNodeRange & node_range = *_fe_problem.mesh().getLocalNodeRange();
  for (ConstNodeRange::const_iterator i = node_range.begin(); i < node_range.end(); ++i)
  {
    // Get data
    const Real value = solution((*i)->dof_number(sys.number(), var.number(), 0));

    elk_data[(*i)->id()] = value;
  }

  // Index into data is SPPARKS node id.
  for (std::multimap<FEMID, SPPARKSID>::const_iterator i = _elk_to_spparks.begin(); i != _elk_to_spparks.end(); ++i)
    data[i->second.id] = elk_data[i->first.id];

  // Copy data across processors
  if (n_processors() > 1)
    sendRecvELKData(elk_data, data);
}

template <typename T>
void
SPPARKSUserObject::setELKData(const std::string & name, unsigned index, MooseVariable & var)
{
  // get SPPARKS data
  std::map<unsigned, T> spparks_data;
  getSPPARKSData(spparks_data, name, index);

  SystemBase & sys = var.sys();
  NumericVector<Number> & solution = sys.solution();

  // Extract ELK data
  ConstNodeRange & node_range = *_fe_problem.mesh().getLocalNodeRange();

  // Set data
  for (ConstNodeRange::const_iterator i = node_range.begin(); i < node_range.end(); ++i)
    solution.set((*i)->dof_number(sys.number(), var.number(), 0), spparks_data[(*i)->id()]);

  solution.close();
}

// template <typename T>
// void
// SPPARKSUserObject::setELKData(T * data, const std::string & name, unsigned index, MooseVariable & var)
// {
//   getSPPARKSDataPointer(data, name, index);
//
//   SystemBase & sys = var.sys();
//   NumericVector<Number> & solution = sys.solution();
//
//   // Extract ELK data
//   std::map<libMesh::dof_id_type, T> spparks_data;
//   ConstNodeRange & node_range = *_fe_problem.mesh().getLocalNodeRange();
//
//   // Index into data is SPPARKS node id.
//   for (std::multimap<FEMID, SPPARKSID>::const_iterator i = _spparks_to_elk.begin(); i != _spparks_to_elk.end(); ++i)
//     spparks_data[i->first.id] = data[i->second.id];
//
//   // set data
//   for (ConstNodeRange::const_iterator i = node_range.begin(); i < node_range.end(); ++i)
//     solution.set((*i)->dof_number(sys.number(), var.number(), 0), spparks_data[(*i)->id()]);
//
//   // Copy data across processors
//   if (n_processors() > 1)
//     sendRecvSPPARKSData(data, spparks_data);
// }

template <typename T>
void
SPPARKSUserObject::sendRecvELKData(const std::map<libMesh::dof_id_type, T> & storage, T * const data)
{
  Parallel::MessageTag comm_tag(101);

  const unsigned num_recvs = _sending_proc_to_spparks_id.size();
  std::vector<Parallel::Request> recv_request(num_recvs);

  // sending proc, vector of ELK values (one value per ELK node)
  std::map<unsigned, std::vector<T> > data_to_me;
  unsigned offset = 0;
  for (std::map<unsigned, std::vector<unsigned> >::const_iterator i = _sending_proc_to_spparks_id.begin();
       i != _sending_proc_to_spparks_id.end(); ++i)
  {
    data_to_me[i->first].resize(i->second.size());
    _communicator.receive(i->first, data_to_me[i->first], recv_request[offset], comm_tag);
    ++offset;
  }

  std::map<unsigned, std::vector<T> > data_from_me; // Processor, list of ELK values
  for (std::map<FEMID, std::vector<unsigned> >::const_iterator i = _elk_to_proc.begin(); i != _elk_to_proc.end(); ++i)
  {
    for (unsigned j = 0; j < i->second.size(); ++j)
      data_from_me[i->second[j]].push_back(storage.find(i->first.id)->second);
  }

  for (typename std::map<unsigned, std::vector<T> >::const_iterator i = data_from_me.begin(); i != data_from_me.end(); ++i)
    _communicator.send(i->first, data_from_me[i->first], comm_tag);

  Parallel::wait(recv_request);

  // Move data into storage
  for (std::map<unsigned, std::vector<unsigned> >::const_iterator i = _sending_proc_to_spparks_id.begin();
       i != _sending_proc_to_spparks_id.end(); ++i)
  {
    const std::vector<unsigned> & id = i->second;
    const std::vector<T> & v = data_to_me[i->first];
    if (id.size() != v.size())
      mooseError("Mismatched communication vectors");

    // data is SPPARKS index, value
    for (unsigned j = 0; j < v.size(); ++j)
      data[id[j]] = v[j];
  }
}

#endif
