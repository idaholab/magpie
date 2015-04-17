#ifndef SPPARKSUSEROBJECT_H
#define SPPARKSUSEROBJECT_H

#include "GeneralUserObject.h"

#include "SPPARKS/src/library.h"

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

  void initSPPARKS();

  char * runSPPARKSCommand( const std::string & cmd )
  {
    std::vector<char> strng(cmd.begin(), cmd.end());
    strng.push_back('\0');
    return spparks_command(_spparks, &strng[0]);
  }

  Real getSPPARKSTime( Real dt )
  {
    return dt / _time_ratio;
  }

  template
  <typename T>
  Real getValue(const std::map<unsigned, std::map<unsigned, T> > & data, unsigned elk_node_id, unsigned index) const
  {
    // Extract the data
    const typename std::map<unsigned, std::map<unsigned, T> >::const_iterator it = data.find(index);
    if (it == data.end())
    {
      std::stringstream err;
      err << "SPPARKSUserObject error: unknown index ";
      err << index << std::endl;
      mooseError(err.str());
    }
    const typename std::map<unsigned, T>::const_iterator it2 = it->second.find(elk_node_id);
    if (it2 == it->second.end())
    {
      std::stringstream err;
      err << "SPPARKSUserObject error: unknown elk node id ";
      err << elk_node_id << std::endl;
      mooseWarning(err.str());
      return 0;
    }
    return it2->second;
  }


  template
  <typename T>
  void getSPPARKSDataPointer( T *& ptr, char * string, unsigned i )
  {
    std::stringstream name;
    name << string;
    name << i;
    std::string s = name.str();
    std::vector<char> strng(s.c_str(), s.c_str()+s.size()+1u);
    getSPPARKSPointer(ptr, &strng[0]);
  }

  template
  <typename T>
  void getSPPARKSPointer( T *& ptr, char * string ) const
  {
    void * p = spparks_extract(_spparks, string);
    if (!p)
    {
      mooseError("SPPARKS returned NULL pointer for " + std::string(string));
    }
    ptr = static_cast<T*>(p);
  }

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

  struct ELKID
  {
    ELKID( libMesh::dof_id_type ident, const Point & p ) :
      id(ident),
      coor(p)
    {}
    libMesh::dof_id_type id;
    Point coor;
    bool operator<(const ELKID & rhs) const
    {
      // return coor < rhs.coor;
      const Real tol = 1e-12;
      if ( !equalCoor( coor(0), rhs.coor(0), tol ) )
      {
        return (coor(0) + tol) < rhs.coor(0);
      }
      if ( !equalCoor( coor(1), rhs.coor(1), tol ) )
      {
        return (coor(1) + tol) < rhs.coor(1);
      }
      return (coor(2) + tol) < rhs.coor(2);
    }
  private:
    bool equalCoor( Real f, Real s, Real tol ) const
    {
      return f < (s+tol) && (f > s-tol);
    }
  };
  typedef ELKID SPPARKSID;

  void getSPPARKSData();
  void setSPPARKSData();

  template
  <typename T>
  void getSPPARKSData( std::map<unsigned, T> & storage, char * string, unsigned index )
  {
    T * data;
    getSPPARKSDataPointer( data, string, index );

    // Copy data from local SPPARKS node to local ELK node
    for (std::multimap<ELKID, SPPARKSID>::const_iterator i = _elk_to_spparks.begin(); i != _elk_to_spparks.end(); ++i)
    {
      // Index into storage is ELK node id.
      storage[i->first.id] = data[i->second.id];
    }

    if ( n_processors() > 1 )
    {
      // Copy data across processors
      sendRecvSPPARKSData( data, storage );
    }

  }

  template
  <typename T>
  void sendRecvSPPARKSData( const T * const data, std::map<unsigned, T> & storage )
  {

    Parallel::MessageTag comm_tag(101);

    const unsigned num_recvs = _sending_proc_to_elk_id.size();
    std::vector<Parallel::Request> recv_request(num_recvs);

    std::map<unsigned, std::vector<T> > data_to_me; // sending proc, vector of SPPARKS values (one value per SPPARKS node)
    unsigned offset = 0;
    for (std::map<unsigned, std::vector<libMesh::dof_id_type> >::const_iterator i = _sending_proc_to_elk_id.begin();
         i != _sending_proc_to_elk_id.end(); ++i)
    {
      data_to_me[i->first].resize( i->second.size() );
      _communicator.receive(i->first, data_to_me[i->first], recv_request[offset], comm_tag);
      ++offset;
    }


    std::map<unsigned, std::vector<T> > data_from_me; // Processor, list of SPPARKS values
    for (std::map<SPPARKSID, std::vector<unsigned> >::const_iterator i = _spparks_to_proc.begin(); i != _spparks_to_proc.end(); ++i)
    {
      for (unsigned j = 0; j < i->second.size(); ++j)
      {
        data_from_me[i->second[j]].push_back( data[i->first.id] );
      }
    }

    for (typename std::map<unsigned, std::vector<T> >::const_iterator i = data_from_me.begin(); i != data_from_me.end(); ++i)
    {
      _communicator.send(i->first, data_from_me[i->first], comm_tag);
    }

    Parallel::wait(recv_request);

    // Move data into storage
    for (std::map<unsigned, std::vector<libMesh::dof_id_type> >::const_iterator i = _sending_proc_to_elk_id.begin();
         i != _sending_proc_to_elk_id.end(); ++i)
    {
      const std::vector<libMesh::dof_id_type> & id = i->second;
      const std::vector<T> & v = data_to_me[i->first];
      for (unsigned j = 0; j < v.size(); ++j)
      {
        // storage is ELK node id, value
        storage[id[j]] = v[j];
      }
    }
  }


  template
  <typename T>
  void setSPPARKSData( T * data, char * string, unsigned index, MooseVariable & var )
  {

    getSPPARKSDataPointer( data, string, index );

    SystemBase & sys = var.sys();
    NumericVector<Number> & solution = sys.solution();

    // Extract ELK data
    std::map<libMesh::dof_id_type, T> elk_data;
    ConstNodeRange & node_range = *_fe_problem.mesh().getLocalNodeRange();
    for ( ConstNodeRange::const_iterator i = node_range.begin(); i < node_range.end(); ++i )
    {
      // Get data
      const Real value = solution( (*i)->dof_number(sys.number(), var.number(), 0) );

      elk_data[(*i)->id()] = value;

    }
    for (std::multimap<ELKID, SPPARKSID>::const_iterator i = _elk_to_spparks.begin(); i != _elk_to_spparks.end(); ++i)
    {
      // Index into data is SPPARKS node id.
      data[i->second.id] = elk_data[i->first.id];
    }

    if ( n_processors() > 1 )
    {
      // Copy data across processors
      sendRecvELKData( elk_data, data );
    }

  }

    // Initiate ELK based on data from SPPARKS, added by YF. 
  void setELKData();

  template
  <typename T>
  void setELKData( char * string, unsigned index, MooseVariable & var )
  {

    // get spparks data 
    std::map<unsigned, T> spparks_data;
    getSPPARKSData( spparks_data, string, index );

    SystemBase & sys = var.sys();
    NumericVector<Number> & solution = sys.solution();

    // Extract ELK data
    ConstNodeRange & node_range = *_fe_problem.mesh().getLocalNodeRange();

    for ( ConstNodeRange::const_iterator i = node_range.begin(); i < node_range.end(); ++i )
    {
      // Set data
      // std::cout << std::endl
      //        << "SetElkData " << spparks_data[(*i)->id()] << std::endl;
      solution.set( (*i)->dof_number(sys.number(), var.number(), 0), spparks_data[(*i)->id()] );
      // solution.set( (*i)->dof_number(sys.number(), var.number(), 0), 0.0 );

    }
    solution.close();
  }
/*
  void setELKData( T * data, char * string, unsigned index, MooseVariable & var )
  {

    getSPPARKSDataPointer( data, string, index );

    SystemBase & sys = var.sys();
    NumericVector<Number> & solution = sys.solution();

    // Extract ELK data
    std::map<libMesh::dof_id_type, T> spparks_data;
    ConstNodeRange & node_range = *_fe_problem.mesh().getLocalNodeRange();

    for (std::multimap<ELKID, SPPARKSID>::const_iterator i = _spparks_to_elk.begin(); i != _spparks_to_elk.end(); ++i)
    {
      // Index into data is SPPARKS node id.
      spparks_data[i->first.id] = data[i->second.id];
    }
    for ( ConstNodeRange::const_iterator i = node_range.begin(); i < node_range.end(); ++i )
    {
      // Set data
      solution.set( (*i)->dof_number(sys.number(), var.number(), 0), spparks_data[(*i)->id()] );

    }

    if ( n_processors() > 1 )
    {
      // Copy data across processors
      sendRecvSPPARKSData( data, spparks_data );
    }

  }
*/

  template
  <typename T>
  void sendRecvELKData( const std::map<libMesh::dof_id_type, T> & storage, T * const data )
  {

    Parallel::MessageTag comm_tag(101);

    const unsigned num_recvs = _sending_proc_to_spparks_id.size();
    std::vector<Parallel::Request> recv_request(num_recvs);

    std::map<unsigned, std::vector<T> > data_to_me; // sending proc, vector of ELK values (one value per ELK node)
    unsigned offset = 0;
    for (std::map<unsigned, std::vector<unsigned> >::const_iterator i = _sending_proc_to_spparks_id.begin();
         i != _sending_proc_to_spparks_id.end(); ++i)
    {
      data_to_me[i->first].resize( i->second.size() );
      _communicator.receive(i->first, data_to_me[i->first], recv_request[offset], comm_tag);
      ++offset;
    }


    std::map<unsigned, std::vector<T> > data_from_me; // Processor, list of ELK values
    for (std::map<ELKID, std::vector<unsigned> >::const_iterator i = _elk_to_proc.begin(); i != _elk_to_proc.end(); ++i)
    {
      for (unsigned j = 0; j < i->second.size(); ++j)
      {
        data_from_me[i->second[j]].push_back( storage.find(i->first.id)->second );
      }
    }

    for (typename std::map<unsigned, std::vector<T> >::const_iterator i = data_from_me.begin(); i != data_from_me.end(); ++i)
    {
      _communicator.send(i->first, data_from_me[i->first], comm_tag);
    }

    Parallel::wait(recv_request);

    // Move data into storage
    for (std::map<unsigned, std::vector<unsigned> >::const_iterator i = _sending_proc_to_spparks_id.begin();
         i != _sending_proc_to_spparks_id.end(); ++i)
    {
      const std::vector<unsigned> & id = i->second;
      const std::vector<T> & v = data_to_me[i->first];
      if ( id.size() != v.size() )
      {
        mooseError("Mismatched communication vectors");
      }
      for (unsigned j = 0; j < v.size(); ++j)
      {
        // data is SPPARKS index, value
        data[id[j]] = v[j];
      }
    }
  }

  std::vector<MooseVariable*> _int_vars;
  std::vector<MooseVariable*> _double_vars;
  std::vector<MooseVariable*> _sol_vars; //added by YF

  // Communication maps
  std::map<SPPARKSID, std::vector<unsigned> > _spparks_to_proc; // SPPARKSID to vector of procs that need the value
  std::map<unsigned, std::vector<libMesh::dof_id_type> > _sending_proc_to_elk_id; // Processor to list of ELK ids
  std::map<ELKID, std::vector<unsigned> > _elk_to_proc; // ELKID to vector of procs that need the value
  std::map<unsigned, std::vector<unsigned> > _sending_proc_to_spparks_id; // Processor to list of SPPARKS ids

  unsigned _num_local_elk_nodes;
  unsigned _num_local_spparks_nodes;

  std::map<SPPARKSID, ELKID> _spparks_to_elk;      // Local SPPARKSID to local ELKID
  std::multimap<ELKID, SPPARKSID> _elk_to_spparks; // Local ELKID to local SPPARKSID

  // Maps from variable index to ELK node id to value
  std::map<unsigned, std::map<unsigned, int> > _int_data_for_elk;
  std::map<unsigned, std::map<unsigned, double> > _double_data_for_elk;

  // Maps from variable index to value (vector index is array index)
  std::map<unsigned, std::vector<int> > _int_data_for_spparks;
  std::map<unsigned, std::vector<double> > _double_data_for_spparks;

  Real _last_time;

};

template<>
InputParameters validParams<SPPARKSUserObject>();

#endif
