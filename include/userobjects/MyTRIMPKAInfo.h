#ifndef MYTRIMPKAINFO_H
#define MYTRIMPKAINFO_H

#include "GeneralUserObject.h"

// forward declarations
class MyTRIMPKAInfo;
class MyTRIMRasterizer;

template<>
InputParameters validParams<MyTRIMPKAInfo>();

class MyTRIMPKAInfo : public GeneralUserObject
{
public:
  MyTRIMPKAInfo(const InputParameters & parameters);

  virtual void residualSetup() {}
  virtual void timestepSetup() {}

  virtual void initialize();
  virtual void threadJoin(const UserObject & y);
  virtual void execute();
  virtual void finalize();

protected:
  ///@{ pack/unpack the _npka_per<mass,zaid,energy> data into a structure suitable for parallel communication
  void serialize(std::string & serialized_buffer);
  void deserialize(std::vector<std::string> & serialized_buffers);
  ///@}

  /// prints statistics about created PKAs
  void printPKAStats() const;

  /// the rasterizer object to pull data from
  const MyTRIMRasterizer & _rasterizer;

  /// map of number of PKAs paired with element ID
  // std::map<dof_id_type, unsigned int> _npka_per_element;

  /// map for PKA mass to number of PKAs with this mass
  std::map<unsigned int, unsigned int> _npka_per_mass;

  /// map for PKA ZAID to number of PKAs ZAID
  std::map<unsigned int, unsigned int> _npka_per_zaid;

  /// @{ histogram for pkas by energy
  std::vector<unsigned int> _npka_per_energy;
  unsigned int _nchannels;
  Real _deltaE;
  /// @}
};

#endif
