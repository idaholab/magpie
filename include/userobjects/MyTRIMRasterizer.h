#ifndef MYTRIMRASTERIZER_H
#define MYTRIMRASTERIZER_H

#include "ElementUserObject.h"

#include <map>
#include <vector>

#include "mytrim/ion.h"

class MyTRIMRasterizer;
class PKAGeneratorBase;

template<>
InputParameters validParams<MyTRIMRasterizer>();

/**
 * This UserObject rasterizes a simulation domain for the MyTRIM library
 */
class MyTRIMRasterizer : public ElementUserObject
{
public:
  MyTRIMRasterizer(const InputParameters & parameters);

  /// determines if a TRIM run is executed during this timestep
  virtual bool executeThisTimestep() const;

  virtual void initialize();
  virtual void execute();
  virtual void threadJoin(const UserObject & y);
  virtual void finalize();

  /// get the concentration array
  const std::vector<Real> & material(const Elem *) const;

  /// get the site volume
  Real siteVolume(const Elem *) const;

  /// get the mass array
  const std::vector<Real> & mass() const { return _trim_mass; }

  /// get the charge array
  const std::vector<Real> & charge() const { return _trim_charge; }

  /// get the PKA list
  const std::vector<MyTRIM_NS::IonBase> & getPKAList() const { return _pka_list; }

  /// get the variable ID of the first coupled variable (to determine the periodicity)
  int periodic() const { return _periodic; }

  /// returns a point with the periodic boundary conditions applied
  Point periodicPoint(const Point &) const;

  /// get the number of elements in the TRIM simulation
  unsigned int nVars() const { return _nvars; }

  /// element averaged data
  struct AveragedData {
    AveragedData(unsigned int nvars = 0) : _elements(nvars, 0.0), _site_volume(0.0) {}

    std::vector<Real> _elements;
    Real _site_volume;
  };

  enum TRIMModuleEnum {
    MYTRIM_CORE = 0,
    MYTRIM_ENERGY_DEPOSITION = 1
  };

  enum Unit {
    ANGSTROM = 0,
    NANOMETER,
    MICROMETER
  };

  TRIMModuleEnum trimModule() const { return _trim_module; }
  Real lengthScale() const { return _length_scale; }

  Real analyticalCutoff() const { return _analytical_cutoff; }

protected:
  ///@{ pack/unpack the _material_map and _pka_list data into a structure suitable for parallel communication
  void serialize(std::string & serialized_buffer);
  void deserialize(std::vector<std::string> & serialized_buffers);
  ///@}

  /// number of coupled variables to map
  const unsigned int _nvars;

  /// dimension of the mesh
  const unsigned int _dim;

  ///@{ Element data
  std::vector<Real> _trim_mass;
  std::vector<Real> _trim_charge;
  ///@}

  /// coupled variable values
  std::vector<const VariableValue *> _var;

  /// lattice site volume material property
  const MaterialProperty<Real> & _site_volume_prop;

  /// @{ PKA generators
  const std::vector<UserObjectName> _pka_generator_names;
  std::vector<const PKAGeneratorBase *> _pka_generators;
  /// @}

  /// @{ material map for the TRIM simulation
  typedef std::map<dof_id_type, AveragedData> MaterialMap;
  MaterialMap _material_map;
  /// @}

  /// variable number to use for minPeriodicDistance calls (i.e. use the periodicity of this variable)
  const int _periodic;

  /// cumulative PKA list
  std::vector<MyTRIM_NS::IonBase> _pka_list;

  /// time for which the next BCMC simulation is responsible (current dt plus skipped steps)
  Real _accumulated_time;

  /// rollback buffer for _accumulated_time if the previous step did not converge
  Real _accumulated_time_old;

  /// timestep interval on which to run BCMC
  const unsigned int _interval;

  /// @{ periodicity data
  Point _min_dim;
  Point _max_dim;
  bool _pbc[LIBMESH_DIM];
  /// @}

  /// the TRIM class to instantiate in the recoil loops
  const TRIMModuleEnum _trim_module;

  ///energy cutoff below which recoils are not followed explicitly but effects are calculated analytically
  const Real _analytical_cutoff;

private:
  bool _execute_this_timestep;

  /// conversion factor from chosen length unit to Angstrom
  Real _length_scale;
};

#endif //MYTRIMRASTERIZER_H
