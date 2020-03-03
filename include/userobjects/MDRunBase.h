/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#pragma once

#include "overlap.hpp"
#include "GeneralUserObject.h"
#include "KDTree.h"
#include "libmesh/bounding_box.h"

class MooseMesh;

/**
 * Base class for molecular dynamics runs in Magpie
 */
class MDRunBase : public GeneralUserObject
{
public:
  static InputParameters validParams();

  MDRunBase(const InputParameters & parameters);

  void initialSetup() override;
  void timestepSetup() override;

  class MDParticles
  {
  public:
    /// particle's position
    std::vector<Point> pos;

    /// the id of particle in the MD calculation
    std::vector<unsigned int> id;

    /// the mesh element the particles are in
    std::vector<unique_id_type> elem_id;

    /// data attached to each particle
    std::vector<std::vector<Real>> properties;

    /// the size of the properties vector
    unsigned int _prop_size;

    /// maps property IDs to position in properties vector
    std::map<unsigned int, unsigned int> _map_props;

    // the index of the radius property in properties vector
    unsigned int _r_index;
  };

  /// access to the MDParticles
  const MDParticles & particles() const { return _md_particles; }

  // check if the stored particles are granular
  bool isGranular() const { return _granular; }

  /// access to the element to particle map
  void elemParticles(unique_id_type elem_id, std::vector<unsigned int> & elem_particles) const;

  /// access the element to granular map
  void granularElementVolumes(unique_id_type elem_id,
                              std::vector<std::pair<unsigned int, Real>> & gran_vol) const;

  /// access to MD particle's properties
  Real particleProperty(unsigned int j, unsigned int prop_id) const;

  /// accessor for md properties that are collected by this UO
  MultiMooseEnum properties() const { return _properties; }

  /// List of quantities to get from MD simulation
  static MultiMooseEnum mdParticleProperties();

  /// helper function to get property index in properties vector
  unsigned int propIndex(unsigned int prop_id) const;

  /// helper function to get property index in properties vector
  unsigned int propIndex(const std::string & prop_name) const;

protected:
  /// call back function to update the particle list
  virtual void updateParticleList() = 0;

  /// updates the KDTree object
  void updateKDTree();

  /// map MDParticles to elements
  void mapMDParticles();

  /// update candidates for
  void updateElementGranularVolumes();

  /// helper function to contruct hexahedron
  OVERLAP::Hexahedron overlapHex(const Elem * elem) const;

  /// helper function to contruct unit hexahedron
  OVERLAP::Hexahedron overlapUnitHex() const;

  /// helper function to contruct tetrahedron
  OVERLAP::Tetrahedron overlapTet(const Elem * elem) const;

  /// helper function to construct unit tetrahedron
  OVERLAP::Tetrahedron overlapUnitTet() const;

  /// Properties that are requested from MD simulation
  MultiMooseEnum _properties;

  /// do the MD particles have extent?
  bool _granular;

  /// The Mesh we're using
  MooseMesh & _mesh;

  /// dimension of the mesh
  const unsigned int _dim;

  /// dimension of the mesh
  const unsigned int _nproc;

  /// the processor bounding box of this processor
  BoundingBox _bbox;

  /// maximum granular radius for parallel bounding boxes
  Real _max_granular_radius;

  /// total number of particles
  unsigned int _n_particles = 0;

  /// number of local particles
  unsigned int _n_local_particles = 0;

  /// stores the location of
  MDParticles _md_particles;

  /// a map from elem unique id to particles in this element
  std::map<unique_id_type, std::vector<unsigned int>> _elem_particles;

  /// a map from element unique id to std::vector of pair(MD particle id, volume of gran. particle in this element)
  std::map<unique_id_type, std::vector<std::pair<unsigned int, Real>>> _elem_granular_volumes;

  /// A KDTree object to handle md_particles
  std::unique_ptr<KDTree> _kd_tree;
};
